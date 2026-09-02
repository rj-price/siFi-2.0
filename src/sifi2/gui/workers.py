"""Background threads: the pipeline run and the bowtie index build.

Upstream connected its one thread with the old-style
``self.connect(obj, QtCore.SIGNAL("threadDone(QString, QString)"), slot)`` and
then called ``thread.start()`` immediately followed by ``thread.wait()``, which
blocks the GUI thread and makes the thread pointless. Both threads here use
new-style :class:`~PyQt5.QtCore.pyqtSignal` and are genuinely asynchronous: the
window stays responsive and reacts to the finished signal.

The interesting one is :class:`PipelineWorker`. Design mode has to ask the user
which hits are the main targets *in the middle of* the run, and a Qt dialog can
only be shown on the GUI thread — that modal call from inside the pipeline is
exactly the coupling Phase 3 broke. The worker therefore emits
:attr:`PipelineWorker.main_targets_requested` and blocks on a
:class:`threading.Event` until the GUI thread answers with
:meth:`PipelineWorker.provide_main_targets`; answering with ``None`` cancels the
run, which the pipeline reports as :class:`~sifi2.pipeline.PipelineCancelled`.
"""

from __future__ import annotations

import threading

from PyQt5 import QtCore

from ..bowtie import BowtieError, build_database
from ..config import SifiConfig
from ..pipeline import PipelineCancelled, QueryResult, SifiPipeline
from ..rnaplfold import RNAplfoldError

__all__ = ["DatabaseWorker", "PipelineWorker"]

#: Failures a worker reports as a message rather than letting them escape as a
#: traceback into a thread nobody is watching.
RUN_ERRORS = (BowtieError, RNAplfoldError, PipelineCancelled, OSError, ValueError)


class PipelineWorker(QtCore.QThread):
    """Run the pipeline over one or more query records off the GUI thread."""

    #: Human-readable progress, one per query record.
    progressed = QtCore.pyqtSignal(str)
    #: ``list[QueryResult]`` — every record that was analysed.
    completed = QtCore.pyqtSignal(object)
    #: A user-facing failure message; the run produced nothing.
    failed = QtCore.pyqtSignal(str)
    #: ``[(hit_name, hit_count)]`` — answer with :meth:`provide_main_targets`.
    main_targets_requested = QtCore.pyqtSignal(object)

    def __init__(self, config: SifiConfig, queries: list[tuple[str, str]], parent=None):
        super().__init__(parent)
        self.config = config
        self.queries = list(queries)
        self._answered = threading.Event()
        self._selection: list[str] | None = None

    # -- GUI thread ----------------------------------------------------
    def provide_main_targets(self, selection: list[str] | None) -> None:
        """Answer :attr:`main_targets_requested`. ``None`` cancels the run."""
        self._selection = selection
        self._answered.set()

    def cancel(self) -> None:
        """Unblock a run waiting on a target selection, cancelling it."""
        self.provide_main_targets(None)

    # -- worker thread -------------------------------------------------
    def _select_main_targets(self, choices: list[tuple[str, int]]) -> list[str] | None:
        self._selection = None
        self._answered.clear()
        self.main_targets_requested.emit(choices)
        self._answered.wait()
        return self._selection

    def run(self) -> None:  # noqa: D102 - QThread's entry point
        selector = self._select_main_targets if self.config.is_design else None
        pipeline = SifiPipeline(self.config, main_target_selector=selector)
        results: list[QueryResult] = []
        try:
            for index, (name, sequence) in enumerate(self.queries, start=1):
                self.progressed.emit(f"Analysing {name} ({index} of {len(self.queries)})…")
                results.append(pipeline.run_query(name, sequence))
        except RUN_ERRORS as error:
            self.failed.emit(str(error))
            return
        self.completed.emit(results)


class DatabaseWorker(QtCore.QThread):
    """Build a bowtie index off the GUI thread, for the database wizard."""

    #: ``(ok, message)``.
    completed = QtCore.pyqtSignal(bool, str)

    def __init__(self, db_name: str, fasta_file: str, db_location: str, parent=None):
        super().__init__(parent)
        self.db_name = db_name
        self.fasta_file = fasta_file
        self.db_location = db_location

    def run(self) -> None:  # noqa: D102 - QThread's entry point
        try:
            build_database(self.db_name, self.fasta_file, self.db_location)
        except (BowtieError, OSError) as error:
            self.completed.emit(False, str(error))
            return
        self.completed.emit(True, "Database successfully created!")
