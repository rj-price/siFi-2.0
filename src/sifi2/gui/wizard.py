"""The "create new database" wizard: a bowtie index built from a FASTA file.

Ported from ``legacy/db_wizard.py``. Three things changed:

* the build runs on :class:`~sifi2.gui.workers.DatabaseWorker` and the wizard
  reacts to its ``completed`` signal, instead of ``start()`` immediately followed
  by ``wait()``, which blocked the GUI for the whole build;
* success is whatever ``bowtie-build`` reported, via
  :func:`sifi2.bowtie.build_database`, rather than the presence of one ``.ebwt``
  file — upstream could not tell a crash from an empty output;
* ``QFileDialog.getOpenFileName`` returns a ``(path, filter)`` tuple in Qt5, so
  upstream's ``.isNull()`` on the result would raise ``AttributeError``.

The index is written into the user's own data directory (the same one the CLI
uses), so no part of this needs elevated privileges.
"""

from __future__ import annotations

import os

from PyQt5 import QtCore, QtWidgets

from .dialogs import show_info_message
from .resources.ui_db_wizard import Ui_wizard
from .workers import DatabaseWorker

__all__ = ["DatabaseWizard"]

#: The page index on which the build runs (``wizardPage2``, the progress page).
BUILD_PAGE = 2


class DatabaseWizard(QtWidgets.QWizard):
    """Choose a FASTA file, name the database, build it."""

    #: Emitted once a database has been built, so the main window can refresh.
    database_created = QtCore.pyqtSignal(str)

    def __init__(self, db_location: str, home_location: str, parent=None):
        super().__init__(parent)
        self.db_location = db_location
        self.home_location = home_location
        self.worker: DatabaseWorker | None = None

        self.ui = Ui_wizard()
        self.ui.setupUi(self)
        self.ui.plainTextEdit.setPlaceholderText(
            "The database is a bowtie index of the sequences you want to search — "
            "a transcriptome, a genome, or the sequences of every organism you care "
            "about off-targets in.\n\nIt is written to:\n" + self.db_location
        )
        self.ui.plainTextEdit.setReadOnly(True)

        self.currentIdChanged.connect(self.page_changed)
        self.ui.commandLinkButton_seq.clicked.connect(self.open_sequence_file)

    # ------------------------------------------------------------------
    def open_sequence_file(self) -> None:
        """Pick the FASTA the index is built from."""
        path, _filter = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open a sequence file", self.home_location, "Sequence (*.fasta *.fa *.fas *.txt)"
        )
        if path and os.path.exists(path):
            self.ui.label_5.setText(path)

    def page_changed(self, page_id: int) -> None:
        if page_id == BUILD_PAGE:
            self.start_build()

    # ------------------------------------------------------------------
    def start_build(self) -> None:
        """Validate the inputs and start the worker."""
        fasta_file = self.ui.label_5.text()
        # Spaces in a database name make a mess of any path it is spliced into.
        db_name = self.ui.lineEdit.text().strip().replace(" ", "_")

        if not os.path.exists(fasta_file):
            self.build_finished(False, "Could not open sequence source file.")
            return
        if not db_name:
            self.build_finished(False, "Please enter a database name.")
            return

        self.button(self.FinishButton).setEnabled(False)
        self.ui.progressBar.setRange(0, 0)  # indeterminate: bowtie-build reports no progress
        self.ui.label_3.setText("Creating bowtie database. Please wait…")
        self.setCursor(QtCore.QCursor(QtCore.Qt.WaitCursor))

        os.makedirs(self.db_location, exist_ok=True)
        self.db_name = db_name
        self.worker = DatabaseWorker(db_name, fasta_file, self.db_location, self)
        self.worker.completed.connect(self.build_finished)
        self.worker.start()

    def build_finished(self, ok: bool, message: str) -> None:
        """Report the outcome and re-enable the wizard."""
        self.ui.progressBar.setRange(0, 100)
        self.ui.progressBar.setValue(100 if ok else 0)
        self.ui.label_3.setText(message if ok else "Database could not be created.")
        self.setCursor(QtCore.QCursor(QtCore.Qt.ArrowCursor))
        self.button(self.FinishButton).setEnabled(True)
        if ok:
            self.database_created.emit(self.db_name)
        else:
            show_info_message(self, message)
