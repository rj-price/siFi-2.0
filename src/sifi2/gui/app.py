"""The main window, and the ``sifi-gui`` entry point.

Ported from ``legacy/main.py``. The window itself is the same one — the layout
comes from the same ``sifi2015.ui``, regenerated with ``pyuic5`` — but what sits
behind it is the ported core: the window builds a
:class:`~sifi2.config.SifiConfig` from its widgets, runs
:class:`~sifi2.pipeline.SifiPipeline` on a worker thread, and hands the results
to :class:`~sifi2.gui.results.ResultsWindow`.

**Nothing here needs administrator rights.** Upstream derived its application
directory as ``QDesktopServices.storageLocation(DataLocation).split('Local')[0] +
'/Local/'``, created ``Bowtie/``, ``RNAplfold/``, ``Images/`` and ``Databases/``
under it, and copied the bundled Windows ``.exe``\\ s into place on every start-up
(``general_helpers.copying_files``) — which, from a ``Program Files`` install, is
a privileged write. This version:

* resolves ``bowtie``, ``bowtie-build`` and ``RNAplfold`` on ``PATH`` with
  :func:`shutil.which` and copies nothing, ever;
* keeps databases in :func:`sifi2.cli.default_db_location` — a ``platformdirs``
  user data directory, the same one the CLI uses, so the two share databases;
* keeps preferences in ``QSettings`` rather than the ``last_db.txt`` and
  ``logging_sifi.txt`` files upstream created beside the binaries;
* holds its images in the compiled Qt resource, so there is no image directory
  to populate.

Two further changes to upstream behaviour:

* ``QDesktopServices.storageLocation`` was **removed in Qt5**;
  ``QStandardPaths.writableLocation`` replaces it (five call sites upstream).
* Multi-record FASTA input is accepted and every record is analysed. Upstream
  refused it — ``"Please enter only one sequence or use the batch mode."`` — and
  the batch mode it pointed at was never written.
"""

from __future__ import annotations

import sys
import webbrowser
from io import StringIO

from PyQt5 import QtCore, QtGui, QtWidgets

from ..bowtie import BowtieError, all_databases, database_exists, delete_databases, resolve_binary
from ..cli import default_db_location
from ..config import Mode, SifiConfig, qt_flag
from ..pipeline import QueryResult
from .dialogs import ask_main_targets, ask_mode, show_info_message
from .resources.ui_sifi2015 import Ui_MainWindow
from .results import ResultsWindow
from .wizard import DatabaseWizard
from .workers import PipelineWorker

__all__ = ["MainWindow", "main"]

DOCUMENTATION_URL = "http://labtools.ipk-gatersleben.de/si-Fi/Quick_help.pdf"
PROJECT_URL = "http://www.snowformatics.com/"

#: The efficiency widgets, disabled together whenever mismatches > 0.
EFFICIENCY_WIDGET_NAMES = [
    "checkBox_end",
    "doubleSpinBox_end",
    "checkBox_access",
    "doubleSpinBox_access",
    "spinBox_xmer",
    "label_5",
    "checkBox_strand",
    "checkBox_terminal",
]


def home_location() -> str:
    """The user's home directory, for file dialogs.

    ``QDesktopServices.storageLocation`` was removed in Qt5.
    """
    return QtCore.QStandardPaths.writableLocation(QtCore.QStandardPaths.HomeLocation)


class AboutWindow(QtWidgets.QWidget):
    """The about box: one image, clickable through to the project page."""

    def __init__(self, parent=None):
        super().__init__(parent, QtCore.Qt.Window)
        self.setWindowTitle("About siFi 2.0")
        label = QtWidgets.QLabel(self)
        label.setPixmap(QtGui.QPixmap(":/Images/Images/about.png"))
        label.mousePressEvent = lambda event: webbrowser.open(PROJECT_URL)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(label)


class MainWindow(QtWidgets.QMainWindow):
    """siFi's main window: sequence in, settings, database, run."""

    def __init__(self, mode: Mode = Mode.DESIGN, db_location: str | None = None, parent=None):
        super().__init__(parent)
        self.db_location = db_location or default_db_location()
        self.settings = QtCore.QSettings("siFi", "sifi2")
        self.worker: PipelineWorker | None = None
        self.result_windows: list[ResultsWindow] = []

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowIcon(QtGui.QIcon(":/Images/Images/siFi21_icon64.png"))
        # The two tabs are fixtures of the window, not documents to close.
        for index in (0, 1):
            self.ui.tabWidget.tabBar().setTabButton(index, QtWidgets.QTabBar.RightSide, None)

        self.create_connects()
        self.set_mode(mode)
        self.refresh_databases()
        self.check_external_binaries()
        self.ui.plainTextEdit_seq.setFocus()

    # ------------------------------------------------------------------
    # Wiring
    # ------------------------------------------------------------------
    @property
    def efficiency_widgets(self) -> list:
        return [getattr(self.ui, name) for name in EFFICIENCY_WIDGET_NAMES]

    def create_connects(self) -> None:
        """Connect the widgets. All new-style signals; upstream used ``SIGNAL()``."""
        self.ui.commandLinkButton_remove.clicked.connect(self.delete_selected_databases)
        self.ui.commandLinkButton_add.clicked.connect(self.start_db_wizard)
        self.ui.actionRNAi_design.triggered.connect(lambda: self.set_mode(Mode.DESIGN))
        self.ui.actionOff_target_prediction_2.triggered.connect(lambda: self.set_mode(Mode.OFFTARGET))
        self.ui.pushButton_run.clicked.connect(self.start_pipeline)
        self.ui.pushButton_default.clicked.connect(self.apply_default_settings)
        self.ui.actionDocumentation.triggered.connect(lambda: webbrowser.open(DOCUMENTATION_URL))
        self.ui.actionAbout.triggered.connect(self.show_about)
        self.ui.pushButton_openDB.clicked.connect(self.open_sequence_file)
        self.ui.comboBox_miss.currentIndexChanged.connect(self.mismatches_changed)

    def check_external_binaries(self) -> None:
        """Warn once, in the status bar, if bowtie or RNAplfold is not on PATH."""
        missing = []
        for name in ("bowtie", "bowtie-build", "RNAplfold"):
            try:
                resolve_binary(name)
            except BowtieError:
                missing.append(name)
        if missing:
            self.statusBar().showMessage(f"Not found on PATH: {', '.join(missing)} — runs will fail until installed.")
        else:
            self.statusBar().showMessage(f"Databases: {self.db_location}")

    # ------------------------------------------------------------------
    # Modes and defaults
    # ------------------------------------------------------------------
    def set_mode(self, mode: Mode) -> None:
        """Switch mode and apply that mode's default settings."""
        self.mode = mode
        for widget in self.efficiency_widgets:
            widget.setEnabled(True)
        design = mode is Mode.DESIGN
        self.ui.checkBox_end.setChecked(design)
        self.ui.checkBox_access.setChecked(design)
        self.ui.checkBox_strand.setChecked(True)
        self.ui.checkBox_terminal.setChecked(True)
        self.ui.spinBox_size.setValue(21)
        self.ui.comboBox_miss.setCurrentIndex(0)
        self.ui.doubleSpinBox_end.setValue(1)
        self.ui.doubleSpinBox_access.setValue(0.1)
        self.ui.spinBox_xmer.setValue(8)
        self.ui.tabWidget.setTabText(0, "RNAi design" if design else "Off-target prediction")
        self.ui.groupBox_1.setTitle("" if design else "Paste RNAi trigger sequence or import a file")
        self.ui.plainTextEdit_seq.setFocus()

    def apply_default_settings(self) -> None:
        self.set_mode(self.mode)

    def mismatches_changed(self, index: int) -> None:
        """Efficiency prediction is only defined for exact siRNA matches.

        The GUI has no "predict efficiency" checkbox: upstream derived
        ``no_efficience`` from whether these widgets were *enabled*, and disabled
        them as soon as mismatches > 0. The CLI's ``--efficiency`` flag is the
        same rule made explicit.
        """
        if index > 0:
            for widget in self.efficiency_widgets:
                widget.setEnabled(False)
        else:
            self.apply_default_settings()

    # ------------------------------------------------------------------
    # Input
    # ------------------------------------------------------------------
    def open_sequence_file(self) -> None:
        """Load a FASTA or plain-text sequence file into the text box."""
        path, _filter = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open a sequence file", home_location(), "Sequence (*.fasta *.fa *.fas *.txt)"
        )
        if not path:
            return
        try:
            with open(path) as handle:
                text = handle.read()
        except OSError as error:
            show_info_message(self, f"Could not open {path}: {error}")
            return
        self.ui.label_1.setText(path)
        self.ui.plainTextEdit_seq.setPlainText(text)

    def read_queries(self) -> list[tuple[str, str]] | None:
        """The text box as ``[(name, sequence)]``, or ``None`` if it is not valid.

        Plain text with no header is wrapped as ``my_sequence``, exactly as
        upstream did; FASTA is parsed as FASTA, and every record is kept.
        """
        from Bio import SeqIO

        from .. import analysis

        text = self.ui.plainTextEdit_seq.toPlainText().strip()
        if not text:
            show_info_message(self, "Please enter a nucleic acid sequence!")
            return None

        if not text.startswith(">"):
            if not analysis.validate_seq(text):
                show_info_message(self, "Please enter a valid nucleic acid sequence!")
                return None
            text = ">my_sequence\n" + text
        elif not analysis.validate_fasta_seq(text):
            show_info_message(self, "Please enter a valid nucleic acid sequence in FASTA format!")
            return None

        queries = [(record.id, str(record.seq)) for record in SeqIO.parse(StringIO(text), "fasta")]
        if not queries:
            show_info_message(self, "No FASTA records found in the sequence box.")
            return None

        names = [name for name, _sequence in queries]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            # Every record gets its own result; two records with one name cannot
            # be told apart in the results window or in an export.
            show_info_message(self, f"Duplicate sequence names: {', '.join(duplicates)}")
            return None
        return queries

    def build_config(self) -> SifiConfig:
        """Read the widgets into a :class:`~sifi2.config.SifiConfig`."""
        ui = self.ui
        # Qt hands over 0/2 for a checkbox, never 0/1.
        return SifiConfig(
            bowtie_db=ui.comboBox_db.currentText(),
            db_location=self.db_location,
            mode=self.mode,
            sirna_size=int(ui.spinBox_size.value()),
            mismatches=int(ui.comboBox_miss.currentText()),
            strand_check=qt_flag(ui.checkBox_strand.checkState()),
            end_check=qt_flag(ui.checkBox_end.checkState()),
            accessibility_check=qt_flag(ui.checkBox_access.checkState()),
            terminal_check=qt_flag(ui.checkBox_terminal.checkState()),
            no_efficience=not all(widget.isEnabled() for widget in self.efficiency_widgets),
            end_stability_threshold=float(ui.doubleSpinBox_end.value()),
            accessibility_threshold=float(ui.doubleSpinBox_access.value()),
            accessibility_window=int(ui.spinBox_xmer.value()),
        )

    # ------------------------------------------------------------------
    # Running
    # ------------------------------------------------------------------
    def start_pipeline(self) -> None:
        """Validate everything, then run the pipeline on a worker thread."""
        if self.worker is not None and self.worker.isRunning():
            show_info_message(self, "A run is already in progress.")
            return

        config = self.build_config()
        if not config.bowtie_db:
            show_info_message(self, "Please create or select a database first.")
            return
        if not database_exists(config.bowtie_db, config.db_location):
            show_info_message(self, f"Database {config.bowtie_db!r} is missing files. Please create it again.")
            return

        queries = self.read_queries()
        if queries is None:
            return

        self.settings.setValue("last_db", config.bowtie_db)
        self.ui.pushButton_run.setEnabled(False)
        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))

        self.worker = PipelineWorker(config, queries, self)
        self.worker.progressed.connect(self.statusBar().showMessage)
        self.worker.main_targets_requested.connect(self.choose_main_targets)
        self.worker.completed.connect(self.run_finished)
        self.worker.failed.connect(self.run_failed)
        self.worker.start()

    def choose_main_targets(self, choices: list[tuple[str, int]]) -> None:
        """Answer the worker's mid-run question. Runs on the GUI thread.

        The pipeline is blocked inside ``get_main_targets`` until this returns —
        which is why it must be answered even when the dialog is cancelled, and
        why cancelling means cancelling the run rather than "no main targets".
        """
        self.worker.provide_main_targets(ask_main_targets(choices, self))

    def run_finished(self, results: list[QueryResult]) -> None:
        """Show the results, and report any query that produced none."""
        self._run_ended()
        for result in results:
            if result.message:
                show_info_message(self, f"{result.query_name}: {result.message}")

        plotted = [result for result in results if result.has_hits]
        if not plotted:
            self.statusBar().showMessage("No results were produced.")
            return

        from .. import plots

        # Global rc state, so the GUI applies it once at the top rather than
        # having the plotting functions mutate it behind the caller's back.
        plots.apply_style()
        window = ResultsWindow(plotted, self.ui.spinBox_size.value(), self.mode is Mode.DESIGN, None)
        window.closed.connect(lambda: self.result_windows.remove(window) if window in self.result_windows else None)
        self.result_windows.append(window)
        window.show()
        self.statusBar().showMessage(f"Finished: {len(plotted)} of {len(results)} query record(s) produced results.")

    def run_failed(self, message: str) -> None:
        self._run_ended()
        self.statusBar().showMessage("Run failed.")
        show_info_message(self, message)

    def _run_ended(self) -> None:
        QtWidgets.QApplication.restoreOverrideCursor()
        self.ui.pushButton_run.setEnabled(True)

    # ------------------------------------------------------------------
    # Databases
    # ------------------------------------------------------------------
    def start_db_wizard(self) -> None:
        self.wizard = DatabaseWizard(self.db_location, home_location(), self)
        self.wizard.database_created.connect(self.refresh_databases)
        self.wizard.show()

    def refresh_databases(self, select: str | None = None) -> None:
        """Repopulate the database combo and the management table."""
        import os

        databases = (
            all_databases(self.db_location)
            if os.path.isdir(self.db_location)
            else {"Database name": [], "Database size (MB)": [], "Created": []}
        )
        names = list(databases["Database name"])

        # Show the database used last (or just created) first.
        preferred = select or self.settings.value("last_db", "")
        if preferred in names:
            names.insert(0, names.pop(names.index(preferred)))
        self.ui.comboBox_db.clear()
        self.ui.comboBox_db.addItems(names)

        columns = ["Database name", "Database size (MB)", "Created"]
        table = self.ui.tableWidget
        table.setSortingEnabled(False)
        table.clear()
        table.setRowCount(len(databases["Database name"]))
        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels(columns)
        for row in range(len(databases["Database name"])):
            for column, key in enumerate(columns):
                item = QtWidgets.QTableWidgetItem(str(databases[key][row]))
                if column == 0:
                    item.setCheckState(QtCore.Qt.Unchecked)
                table.setItem(row, column, item)
        table.setSortingEnabled(True)

    def checked_databases(self) -> list[str]:
        """The database names ticked in the management table."""
        table = self.ui.tableWidget
        checked = []
        for row in range(table.rowCount()):
            item = table.item(row, 0)
            if item is not None and item.checkState() == QtCore.Qt.Checked:
                checked.append(item.text())
        return checked

    def delete_selected_databases(self) -> None:
        """Delete the ticked databases, after confirming."""
        to_delete = self.checked_databases()
        if not to_delete:
            show_info_message(self, "Please select a database!")
            return
        confirm = QtWidgets.QMessageBox.question(
            self,
            "Delete database?",
            "Are you sure that you want to delete the database(s)?\n\n" + "\n".join(to_delete),
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.No,
        )
        if confirm != QtWidgets.QMessageBox.Yes:
            return
        deleted = delete_databases(to_delete, self.db_location)
        self.refresh_databases()
        if len(deleted) == len(to_delete):
            show_info_message(self, "Database successfully deleted!")
        else:
            missed = sorted(set(to_delete) - set(deleted))
            show_info_message(self, f"Could not delete: {', '.join(missed)}")

    # ------------------------------------------------------------------
    def show_about(self) -> None:
        self.about_window = AboutWindow(self)
        self.about_window.show()

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt's spelling
        """Let a running pipeline out of its wait before the window goes away."""
        if self.worker is not None and self.worker.isRunning():
            self.worker.cancel()
            self.worker.wait(2000)
        for window in list(self.result_windows):
            window.close()
        event.accept()


def main(argv: list[str] | None = None, db_location: str | None = None) -> int:
    """Start the GUI. Returns the Qt exit code.

    ``argv`` is handed to Qt, so the usual platform arguments work —
    ``-platform offscreen`` in particular, which is how the tests run this on a
    node with no display.
    """
    # Reuse an existing instance if one is already up: Qt allows only one per
    # process, and the tests drive this entry point from a live application.
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(
        sys.argv if argv is None else [sys.argv[0], *argv]
    )
    app.setApplicationName("siFi 2.0")
    # Upstream forced the "gtk" style, which does not exist in Qt5; the platform
    # default is the right choice anyway.
    app.setPalette(app.style().standardPalette())

    mode = ask_mode()
    if mode is None:
        # Cancelling the mode question used to return 0 from `exec_()`, and 0 was
        # design mode, so closing this dialog silently started a design run
        # (PLAN.md Phase 6, handed to Phase 7).
        return 0

    window = MainWindow(mode, db_location=db_location)
    available_height = app.primaryScreen().availableGeometry().height()
    if available_height < 830:
        window.resize(640, available_height - 50)
    window.show()
    return app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
