"""Phase 7: the PyQt5 GUI.

The GUI is tested through the offscreen Qt platform plugin, so this file runs on
a headless cluster node with no display. Every test skips if PyQt5 is absent —
it is an optional dependency (``pip install sifi2[gui]``), and the core must keep
working without it.

What is worth testing here is not the layout, which comes from a generated
``.ui`` file, but the seams: the widgets-to-:class:`~sifi2.config.SifiConfig`
mapping, the sequence-box parsing, the database table, the mid-run main-target
handshake between the worker thread and the GUI thread, and the four upstream
defects the port fixes in Qt code (the mode dialog's cancel, the doubled export
filename, the ``getSaveFileName`` tuple, and the refusal of multi-record FASTA).
"""

import json
import os
import shutil
import time

import pytest
from conftest import BASELINE, DATA

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt5", reason="the GUI is an optional extra; install sifi2[gui]")

from PyQt5 import QtCore, QtWidgets  # noqa: E402

from sifi2.bowtie import INDEX_EXTENSIONS  # noqa: E402
from sifi2.config import Mode  # noqa: E402
from sifi2.gui import app as gui_app  # noqa: E402
from sifi2.gui import dialogs, results, workers  # noqa: E402
from sifi2.pipeline import QueryResult  # noqa: E402

needs_binaries = pytest.mark.skipif(
    shutil.which("bowtie") is None or shutil.which("bowtie-build") is None or shutil.which("RNAplfold") is None,
    reason="needs bowtie, bowtie-build and RNAplfold on PATH",
)


@pytest.fixture(scope="session")
def qapp():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield application


@pytest.fixture
def quiet(monkeypatch):
    """Collect the messages that would have been modal information boxes."""
    messages = []
    for module in (gui_app, dialogs, results):
        monkeypatch.setattr(
            module,
            "show_info_message",
            lambda parent, message, title="Information": messages.append(message),
            raising=False,
        )
    return messages


@pytest.fixture
def window(qapp, tmp_path, quiet):
    main_window = gui_app.MainWindow(Mode.DESIGN, db_location=str(tmp_path / "dbs"))
    yield main_window
    main_window.close()


def fake_index(db_location, name):
    """The six files ``bowtie-build`` writes, with no bowtie involved."""
    db_location.mkdir(parents=True, exist_ok=True)
    for extension in INDEX_EXTENSIONS:
        (db_location / (name + extension)).write_bytes(b"x" * 2_000_000)


def baseline_result(name="offtarget_TUB3_fragment.json"):
    records = json.loads((BASELINE / name).read_text())
    sequence = "".join(line.strip() for line in (DATA / "query.fasta").read_text().splitlines()[1:])
    from sifi2 import analysis

    return QueryResult(
        query_name="TUB3_fragment",
        query_sequence=sequence,
        records=records,
        table_data=analysis.get_table_data(records),
        main_targets=[row[0] for row in analysis.get_table_data(records)],
    )


# ----------------------------------------------------------------------
# The mode question
# ----------------------------------------------------------------------
def test_cancelling_the_mode_dialog_does_not_silently_choose_design(qapp):
    """``PLAN.md`` Phase 6's last defect, handed to Phase 7.

    Upstream did ``self.mode = msg_box.exec_()``, and closing a ``QMessageBox``
    with the window's X returns ``0`` — which was design mode. Cancelling the
    question therefore started a design run.
    """
    dialog = dialogs.ModeDialog()
    dialog.reject()
    assert dialog.selected_mode() is None

    closed = dialogs.ModeDialog()
    closed.close()
    assert closed.selected_mode() is None


@pytest.mark.parametrize("mode", [Mode.DESIGN, Mode.OFFTARGET])
def test_the_mode_dialog_returns_the_chosen_mode(qapp, mode):
    dialog = dialogs.ModeDialog()
    dialog._choose(mode)
    assert dialog.selected_mode() is mode
    assert dialog.result() == QtWidgets.QDialog.Accepted


# ----------------------------------------------------------------------
# Main-target selection
# ----------------------------------------------------------------------
def test_target_dialog_returns_the_ticked_names(qapp):
    choices = [("TUB4", 500), ("TUB3", 400), ("TUB2", 12)]
    dialog = dialogs.TargetSelectionDialog(choices)
    dialog.tree.topLevelItem(0).setCheckState(0, QtCore.Qt.Checked)
    dialog.tree.topLevelItem(2).setCheckState(0, QtCore.Qt.Checked)
    assert dialog.selected_targets() == ["TUB4", "TUB2"]

    dialog.select_all()
    assert dialog.selected_targets() == ["TUB4", "TUB3", "TUB2"]


def test_target_dialog_handles_no_targets(qapp):
    """Design mode with no hits still scores the siRNAs (Phase 6 defect 2), so
    an empty choice list is an outcome, not an error."""
    dialog = dialogs.TargetSelectionDialog([])
    assert dialog.tree is None
    assert dialog.selected_targets() == []


# ----------------------------------------------------------------------
# Widgets to SifiConfig
# ----------------------------------------------------------------------
def test_design_defaults_match_the_cli(window):
    config = window.build_config()
    assert config.mode is Mode.DESIGN
    assert (config.sirna_size, config.mismatches) == (21, 0)
    assert config.end_check and config.accessibility_check
    assert config.strand_check and config.terminal_check
    assert (config.end_stability_threshold, config.accessibility_threshold, config.accessibility_window) == (
        1.0,
        0.1,
        8,
    )
    assert config.no_efficience is False


def test_offtarget_defaults_turn_two_rules_off(window):
    window.set_mode(Mode.OFFTARGET)
    config = window.build_config()
    assert config.mode is Mode.OFFTARGET
    assert not config.end_check and not config.accessibility_check
    assert config.strand_check and config.terminal_check


def test_checkbox_states_reach_the_config_as_bools(window):
    """Qt hands over 0/2, and ``SifiConfig`` holds real booleans."""
    window.ui.checkBox_strand.setCheckState(QtCore.Qt.Unchecked)
    window.ui.checkBox_terminal.setCheckState(QtCore.Qt.Checked)
    config = window.build_config()
    assert config.strand_check is False
    assert config.terminal_check is True


def test_mismatches_disable_efficiency_prediction(window):
    """The GUI never had a "predict efficiency" checkbox: ``no_efficience`` was
    derived from whether the efficiency widgets were enabled, and mismatches > 0
    disabled all four. The CLI's ``--efficiency`` is the same rule, explicit."""
    window.ui.comboBox_miss.setCurrentIndex(2)
    assert not any(widget.isEnabled() for widget in window.efficiency_widgets)
    config = window.build_config()
    assert config.mismatches == 2
    assert config.no_efficience is True

    window.ui.comboBox_miss.setCurrentIndex(0)
    assert all(widget.isEnabled() for widget in window.efficiency_widgets)
    assert window.build_config().no_efficience is False


# ----------------------------------------------------------------------
# The sequence box
# ----------------------------------------------------------------------
def test_plain_sequence_is_wrapped_as_a_fasta_record(window):
    window.ui.plainTextEdit_seq.setPlainText("ACGTACGTAC\nGTACGTACGT")
    assert window.read_queries() == [("my_sequence", "ACGTACGTACGTACGTACGT")]


def test_multi_record_fasta_is_accepted(window):
    """Upstream refused this — "Please enter only one sequence or use the batch
    mode" — and the batch mode it pointed at was never written."""
    window.ui.plainTextEdit_seq.setPlainText(">one\nACGTACGT\n>two\nTTTTGGGG\n>three\nCCCCAAAA")
    assert window.read_queries() == [("one", "ACGTACGT"), ("two", "TTTTGGGG"), ("three", "CCCCAAAA")]


def test_invalid_sequences_are_rejected(window, quiet):
    window.ui.plainTextEdit_seq.setPlainText("this is not a sequence!")
    assert window.read_queries() is None
    assert "valid nucleic acid sequence" in quiet[-1]


def test_empty_sequence_box_is_rejected(window, quiet):
    window.ui.plainTextEdit_seq.setPlainText("   \n ")
    assert window.read_queries() is None
    assert quiet


def test_duplicate_record_names_are_rejected(window, quiet):
    """Two records with one name cannot be told apart in the results window or
    in an export — the same check the CLI makes on output filenames."""
    window.ui.plainTextEdit_seq.setPlainText(">dup\nACGTACGT\n>dup\nTTTTGGGG")
    assert window.read_queries() is None
    assert "Duplicate sequence names: dup" in quiet[-1]


# ----------------------------------------------------------------------
# Databases
# ----------------------------------------------------------------------
def test_database_table_lists_indices_and_reports_selection(window, tmp_path):
    fake_index(tmp_path / "dbs", "alpha")
    fake_index(tmp_path / "dbs", "beta")
    window.refresh_databases()

    assert window.ui.tableWidget.rowCount() == 2
    assert sorted(window.ui.comboBox_db.itemText(i) for i in range(window.ui.comboBox_db.count())) == ["alpha", "beta"]
    assert window.checked_databases() == []

    rows = {window.ui.tableWidget.item(row, 0).text(): row for row in range(window.ui.tableWidget.rowCount())}
    window.ui.tableWidget.item(rows["beta"], 0).setCheckState(QtCore.Qt.Checked)
    assert window.checked_databases() == ["beta"]


def test_the_last_used_database_is_listed_first(window, tmp_path):
    fake_index(tmp_path / "dbs", "alpha")
    fake_index(tmp_path / "dbs", "zulu")
    window.settings.setValue("last_db", "zulu")
    window.refresh_databases()
    assert window.ui.comboBox_db.itemText(0) == "zulu"


def test_deleting_databases_asks_first(window, tmp_path, monkeypatch, quiet):
    fake_index(tmp_path / "dbs", "alpha")
    window.refresh_databases()
    window.ui.tableWidget.item(0, 0).setCheckState(QtCore.Qt.Checked)

    monkeypatch.setattr(QtWidgets.QMessageBox, "question", staticmethod(lambda *a, **k: QtWidgets.QMessageBox.No))
    window.delete_selected_databases()
    assert (tmp_path / "dbs" / "alpha.1.ebwt").exists()

    monkeypatch.setattr(QtWidgets.QMessageBox, "question", staticmethod(lambda *a, **k: QtWidgets.QMessageBox.Yes))
    window.ui.tableWidget.item(0, 0).setCheckState(QtCore.Qt.Checked)
    window.delete_selected_databases()
    assert not (tmp_path / "dbs" / "alpha.1.ebwt").exists()
    assert window.ui.tableWidget.rowCount() == 0


def test_running_without_a_database_says_so(window, quiet):
    window.ui.plainTextEdit_seq.setPlainText("ACGTACGTACGTACGTACGTACGT")
    window.start_pipeline()
    assert "database" in quiet[-1].lower()
    assert window.worker is None


# ----------------------------------------------------------------------
# The results window
# ----------------------------------------------------------------------
@pytest.fixture
def results_window(qapp):
    from sifi2 import plots

    plots.apply_style()
    window = results.ResultsWindow([baseline_result()], sirna_size=21, is_design=False)
    yield window
    window.close()


def test_results_window_draws_a_canvas_per_query(qapp):
    one = baseline_result()
    two = QueryResult(
        query_name="second",
        query_sequence=one.query_sequence,
        records=one.records,
        table_data=one.table_data,
        main_targets=one.main_targets,
    )
    window = results.ResultsWindow([one, two], sirna_size=21, is_design=False)
    try:
        assert window.selector.count() == 2
        assert window._chooser_widget.isVisibleTo(window)
        first_canvas = window._canvas
        window.selector.setCurrentIndex(1)
        assert window.current.query_name == "second"
        assert window._canvas is not first_canvas
    finally:
        window.close()


def test_a_single_query_hides_the_chooser(results_window):
    assert not results_window._chooser_widget.isVisibleTo(results_window)
    assert results_window._canvas is not None


def test_export_extension_is_appended_not_doubled():
    """Upstream wrote ``filename += filename + '.png'`` in four places, which
    saves "results" as "resultsresults.png"."""
    assert results._with_suffix("/tmp/plot", ".png") == "/tmp/plot.png"
    assert results._with_suffix("/tmp/plot.png", ".png") == "/tmp/plot.png"
    assert results._with_suffix("/tmp/plot.PNG", ".png") == "/tmp/plot.PNG"


def test_exports_write_real_files(results_window, tmp_path, monkeypatch, quiet):
    written = []

    def save_path(caption, suffix):
        path = str(tmp_path / ("export" + suffix))
        written.append(path)
        return path

    monkeypatch.setattr(results_window, "_save_path", save_path)

    results_window.export_image()
    results_window.export_table()
    assert (tmp_path / "export.png").stat().st_size > 0
    table = (tmp_path / "export.csv").read_text().splitlines()
    assert table[0] == "Targets;Total siRNA hits;Efficient siRNA hits"
    assert len(table) == 1 + len(results_window.current.table_data)


def test_genbank_export_is_reachable_in_design_mode(qapp, tmp_path, monkeypatch, quiet):
    """It was dead code upstream: ``imageviewer`` was only ever constructed in a
    commented-out block, so nothing could call it. ``PLAN.md`` Phase 7 asked
    whether to keep it — it is kept, and now reachable."""
    window = results.ResultsWindow([baseline_result("design_TUB3_fragment.json")], sirna_size=21, is_design=True)
    try:
        actions = [
            action.text() for menu in window.menuBar().findChildren(QtWidgets.QMenu) for action in menu.actions()
        ]
        assert "GenBank file" in actions
        monkeypatch.setattr(window, "_save_path", lambda caption, suffix: str(tmp_path / ("export" + suffix)))
        window.export_gbk()
        text = (tmp_path / "export.gbk").read_text()
        assert text.startswith("LOCUS")
        assert "MT_" in text or "OT_" in text
    finally:
        window.close()


def test_offtarget_mode_has_no_genbank_export(results_window):
    actions = [
        action.text() for menu in results_window.menuBar().findChildren(QtWidgets.QMenu) for action in menu.actions()
    ]
    assert "GenBank file" not in actions
    assert "Image file" in actions and "Table (CSV)" in actions


# ----------------------------------------------------------------------
# The worker thread
# ----------------------------------------------------------------------
def test_the_worker_blocks_until_the_main_targets_come_back(qapp):
    """The pipeline asks for main targets from the worker thread, and only the
    GUI thread may show the dialog — so the worker emits and waits."""
    from sifi2.config import SifiConfig

    worker = workers.PipelineWorker(SifiConfig(mode=Mode.DESIGN), [])
    answers = []

    def answer(choices):
        answers.append(choices)
        worker.provide_main_targets(["TUB3"])

    # Emitted and handled on this thread, so the handler runs inside `emit()` and
    # the worker's wait is already satisfied by the time it starts waiting.
    worker.main_targets_requested.connect(answer)
    assert worker._select_main_targets([("TUB3", 10)]) == ["TUB3"]
    assert answers == [[("TUB3", 10)]]

    worker.main_targets_requested.disconnect(answer)
    worker.main_targets_requested.connect(lambda choices: worker.cancel())
    assert worker._select_main_targets([("TUB3", 10)]) is None


@needs_binaries
def test_a_full_gui_run_reproduces_the_cli_baseline(qapp, window, tmp_path, monkeypatch, quiet):
    """The GUI, the CLI and the end-to-end baseline must agree exactly: they are
    the same pipeline, and the GUI only fills in the config."""
    from sifi2 import bowtie

    db_location = tmp_path / "dbs"
    db_location.mkdir(parents=True, exist_ok=True)
    bowtie.build_database("testdb", str(DATA / "reference.fasta"), str(db_location))
    window.refresh_databases()
    window.ui.comboBox_db.setCurrentIndex(window.ui.comboBox_db.findText("testdb"))
    window.ui.plainTextEdit_seq.setPlainText((DATA / "query.fasta").read_text())

    # Stand in for the modal dialog: every hit is a main target, as
    # `sifi design --all-targets-main` does.
    monkeypatch.setattr(gui_app, "ask_main_targets", lambda choices, parent=None: [name for name, _ in choices])

    window.start_pipeline()
    # The worker's signals are queued to this thread, and the main-target
    # question blocks the run until this loop delivers it — a bare `wait()` here
    # would deadlock exactly as upstream's `start()`/`wait()` pair did.
    deadline = time.monotonic() + 600
    while window.worker.isRunning() and time.monotonic() < deadline:
        qapp.processEvents()
        window.worker.wait(50)
    assert not window.worker.isRunning(), "the pipeline run timed out"
    qapp.processEvents()

    assert window.result_windows, "a results window should have opened"
    results_window = window.result_windows[0]
    records = results_window.results[0].records
    assert records == json.loads((BASELINE / "design_TUB3_fragment.json").read_text())
    results_window.close()


# ----------------------------------------------------------------------
# Packaging
# ----------------------------------------------------------------------
def test_generated_ui_imports_its_resource_module_relatively():
    """``pyuic5`` without ``--from-imports`` emits a bare ``import
    sifi_2015_rc``, which under Python 3 resolves against ``sys.path`` — and
    upstream had a 0-byte ``sifi_2015_rc.py`` at its repo root that would satisfy
    it, register nothing, and render every icon blank (``PLAN.md`` Phase 3)."""
    from pathlib import Path

    from sifi2.gui import resources

    for name in ("ui_sifi2015.py", "ui_db_wizard.py"):
        text = (Path(resources.__file__).parent / name).read_text()
        assert "from . import sifi_2015_rc" in text
        assert "\nimport sifi_2015_rc" not in text


def test_the_resource_images_are_registered(qapp):
    from PyQt5 import QtGui

    for path in (":/Images/Images/about.png", ":/Images/Images/siFi21_icon64.png", ":/Images/Images/add.png"):
        assert not QtGui.QPixmap(path).isNull(), path


def test_the_entry_point_starts_nothing_when_the_mode_question_is_cancelled(qapp, monkeypatch):
    """`main` must not fall back to a mode of its own choosing."""
    monkeypatch.setattr(gui_app, "ask_mode", lambda parent=None: None)
    monkeypatch.setattr(gui_app.MainWindow, "__init__", lambda *a, **k: pytest.fail("no window should be built"))
    assert gui_app.main([]) == 0


def test_sifi_gui_subcommand_defers_the_import(monkeypatch):
    """`sifi gui` is the only route from the CLI into Qt, and it imports the GUI
    inside the handler — a headless `sifi design` never touches PyQt5."""
    from sifi2 import cli

    called = {}

    def fake_main(argv=None, db_location=None):
        called["db"] = db_location
        return 0

    monkeypatch.setattr(gui_app, "main", fake_main)
    assert cli.main(["gui", "--db-location", "/tmp/dbs"]) == 0
    assert called["db"] == "/tmp/dbs"
