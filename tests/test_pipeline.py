"""Phase 3: ``sifi2.pipeline`` — record construction, and the broken Qt coupling.

The load-bearing test here is ``test_data_to_json_reproduces_the_captured_json``:
``tests/data/pipeline_design.json`` is real output from the Python 2 original, so
reproducing it record-for-record exercises windowing, bowtie parsing, RNAplfold
indexing and the whole efficiency chain at once.
"""

import subprocess
import sys

import pytest
from conftest import DATA

from sifi2 import analysis, bowtie, rnaplfold, sirna
from sifi2.config import Mode, SifiConfig, qt_flag
from sifi2.pipeline import PipelineCancelled, SifiPipeline

SIRNA_SIZE = 21
MAIN_TARGETS = ["NM_125665.4"]
EXPECTED = analysis.load_records(DATA / "pipeline_design.json")


def bowtie_lines():
    return (DATA / "bowtie_output_mm2.txt").read_text().splitlines(keepends=True)


def query():
    lines = (DATA / "query.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines[1:] if not line.startswith(">"))


def design_records():
    """Rebuild the captured JSON with the ported pipeline."""
    pipe = SifiPipeline(SifiConfig(mode=Mode.DESIGN))
    sirnas = sirna.create_sirnas(query(), SIRNA_SIZE)
    bowtie_data = bowtie.bowtie_to_lst(bowtie_lines())
    lunp = rnaplfold.load_lunp(DATA / "TUB3_fragment_lunp", SIRNA_SIZE)
    return pipe.data_to_json("TUB3_fragment", bowtie_data, False, lunp, MAIN_TARGETS, sirnas)


RECORDS = design_records()


PLUS_RECORDS = [record for record in RECORDS if record["strand"] == "+"]


def test_data_to_json_reproduces_the_captured_json():
    """The captured Python 2 JSON holds plus-strand hits only, since the original
    discarded the rest (Phase 6 defect 1, fixed). Its records must still be
    reproduced exactly — the fix adds records, it must not alter any."""
    assert len(PLUS_RECORDS) == len(EXPECTED)
    for produced, expected in zip(PLUS_RECORDS, EXPECTED, strict=True):
        assert produced == expected


def test_minus_strand_hits_are_reported():
    """PLAN.md Phase 6 defect 1, fixed: 1430 bowtie rows in, 1430 records out,
    where the original produced 1370. The 60 new records are the minus-strand
    hits its ``if strand == '+':`` guard discarded."""
    rows = bowtie.bowtie_to_lst(bowtie_lines())
    assert len(rows) == 1430
    assert len(RECORDS) == 1430
    assert len(PLUS_RECORDS) == 1370
    assert {record["strand"] for record in RECORDS} == {"+", "-"}


def test_a_minus_strand_record_carries_the_sirna_not_its_reverse_complement():
    """bowtie reports the read as it aligned, so column 5 of a minus-strand row
    is the reverse complement of the siRNA. Scoring that would score a sequence
    the construct never contains, so the record carries the siRNA itself."""
    by_position = dict(sirna.create_sirnas(query(), SIRNA_SIZE))
    minus = [record for record in RECORDS if record["strand"] == "-"]
    assert minus
    for record in minus:
        assert record["sirna_sequence"] == by_position[record["sirna_name"]]


def test_design_mode_without_hits_still_yields_nothing():
    """PLAN.md Phase 6 defect 2: the no-hits path sets ``strand = None`` and is
    still discarded, so the efficiency-only plot the code provides for is dead.
    Fixed in its own commit."""
    pipe = SifiPipeline(SifiConfig(mode=Mode.DESIGN))
    sirnas = sirna.create_sirnas(query(), SIRNA_SIZE)
    lunp = rnaplfold.load_lunp(DATA / "TUB3_fragment_lunp", SIRNA_SIZE)
    records = pipe.data_to_json("TUB3_fragment", [list(pair) for pair in sirnas], True, lunp, None, sirnas)
    assert records == []


def test_off_target_mode_leaves_is_off_target_unset():
    pipe = SifiPipeline(SifiConfig(mode=Mode.OFFTARGET))
    sirnas = sirna.create_sirnas(query(), SIRNA_SIZE)
    lunp = rnaplfold.load_lunp(DATA / "TUB3_fragment_lunp", SIRNA_SIZE)
    rows = bowtie.bowtie_to_lst(bowtie_lines())
    records = [r for r in pipe.data_to_json("TUB3_fragment", rows, False, lunp, None, sirnas) if r["strand"] == "+"]
    assert {record["is_off_target"] for record in records} == {None}
    # Everything else must be identical to the design-mode run.
    for produced, expected in zip(records, EXPECTED, strict=True):
        assert {k: v for k, v in produced.items() if k != "is_off_target"} == {
            k: v for k, v in expected.items() if k != "is_off_target"
        }


def test_main_targets_decide_which_hits_are_off_targets():
    by_name = {}
    for record in RECORDS:
        by_name.setdefault(record["hit_name"], set()).add(record["is_off_target"])
    assert by_name["NM_125665.4"] == {False}
    for hit_name, flags in by_name.items():
        if hit_name not in MAIN_TARGETS:
            assert flags == {True}, hit_name


def test_first_two_sirnas_have_no_dangling_end_partner():
    """siRNA n takes siRNA n-2 as its dangling end, so siRNA 1 and 2 fall back to
    the plain three-nucleotide calculation."""
    pipe = SifiPipeline(SifiConfig(mode=Mode.DESIGN))
    sirnas = sirna.create_sirnas(query(), SIRNA_SIZE)
    lunp = rnaplfold.load_lunp(DATA / "TUB3_fragment_lunp", SIRNA_SIZE)
    rows = [row for row in bowtie.bowtie_to_lst(bowtie_lines()) if row[0] in ("sirna1", "sirna2", "sirna3")]
    records = pipe.data_to_json("TUB3_fragment", rows, False, lunp, MAIN_TARGETS, sirnas)
    assert records == [r for r in EXPECTED if r["sirna_name"] in ("sirna1", "sirna2", "sirna3")]


# ----------------------------------------------------------------------
# Main-target selection: the injected callable that replaced the Qt dialog
# ----------------------------------------------------------------------
def test_get_main_targets_calls_the_selector_with_counts_descending():
    seen = []

    def selector(choices):
        seen.append(choices)
        return [choices[0][0]]

    pipe = SifiPipeline(SifiConfig(mode=Mode.DESIGN), main_target_selector=selector)
    rows = bowtie.bowtie_to_lst(bowtie_lines())
    assert pipe.get_main_targets(rows) == [seen[0][0][0]]
    counts = [count for _name, count in seen[0]]
    assert counts == sorted(counts, reverse=True)
    assert dict(seen[0])["NM_125665.4"] > 0


def test_a_cancelled_selection_raises():
    pipe = SifiPipeline(SifiConfig(mode=Mode.DESIGN), main_target_selector=lambda choices: None)
    with pytest.raises(PipelineCancelled):
        pipe.get_main_targets([["sirna1", "+", "NM_125665.4"]])


def test_design_mode_without_a_selector_raises():
    pipe = SifiPipeline(SifiConfig(mode=Mode.DESIGN))
    with pytest.raises(PipelineCancelled, match="main_target_selector"):
        pipe.get_main_targets([["sirna1", "+", "NM_125665.4"]])


# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
def test_config_defaults_are_the_guis_own():
    config = SifiConfig()
    assert (config.sirna_size, config.mismatches, config.accessibility_window) == (21, 0, 8)
    assert (config.end_stability_threshold, config.accessibility_threshold) == (1.0, 0.1)
    assert (config.winsize, config.span, config.temperature) == (80, 40, 22)
    assert (config.overhang, config.end_nucleotides, config.sirna_start_position) == (2, 3, 0)
    assert config.is_design


def test_qt_checkbox_states_normalise_to_bool():
    """Qt hands over 0/2, not 0/1."""
    assert qt_flag(2) is True
    assert qt_flag(0) is False
    assert Mode.from_legacy(0) is Mode.DESIGN
    assert Mode.from_legacy(1) is Mode.OFFTARGET


# ----------------------------------------------------------------------
# Headless
# ----------------------------------------------------------------------
def test_pipeline_imports_without_pyqt():
    """The regression test that the Qt coupling stayed broken (PLAN.md
    Verification). PyQt5 is deliberately absent from the sifi2 env, so this also
    fails loudly if something re-introduces the dependency."""
    code = (
        "import sys;"
        "sys.modules['PyQt5'] = None; sys.modules['PyQt4'] = None;"
        "import sifi2.pipeline, sifi2.efficiency, sifi2.analysis, sifi2.bowtie;"
        "print('ok')"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "ok" in result.stdout


def test_no_qt_import_anywhere_in_the_package():
    import re
    from pathlib import Path

    import sifi2

    qt_import = re.compile(r"^\s*(import|from)\s+(PyQt\d|popup|show_plot)\b", re.MULTILINE)
    for path in Path(sifi2.__file__).parent.glob("*.py"):
        assert not qt_import.search(path.read_text()), path
