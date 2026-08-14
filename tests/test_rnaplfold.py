"""Phase 3: ``sifi2.rnaplfold`` — the committed ``_lunp`` table, and a real run.

``tests/data/TUB3_fragment_lunp`` was produced by ``tests/py2_capture/make_inputs.py``
with the pipeline's own parameters, so re-running RNAplfold here must reproduce it.
"""

import numpy as np
import pytest
from conftest import DATA, load_golden

from sifi2 import rnaplfold

SIRNA_SIZE = 21
WINDOW = 8
LUNP = rnaplfold.load_lunp(DATA / "TUB3_fragment_lunp", SIRNA_SIZE)


def query():
    lines = (DATA / "query.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines[1:] if not line.startswith(">"))


def test_load_lunp_drops_the_incomplete_leading_rows():
    """The first ``sirna_size - 1`` positions have no complete window."""
    raw = np.loadtxt(DATA / "TUB3_fragment_lunp", dtype="str")
    assert len(LUNP) == len(raw) - (SIRNA_SIZE - 1)
    assert LUNP.shape[1] == SIRNA_SIZE + 1
    # Row 0 is now query position 21, the first with a full 21 nt window.
    assert LUNP[0, 0] == str(SIRNA_SIZE)


def test_accessibility_values_match_the_pipeline_json():
    """Every ``accessibility_value`` in the captured design JSON must come back
    out of the same table, by the same indexing."""
    from sifi2 import analysis

    records = analysis.load_records(DATA / "pipeline_design.json")
    for record in records:
        assert rnaplfold.accessibility_value(LUNP, record["sirna_name"], WINDOW) == record["accessibility_value"]


def test_accessibility_value_reads_the_requested_window():
    values = LUNP[0, :].astype(float).tolist()
    assert rnaplfold.accessibility_value(LUNP, "sirna1", WINDOW) == values[WINDOW]
    assert rnaplfold.accessibility_value(LUNP, "sirna1", 1) == values[1]


def test_run_rnaplfold_reproduces_the_committed_table(tmp_path):
    produced = rnaplfold.run_rnaplfold(
        "TUB3_fragment", query(), tmp_path, sirna_size=SIRNA_SIZE, winsize=80, span=40, temperature=22
    )
    assert produced.tolist() == LUNP.tolist()


def test_run_rnaplfold_raises_on_a_bad_binary(tmp_path):
    from sifi2.bowtie import BowtieError

    with pytest.raises(BowtieError, match="not found"):
        rnaplfold.run_rnaplfold("q", "ACGT" * 30, tmp_path, rnaplfold_path=str(tmp_path))


def test_golden_lunp_is_the_one_the_fixtures_used():
    """Guards against the committed table being regenerated with other settings."""
    assert load_golden("analysis.json")["sirna_size"] == SIRNA_SIZE
