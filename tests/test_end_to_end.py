"""Phase 5: the end-to-end baseline.

The per-function golden fixtures in ``tests/golden/`` carry the fidelity
guarantee against the Python 2 original — a whole-pipeline comparison against it
was never possible, since its design mode blocks on a Qt dialog and the binaries
it shipped with are Windows ``.exe``\\ s. What is possible, and what this file
does, is to pin the **port's own** whole-pipeline output: a bowtie index built
from the committed reference FASTA, both modes run over the committed query, and
the resulting records compared record-for-record against ``tests/baseline/``.

That makes the Phase 6 defect fixes show up here as a reviewed diff, exactly as
they will in ``tests/golden/`` — see ``tests/baseline/README.md`` for how to
regenerate these two files.
"""

import json
import shutil

import pytest
from conftest import BASELINE, DATA

from sifi2 import bowtie
from sifi2.cli import main

needs_binaries = pytest.mark.skipif(
    shutil.which("bowtie") is None or shutil.which("bowtie-build") is None or shutil.which("RNAplfold") is None,
    reason="needs bowtie, bowtie-build and RNAplfold on PATH",
)

pytestmark = needs_binaries


@pytest.fixture(scope="module")
def test_db(tmp_path_factory):
    db_location = tmp_path_factory.mktemp("dbs")
    bowtie.build_database("testdb", str(DATA / "reference.fasta"), str(db_location))
    return str(db_location)


def run(mode, test_db, outdir, *extra):
    argv = [
        mode, "--query", str(DATA / "query.fasta"), "--db", "testdb",
        "--db-location", test_db, "--outdir", str(outdir), "--quiet", *extra,
    ]  # fmt: skip
    assert main(argv) == 0
    return json.loads((outdir / "TUB3_fragment.json").read_text())


def load_baseline(name):
    with open(BASELINE / name) as handle:
        return json.load(handle)


def test_offtarget_run_reproduces_the_baseline(test_db, tmp_path):
    records = run("offtarget", test_db, tmp_path)
    assert records == load_baseline("offtarget_TUB3_fragment.json")


def test_design_run_reproduces_the_baseline(test_db, tmp_path):
    records = run("design", test_db, tmp_path, "--all-targets-main")
    assert records == load_baseline("design_TUB3_fragment.json")


def test_the_two_modes_differ_only_in_efficiency_scoring_and_target_attribution(test_db, tmp_path):
    """Design mode turns the end-stability and accessibility rules on, and
    attributes every hit to a main target or an off-target; off-target mode
    reports the same hits with neither."""
    design = run("design", test_db, tmp_path / "design", "--all-targets-main")
    offtarget = run("offtarget", test_db, tmp_path / "offtarget")

    keys = ("sirna_name", "hit_name", "reference_strand_pos", "strand", "mismatches")
    assert [tuple(record[key] for key in keys) for record in design] == [
        tuple(record[key] for key in keys) for record in offtarget
    ]
    assert {record["is_off_target"] for record in design} == {False}
    assert {record["is_off_target"] for record in offtarget} == {None}
    # Design applies two rules off-target mode does not, so it is stricter.
    assert sum(record["is_efficient"] for record in design) < sum(record["is_efficient"] for record in offtarget)


def test_every_output_file_is_written(test_db, tmp_path):
    run("offtarget", test_db, tmp_path)
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "TUB3_fragment.json",
        "TUB3_fragment.png",
        "TUB3_fragment.tsv",
    ]
