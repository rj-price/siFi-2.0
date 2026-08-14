"""Phase 3: ``sifi2.bowtie`` — output parsing against the golden fixture, plus
the binary resolution and database management that replaced ``os.chdir()``."""

import os

import pytest
from conftest import DATA, load_golden

from sifi2 import bowtie

GOLDEN = load_golden("bowtie_parse.json")


def test_bowtie_to_lst_reproduces_the_golden_rows():
    lines = (DATA / "bowtie_output_mm2.txt").read_text().splitlines(keepends=True)
    assert len(lines) == GOLDEN["input_lines"]
    assert bowtie.bowtie_to_lst(lines) == GOLDEN["rows"]


def test_both_row_shapes_are_covered():
    """Exact hits leave the mismatch descriptor empty, so ``strip()`` drops the
    trailing tab and the row arrives with seven fields, not eight."""
    lines = (DATA / "bowtie_output_mm2.txt").read_text().splitlines(keepends=True)
    histogram = {}
    for row in bowtie.bowtie_to_lst(lines):
        histogram[str(len(row))] = histogram.get(str(len(row)), 0) + 1
    assert histogram == GOLDEN["field_count_histogram"]
    assert set(histogram) == {"8"}, "every parsed row is padded back to 8 fields"


def test_missing_mismatch_descriptor_becomes_integer_zero():
    exact = " sirna1\t+\tNM_125665.4\t300\tTTCTGGG\tIIIIIII\t1\t\n"
    mismatched = " sirna1\t+\tNM_001203444.1\t143\tTTCTGGG\tIIIIIII\t1\t11:A>G\n"
    assert bowtie.bowtie_to_lst([exact])[0][7] == 0
    assert bowtie.bowtie_to_lst([mismatched])[0][7] == "11:A>G"


def test_both_strands_are_present_in_the_fixture():
    """The reference deliberately contains a reverse-complemented paralogue, so
    the minus-strand rows Phase 6 defect 1 discards are really in the input."""
    strands = {row[1] for row in GOLDEN["rows"]}
    assert strands == {"+", "-"}


def test_resolve_binary_finds_bowtie_on_path():
    assert os.path.basename(bowtie.resolve_binary("bowtie")).startswith("bowtie")


def test_resolve_binary_honours_an_explicit_path(tmp_path):
    real = bowtie.resolve_binary("bowtie")
    assert bowtie.resolve_binary("bowtie", os.path.dirname(real)) == real
    assert bowtie.resolve_binary("bowtie", real) == real
    with pytest.raises(bowtie.BowtieError, match="not found"):
        bowtie.resolve_binary("bowtie", str(tmp_path))
    with pytest.raises(bowtie.BowtieError, match="not found on PATH"):
        bowtie.resolve_binary("definitely-not-a-real-binary")


def test_build_run_and_delete_a_real_database(tmp_path):
    """End-to-end over the real binaries: build an index, align, tear it down."""
    db_location = tmp_path / "db"
    db_location.mkdir()
    bowtie.build_database("testdb", DATA / "reference.fasta", db_location)
    assert bowtie.database_exists("testdb", db_location)

    summary = bowtie.all_databases(db_location)
    assert summary["Database name"] == ["testdb"]
    assert len(summary["Created"]) == 1
    assert isinstance(summary["Database size (MB)"][0], int)

    rows = bowtie.run_bowtie(DATA / "query.fasta", db_location / "testdb", 0, tmp_path / "out.txt", None)
    # The query is a fragment of one of the reference records, so it must hit.
    assert rows and rows[0].split("\t")[2].startswith("NM_")

    assert bowtie.delete_databases(["testdb"], db_location) == ["testdb"]
    assert not bowtie.database_exists("testdb", db_location)
    assert bowtie.delete_databases(["testdb"], db_location) == []


def test_run_bowtie_raises_instead_of_reporting_no_hits(tmp_path):
    """The original could not tell a crash from an empty result, because
    ``mkstemp`` had already created the output file it tested for."""
    with pytest.raises(bowtie.BowtieError, match="bowtie failed"):
        bowtie.run_bowtie(DATA / "query.fasta", tmp_path / "nosuchdb", 0, tmp_path / "out.txt", None)


def test_build_database_raises_on_a_missing_fasta(tmp_path):
    with pytest.raises(bowtie.BowtieError, match="bowtie-build failed"):
        bowtie.build_database("testdb", tmp_path / "nosuch.fasta", tmp_path)
