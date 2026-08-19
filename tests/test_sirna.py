"""Phase 3: ``sifi2.sirna`` against ``tests/golden/sirna_windows.json``.

The fixture is the exact content of the two temp files the original wrote, so
the windowing, the naming and the file formats are all pinned by it.
"""

from pathlib import Path

from conftest import load_golden

from sifi2 import sirna

GOLDEN = load_golden("sirna_windows.json")
SIRNA_SIZE = GOLDEN["sirna_size"]


def query_sequence():
    """The query read without Biopython, exactly as capture.py read it."""
    from conftest import DATA

    lines = (DATA / "query.fasta").read_text().splitlines()
    return "".join(line.strip() for line in lines[1:] if not line.startswith(">"))


def test_create_sirnas_reproduces_the_golden_windows():
    sirnas = sirna.create_sirnas(query_sequence(), SIRNA_SIZE, GOLDEN["sirna_start_position"])
    assert len(sirnas) == GOLDEN["count"]
    # The fixture's tab rows keep their trailing newline in the sequence field.
    expected = [(name, seq.rstrip("\n")) for name, seq in GOLDEN["rows"]]
    assert sirnas == expected


def test_windows_are_contiguous_and_named_by_position():
    seq = query_sequence()
    sirnas = sirna.create_sirnas(seq, SIRNA_SIZE)
    assert len(sirnas) == len(seq) - SIRNA_SIZE + 1
    for name, window in sirnas:
        position = int(name[len("sirna") :])
        assert window == seq[position - 1 : position - 1 + SIRNA_SIZE]
        assert len(window) == SIRNA_SIZE


def test_short_sequences_yield_nothing():
    assert sirna.create_sirnas("ACGT", 21) == []
    assert len(sirna.create_sirnas("A" * 21, 21)) == 1


def test_write_multi_fasta_matches_the_golden_files(tmp_path):
    sirnas = sirna.create_sirnas(query_sequence(), SIRNA_SIZE)
    fasta_path, tab_path = sirna.write_multi_fasta(sirnas, tmp_path / "sirnas.fasta", tmp_path / "sirnas.tab")
    assert Path(fasta_path).read_text() == GOLDEN["fasta"]
    assert Path(tab_path).read_text() == GOLDEN["tab"]


def test_fasta_headers_keep_the_space_after_the_marker(tmp_path):
    """``"> sirnaN"`` is what bowtie was run against when the fixtures were
    captured; tidying it changes the read names in bowtie's output."""
    fasta_path, _ = sirna.write_multi_fasta([("sirna1", "ACGT")], tmp_path / "f.fasta", tmp_path / "f.tab")
    assert Path(fasta_path).read_text() == "> sirna1\nACGT\n"


def test_write_single_fasta_has_no_such_space(tmp_path):
    path = sirna.write_single_fasta("TUB3_fragment", "ACGT", tmp_path / "q.fasta")
    assert Path(path).read_text() == ">TUB3_fragment\nACGT\n"


def test_read_tab_file_round_trips(tmp_path):
    sirnas = sirna.create_sirnas(query_sequence(), SIRNA_SIZE)
    _, tab_path = sirna.write_multi_fasta(sirnas, tmp_path / "f.fasta", tmp_path / "f.tab")
    assert sirna.read_tab_file(tab_path) == GOLDEN["rows"]


def test_reverse_complement_matches_biopython():
    """Pinned against the ``Seq(...).reverse_complement()`` strings captured
    under Python 2 with Biopython 1.76."""
    for case in load_golden("thermo_dangling.json")["cases"]:
        if case["sirna_sequence_n2"] is None:
            continue
        assert sirna.reverse_complement(case["sirna_sequence_n2"]) == case["reverse_complement"]


def test_reverse_complement_handles_ambiguity_and_case():
    assert sirna.reverse_complement("ACGT") == "ACGT"
    assert sirna.reverse_complement("AAAC") == "GTTT"
    assert sirna.reverse_complement("acgt") == "acgt"
    assert sirna.reverse_complement("RYSWKMN") == "NKMWSRY"


def test_lowercase_input_is_upper_cased():
    """PLAN.md Phase 6 defect 4, fixed: the original called ``sirna.upper()`` and
    threw the result away, so a lowercase or soft-masked FASTA reached the
    terminal-nucleotide rule in lower case and failed every one of its
    comparisons against ``"A"``/``"T"``/``"G"``/``"C"``."""
    seq = query_sequence()
    assert sirna.create_sirnas(seq.lower(), SIRNA_SIZE) == sirna.create_sirnas(seq, SIRNA_SIZE)
    assert sirna.create_sirnas(seq.lower(), SIRNA_SIZE)[0][1].isupper()
