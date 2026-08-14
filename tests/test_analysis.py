"""Phase 3: ``sifi2.analysis`` against ``tests/golden/analysis.json``.

The fixture aggregates ``tests/data/pipeline_design.json``, which is real
``data_to_json`` output from the Python 2 original for design mode with
``main_targets=['NM_125665.4']``.
"""

from collections import Counter

from conftest import DATA, load_golden

from sifi2 import analysis

GOLDEN = load_golden("analysis.json")
SIRNA_SIZE = GOLDEN["sirna_size"]
RECORDS = analysis.load_records(DATA / "pipeline_design.json")


def test_records_load():
    assert len(RECORDS) == GOLDEN["record_count"]


def test_get_table_data():
    assert sorted(analysis.get_table_data(RECORDS)) == sorted(GOLDEN["table_data"])


def test_get_table_data_counts_are_consistent():
    for hit_name, hits, efficient in analysis.get_table_data(RECORDS):
        assert efficient <= hits
        assert hits == sum(1 for record in RECORDS if record["hit_name"] == hit_name)


def test_get_target_data():
    off_target_dict, main_target_dict, efficient_dict, main_hits_histo = analysis.get_target_data(RECORDS, SIRNA_SIZE)
    assert {k: sorted(v) for k, v in off_target_dict.items()} == GOLDEN["off_target_dict"]
    assert {k: sorted(v) for k, v in main_target_dict.items()} == GOLDEN["main_target_dict"]
    assert {k: list(v) for k, v in efficient_dict.items()} == GOLDEN["efficient_dict"]
    # Same multiset, deliberately different order: see below.
    assert sorted(main_hits_histo) == sorted(GOLDEN["main_hits_histo"])


def test_main_hits_histo_is_ascending_not_set_ordered():
    """The original extended this list from a ``set``, so its order was CPython's
    set-iteration order, which is not the same in Python 2 and 3. The port emits
    ascending positions instead; the only consumer is a histogram."""
    _, _, _, main_hits_histo = analysis.get_target_data(RECORDS, SIRNA_SIZE)
    assert main_hits_histo != GOLDEN["main_hits_histo"]
    assert Counter(main_hits_histo) == Counter(GOLDEN["main_hits_histo"])
    chunks = [main_hits_histo[i : i + SIRNA_SIZE] for i in range(0, len(main_hits_histo), SIRNA_SIZE)]
    assert all(chunk == sorted(chunk) for chunk in chunks)


def test_main_target_dict_still_aliases_one_set():
    """PLAN.md Phase 6 defect 4: one accumulating set is shared across all keys,
    so each target is credited with every position seen so far, not its own. The
    real fixture has a single main target, so this needs two of them to show.
    Kept until the fix lands as its own commit."""
    records = [
        {
            "sirna_name": "sirna1",
            "sirna_position": 1,
            "is_efficient": False,
            "is_off_target": False,
            "hit_name": "target_a",
        },
        {
            "sirna_name": "sirna50",
            "sirna_position": 50,
            "is_efficient": False,
            "is_off_target": False,
            "hit_name": "target_b",
        },
    ]
    _, main_target_dict, _, _ = analysis.get_target_data(records, SIRNA_SIZE)
    assert set(main_target_dict) == {"target_a", "target_b"}
    # target_b hits at 50-70 only, but is credited with target_a's 1-21 too.
    assert main_target_dict["target_a"] == set(range(1, 22))
    assert main_target_dict["target_b"] == set(range(1, 22)) | set(range(50, 71))


def test_group_ranges():
    assert analysis.group_ranges([1, 2, 3, 7, 8, 20]) == [(1, 3), (7, 8), (20, 20)]
    assert analysis.group_ranges([]) == []
    assert analysis.group_ranges([5]) == [(5, 5)]


def test_validate_seq():
    assert analysis.validate_seq("ACGT acgt\nNRY") is True
    assert analysis.validate_seq("") is True
    assert analysis.validate_seq("ACGT!") is False


def test_validate_fasta_seq():
    assert analysis.validate_fasta_seq(">one\nACGT\n") == 1
    assert analysis.validate_fasta_seq(">one\nACGT\n>two\nACGT\n") == 2
    assert analysis.validate_fasta_seq("ACGT") is False


def test_create_gbk_writes_features(tmp_path):
    """``Bio.Alphabet`` is gone; GenBank output now needs ``molecule_type``."""
    out = tmp_path / "out.gbk"
    analysis.create_gbk(
        {"NM_125665.4": {1, 2, 3, 10, 11}},
        {"NM_125664.4": {5, 6, 7}},
        "TUB3_fragment",
        "ACGT" * 10,
        out,
    )
    text = out.read_text()
    # Biopython replaces the space in the feature key with an underscore.
    assert "MT_NM_125665.4" in text
    assert "OT_NM_125664.4" in text
    assert "acgtacgt" in text.lower()
