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


def expected_positions(hit_name, off_target):
    """The positions a target really covers, derived straight from the records."""
    positions = set()
    for record in RECORDS:
        if record["hit_name"] == hit_name and bool(record["is_off_target"]) == off_target:
            positions |= set(range(record["sirna_position"], record["sirna_position"] + SIRNA_SIZE))
    return positions


def test_get_target_data():
    off_target_dict, main_target_dict, efficient_dict, main_hits_histo = analysis.get_target_data(RECORDS, SIRNA_SIZE)
    # Per target, not the running union the original assigned to every key
    # (PLAN.md Phase 6 defect 4). Same keys as the golden, different values.
    assert set(off_target_dict) == set(GOLDEN["off_target_dict"])
    assert set(main_target_dict) == set(GOLDEN["main_target_dict"])
    for hit_name, positions in off_target_dict.items():
        assert positions == expected_positions(hit_name, off_target=True)
    all_off_target = set().union(*off_target_dict.values())
    for hit_name, positions in main_target_dict.items():
        assert positions == expected_positions(hit_name, off_target=False) - all_off_target
    # The union is what the plot shades, and the original's keys were snapshots
    # of the running union, so their union is the same set.
    assert all_off_target == {position for positions in GOLDEN["off_target_dict"].values() for position in positions}
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


def test_each_target_is_credited_with_its_own_positions_only():
    """PLAN.md Phase 6 defect 4, fixed: the original assigned every key the
    running union of everything seen so far, so each target was credited with
    every other target's positions. The real fixture has a single main target,
    so this needs two of them to show."""
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
    # Before the fix target_b was credited with target_a's 1-21 as well.
    assert main_target_dict["target_a"] == set(range(1, 22))
    assert main_target_dict["target_b"] == set(range(50, 71))


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


def record(sirna_position, hit_name, is_off_target, is_efficient=False):
    return {
        "sirna_name": f"sirna{sirna_position}",
        "sirna_position": sirna_position,
        "is_efficient": is_efficient,
        "is_off_target": is_off_target,
        "hit_name": hit_name,
    }


def test_off_target_positions_are_per_target_too():
    """The same defect sat on ``off_target_dict``: every key held whatever the
    running union happened to be when that target's last record was read."""
    records = [record(1, "target_a", True), record(50, "target_b", True)]
    off_target_dict, _, _, _ = analysis.get_target_data(records, SIRNA_SIZE)
    assert off_target_dict["target_a"] == set(range(1, 22))
    assert off_target_dict["target_b"] == set(range(50, 71))


def test_a_sirna_hitting_two_main_targets_is_credited_to_both_but_counted_once():
    """The original's ``main_target_ready`` guard skipped the whole branch on the
    second target, so only the first was credited. The histogram still counts the
    siRNA once, which is what that guard was for."""
    records = [record(1, "target_a", False), record(1, "target_b", False)]
    _, main_target_dict, _, main_hits_histo = analysis.get_target_data(records, SIRNA_SIZE)
    assert main_target_dict == {"target_a": set(range(1, 22)), "target_b": set(range(1, 22))}
    assert main_hits_histo == list(range(1, 22))
