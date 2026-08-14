"""Phase 3: ``sifi2.efficiency`` against the goldens captured from the original.

``check_efficient`` is pinned exhaustively (2**3 inputs x 2**3 flags = 64 cases),
``calculate_efficiency`` across every terminal-rule branch, and the four
free-energy wrappers against the ``thermo_dangling`` fixture. All comparisons are
``repr()`` string comparisons, so the floats must match bit-for-bit.
"""

import pytest
from conftest import check_result, load_golden

from sifi2 import efficiency
from sifi2.config import SifiConfig

DANGLING = load_golden("thermo_dangling.json")
CHECK_EFFICIENT = load_golden("check_efficient.json")
CALCULATE_EFFICIENCY = load_golden("calculate_efficiency.json")


def config(**flags):
    """A config with the GUI's own defaults, as ``capture.py`` built it."""
    return SifiConfig(
        sirna_size=DANGLING["sirna_size"],
        end_nucleotides=DANGLING["end_nucleotides"],
        overhang=DANGLING["overhang"],
        end_stability_threshold=float(DANGLING["end_stability_treshold"]),
        **flags,
    )


def dangling_cases():
    return [(c["sirna_sequence"], c) for c in DANGLING["cases"]]


IDS = [name for name, _ in dangling_cases()]


@pytest.mark.parametrize("name,case", dangling_cases(), ids=IDS)
def test_free_energy3(name, case):
    """The no-dangling-ends path, used for the first two siRNAs of a query."""
    check_result(case["free_energy3"], efficiency.free_energy3, config(), case["sirna_sequence"])


@pytest.mark.parametrize("name,case", dangling_cases(), ids=IDS)
def test_free_energy_dangling_ends(name, case):
    if case["sirna_sequence_n2"] is None:
        pytest.skip("no n-2 neighbour: this case takes the free_energy3 path")
    check_result(
        case["free_energy_dangling_ends"],
        efficiency.free_energy_dangling_ends,
        config(),
        case["sirna_sequence"],
        case["sirna_sequence_n2"],
    )


@pytest.mark.parametrize("name,case", dangling_cases(), ids=IDS)
def test_strand_selection(name, case):
    check_result(
        case["strand_selection"],
        efficiency.strand_selection,
        config(),
        case["sirna_sequence"],
        case["sirna_sequence_n2"],
    )


@pytest.mark.parametrize("name,case", dangling_cases(), ids=IDS)
def test_end_stability(name, case):
    check_result(
        case["end_stability"],
        efficiency.end_stability,
        config(),
        case["sirna_sequence"],
        case["sirna_sequence_n2"],
    )


def test_pair_probability_is_inclusive_at_the_threshold():
    cfg = config()
    assert cfg.accessibility_threshold == 0.1
    assert efficiency.pair_probability(cfg, 0.1) is True
    assert efficiency.pair_probability(cfg, 0.09999) is False
    assert efficiency.pair_probability(cfg, 1.0) is True


@pytest.mark.parametrize("case", CHECK_EFFICIENT["cases"], ids=lambda c: str(sorted(c["flags"].items())))
def test_check_efficient_exhaustively(case):
    """All 64 combinations of the three rules and the three enabling flags."""
    inputs = case["inputs"]
    check_result(
        case["result"],
        efficiency.check_efficient,
        config(**case["flags"]),
        inputs["strand_selection"],
        inputs["end_stability"],
        inputs["target_site_accessibility"],
    )


def test_check_efficient_covers_all_64_combinations():
    assert len(CHECK_EFFICIENT["cases"]) == 64


def test_all_rules_disabled_means_efficient():
    cfg = config(strand_check=False, end_check=False, accessibility_check=False)
    assert efficiency.check_efficient(cfg, False, False, False) is True


@pytest.mark.parametrize(
    "case", CALCULATE_EFFICIENCY["cases"], ids=lambda c: f"{c['sirna_name']}-{sorted(c['flags'].items())}"
)
def test_calculate_efficiency(case):
    cfg = config(
        **case["flags"],
        accessibility_threshold=float(CALCULATE_EFFICIENCY["target_site_accessibility_treshold"]),
    )
    check_result(
        case["result"],
        efficiency.calculate_efficiency,
        cfg,
        case["sirna_sequence"],
        case["sirna_sequence_n2"],
        float(case["lunp_data_xmer"]),
    )


def test_no_efficience_short_circuits():
    cfg = config(no_efficience=True)
    assert efficiency.calculate_efficiency(cfg, "A" * 21, None, 1.0) == (
        False,
        None,
        None,
        None,
        None,
    )


def test_terminal_rule_branches_are_all_exercised_by_the_fixture():
    """The rule branches on the sense strand's position 19 and position 2."""
    branches = set()
    for case in CALCULATE_EFFICIENCY["cases"]:
        seq = case["sirna_sequence"]
        branches.add((seq[21 - 3] in "AT", seq[1] in "AT", seq[1] in "GC"))
    assert len(branches) >= 3
