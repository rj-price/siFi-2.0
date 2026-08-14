"""Phase 1 guard rails for the golden fixtures.

These do not test any ported code -- there isn't any yet. They test that the
fixtures captured from the Python 2 original in tests/golden/ are well-formed and
usable as a fidelity contract under Python 3, which is the whole premise of the
port:

  * every float was recorded as a repr() string, and Python 3 round-trips that
    string back to the identical float (Python 2.7 and Python 3 both use the
    shortest-round-trip repr, but the port depends on it, so it is pinned here);
  * the fixtures cover what PLAN.md Phase 1 says they cover;
  * the defects PLAN.md Phase 6 will later fix are present in the fixtures, so
    that fixing them shows up as a deliberate golden diff rather than a surprise.
"""

import ast
import json
from pathlib import Path

import pytest

GOLDEN = Path(__file__).parent / "golden"
DATA = Path(__file__).parent / "data"

FIXTURES = [
    "thermo_nn.json",
    "thermo_dangling.json",
    "salt_correction.json",
    "check_efficient.json",
    "calculate_efficiency.json",
    "sirna_windows.json",
    "bowtie_parse.json",
    "analysis.json",
]

# sirna_windows and analysis are pure sequence/position data: no float ever
# reaches them, so there is nothing for the round-trip guard to check.
FLOAT_FIXTURES = [
    "thermo_nn.json",
    "thermo_dangling.json",
    "salt_correction.json",
    "calculate_efficiency.json",
]


def load(name):
    with open(GOLDEN / name) as fh:
        return json.load(fh)


def iter_strings(obj):
    """Yield every string anywhere in a nested structure."""
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for value in obj.values():
            yield from iter_strings(value)
    elif isinstance(obj, list):
        for value in obj:
            yield from iter_strings(value)


@pytest.mark.parametrize("name", FIXTURES)
def test_fixture_exists_and_is_documented(name):
    data = load(name)
    assert data["description"], f"{name} must carry a description"


@pytest.mark.parametrize("name", FIXTURES)
def test_no_bare_json_floats(name):
    """Floats must be repr() strings, never JSON numbers -- JSON numbers would
    lose the exactness the port is being held to."""

    def walk(obj, path="$"):
        if isinstance(obj, float):
            pytest.fail(f"{name}: bare JSON float at {path}: {obj!r}")
        elif isinstance(obj, dict):
            for k, v in obj.items():
                walk(v, f"{path}.{k}")
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                walk(v, f"{path}[{i}]")

    walk(load(name))


@pytest.mark.parametrize("name", FLOAT_FIXTURES)
def test_recorded_reprs_round_trip_under_python3(name):
    """Every recorded repr() must evaluate under Python 3 and re-repr identically.

    This is the load-bearing assumption of the whole fixture scheme: a Python 2.7
    repr of a float, a bool, or a tuple of them must be exactly what Python 3
    produces for the same value.
    """
    checked = 0
    for text in iter_strings(load(name)):
        try:
            value = ast.literal_eval(text)
        except (ValueError, SyntaxError):
            continue  # a plain string (a sequence, a hit name, a description)
        if not isinstance(value, float | bool | tuple | int) or isinstance(value, str):
            continue
        assert repr(value) == text, f"{name}: {text!r} re-reprs as {repr(value)!r}"
        checked += 1
    assert checked > 0, f"{name}: found nothing to round-trip"


def test_thermo_nn_covers_all_64_trinucleotides():
    seqs = {c["seq"] for c in load("thermo_nn.json")["cases"]}
    expected = {a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"}
    assert expected <= seqs
    assert len(expected) == 64


def test_check_efficient_is_exhaustive():
    cases = load("check_efficient.json")["cases"]
    assert len(cases) == 64
    seen = {
        (
            c["flags"]["strand_check"],
            c["flags"]["end_check"],
            c["flags"]["accessibility_check"],
            c["inputs"]["strand_selection"],
            c["inputs"]["end_stability"],
            c["inputs"]["target_site_accessibility"],
        )
        for c in cases
    }
    assert len(seen) == 64, "all 2**3 inputs x 2**3 flags must be distinct"
    # The original never returns None: the all-flags-off branch hard-codes True.
    assert all(c["result"]["value"] in ("True", "False") for c in cases)


def test_salt_correction_covers_every_method_including_the_gc_ones():
    cases = load("salt_correction.json")["cases"]
    methods = {c["kwargs"]["method"] for c in cases}
    assert set(range(0, 9)) <= methods
    # Methods 6 and 7 are the SeqUtils.GC ones -- the 100x trap in PLAN.md Risks.
    for method in (6, 7):
        assert any(c["kwargs"]["method"] == method and "value" in c["result"] for c in cases), (
            f"method {method} must have at least one non-error case to pin"
        )


def test_calculate_efficiency_covers_every_terminal_rule_branch():
    cases = load("calculate_efficiency.json")["cases"]
    branches = set()
    for c in cases:
        s = c["sirna_sequence"]
        branches.add((s[21 - 3] in "AT", s[1] in "AT", s[1] in "GC"))
    # A/T at S19 with A/T at S1; A/T at S19 without; not-A/T at S19 with G/C at
    # S1; and not-A/T at S19 without -- the four leaves of the terminal rule.
    assert len(branches) >= 4, branches
    assert any(c["flags"]["no_efficience"] for c in cases)
    assert any(not c["flags"]["terminal_check"] for c in cases)


def test_bowtie_parse_covers_both_row_shapes_and_both_strands():
    data = load("bowtie_parse.json")
    rows = data["rows"]
    assert len(rows) == data["input_lines"]
    assert all(len(r) == 8 for r in rows), "every row is padded to 8 fields"
    # Exact hits lose their empty mismatch column to .strip() and get int 0.
    assert any(r[7] == 0 for r in rows)
    assert any(isinstance(r[7], str) and r[7] for r in rows)
    strands = {r[1] for r in rows}
    assert strands == {"+", "-"}, "fixture must exercise minus-strand hits"


def test_sirna_windows_match_the_query():
    data = load("sirna_windows.json")
    rows = data["rows"]
    assert data["count"] == data["query_length"] - data["sirna_size"] + 1
    assert len(rows) == data["count"]
    assert rows[0][0] == "sirna1"
    assert all(len(r[1].strip()) == data["sirna_size"] for r in rows)
    # The upstream FASTA header really is '> sirnaN', with a space.
    assert data["fasta"].startswith("> sirna1\n")


def test_analysis_pins_the_minus_strand_defect():
    """PLAN.md Phase 6 defect 1: data_to_json builds its record inside
    `if strand == '+':`, so minus-strand bowtie hits are silently discarded.

    The fixtures must show the defect, so that Phase 6's fix produces a visible,
    attributable golden diff rather than a silent change.
    """
    analysis = load("analysis.json")
    bowtie = load("bowtie_parse.json")
    plus = sum(1 for r in bowtie["rows"] if r[1] == "+")
    minus = sum(1 for r in bowtie["rows"] if r[1] == "-")
    assert minus > 0
    assert analysis["record_count"] == plus, (
        "the golden must record the buggy behaviour: "
        f"{minus} minus-strand hits dropped from {bowtie['input_lines']} bowtie rows"
    )


def test_analysis_inputs_are_committed():
    for name in (
        "query.fasta",
        "query_multi.fasta",
        "reference.fasta",
        "bowtie_output_mm2.txt",
        "TUB3_fragment_lunp",
        "pipeline_design.json",
    ):
        assert (DATA / name).is_file(), f"tests/data/{name} is missing"
