"""Phase 2: ``sifi2.thermo`` against the goldens captured from the Python 2 original.

Every case is compared as ``repr()`` strings, exactly as captured, so a match is
bit-for-bit and not "close enough" — silent numerical drift is the main risk of
this port (``PLAN.md`` Risks).
"""

import json
from pathlib import Path

import pytest

from sifi2 import thermo

GOLDEN = Path(__file__).parent / "golden"


def load(name):
    with open(GOLDEN / name) as fh:
        return json.load(fh)


def check(expected, fn, *args, **kwargs):
    """Run ``fn`` and assert it reproduces a captured ``{"value"|"error"}`` result."""
    if "error" in expected:
        with pytest.raises(Exception) as excinfo:  # noqa: B017 - the type is in the fixture
            fn(*args, **kwargs)
        assert f"{type(excinfo.value).__name__}: {excinfo.value}" == expected["error"]
    else:
        assert repr(fn(*args, **kwargs)) == expected["value"]


def nn_cases():
    return [(c["seq"], c["result"]) for c in load("thermo_nn.json")["cases"]]


def dangling_cases():
    out = []
    for case in load("thermo_dangling.json")["cases"]:
        for strand in ("sense", "antisense"):
            if strand in case:
                sub = case[strand]
                out.append((f"{case['sirna_sequence']}-{strand}", sub))
    return out


def salt_cases():
    return [(c["kwargs"], c["result"]) for c in load("salt_correction.json")["cases"]]


@pytest.mark.parametrize("seq,expected", nn_cases(), ids=[s for s, _ in nn_cases()])
def test_calculate_free_energy_default_call(seq, expected):
    """The single-argument call shape: all 64 trinucleotides plus real 21mers."""
    check(expected, thermo.calculate_free_energy, seq)


@pytest.mark.parametrize("name,case", dangling_cases(), ids=[n for n, _ in dangling_cases()])
def test_calculate_free_energy_dangling_ends(name, case):
    """The pipeline's second call shape: ``c_seq=..., shift=1``."""
    check(
        case["result"],
        thermo.calculate_free_energy,
        case["seq"],
        check=True,
        strict=True,
        c_seq=case["c_seq"],
        shift=case["shift"],
    )


@pytest.mark.parametrize("kwargs,expected", salt_cases())
def test_salt_correction(kwargs, expected):
    """Methods 0-8 across four sequences and five ion combinations."""
    kw = {k: (v if k in ("method", "seq") else float(v)) for k, v in kwargs.items()}
    check(expected, thermo.salt_correction, **kw)


def test_dna_tmm1_leading_space_keys_are_preserved():
    """PLAN.md Phase 6 defect 3: two DNA_TMM1 keys carry a leading space, so the
    exact-string lookup can never reach them. Keep the typo until Phase 6 removes
    it deliberately, with the golden diff attached."""
    assert " CC/GC" in thermo.DNA_TMM1
    assert " GG/CA" in thermo.DNA_TMM1
    assert "CC/GC" not in thermo.DNA_TMM1
    assert "GG/CA" not in thermo.DNA_TMM1


def test_sifi_defaults_are_not_biopythons():
    """The siFi-specific table selection: RNA nearest neighbours with DNA
    mismatch tables and RNA dangling ends, at Na=20 K=50 saltcorr=5."""
    import inspect

    defaults = {
        name: param.default for name, param in inspect.signature(thermo.calculate_free_energy).parameters.items()
    }
    assert defaults["nn_table"] is thermo.RNA_NN3
    assert defaults["tmm_table"] is thermo.DNA_TMM1
    assert defaults["imm_table"] is thermo.DNA_IMM1
    assert defaults["de_table"] is thermo.RNA_DE2
    assert (defaults["Na"], defaults["K"], defaults["saltcorr"]) == (20, 50, 5)


def test_gc_percent_matches_biopython_1_76():
    """``_gc_percent`` must stay a percentage computed as ``gc * 100.0 / len``:
    the salt-correction expressions divide by 100 again, and going via a fraction
    is not always the same float."""
    seq = "GGGATGGCTCAAAGGCGTAGT"  # 12 G/C in 21 nt
    assert thermo._gc_percent(seq) == 12 * 100.0 / 21
    # 3 G/C in 9 nt is one of the combinations where the two routes disagree.
    assert thermo._gc_percent("GCGAAAAAA") / 100 != 3 / 9
    assert thermo._gc_percent("") == 0.0
    assert thermo._gc_percent("AAAA") == 0.0


def test_complement_and_back_transcribe():
    assert thermo._complement("ACGT") == "TGCA"
    assert thermo._complement("ACGU") == "UGCA"  # RNA table when U is present
    assert thermo._complement("RYKMBVDHNSW") == "YRMKVBHDNSW"
    assert thermo._back_transcribe("ACGU") == "ACGT"
    with pytest.raises(ValueError, match="Mixed RNA/DNA"):
        thermo._complement("ACGTU")


def test_check_rejects_an_unknown_method():
    """The original left ``baseset`` unbound here and died with UnboundLocalError."""
    assert thermo._check("acg u", "Tm_NN") == "ACGT"
    assert thermo._check("acg u", "Tm_Wallace") == "ACGT"
    with pytest.raises(ValueError, match="unknown method"):
        thermo._check("ACGT", "Tm_Bogus")
