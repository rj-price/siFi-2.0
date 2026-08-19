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


def test_dna_tmm1_has_no_unreachable_keys():
    """PLAN.md Phase 6 defect 3, fixed: two DNA_TMM1 keys carried a leading space
    upstream, so the exact-string lookup could never reach them. Their values are
    unchanged, and the table still has one entry per key."""
    assert [key for key in thermo.DNA_TMM1 if key != key.strip()] == []
    assert thermo.DNA_TMM1["CC/GC"] == (-2.1, -5.1)
    assert thermo.DNA_TMM1["GG/CA"] == (-4.6, -11.4)
    assert len(thermo.DNA_TMM1) == 48


def test_the_terminal_mismatch_table_is_unreachable_from_sifis_own_call_shapes():
    """Why fixing those two keys moves no published number: the table is
    consulted only for a *terminal mismatch*, and siFi never presents one. Its
    dangling-end duplex pairs siRNA n with the reverse complement of siRNA n-2,
    which overlap by 19 nt and so pair exactly."""
    from sifi2 import efficiency, sirna
    from sifi2.config import SifiConfig

    class Watched(dict):
        def __init__(self, base):
            super().__init__(base)
            self.matched = 0

        def __contains__(self, key):
            found = dict.__contains__(self, key)
            self.matched += found
            return found

    watched = Watched(thermo.DNA_TMM1)
    defaults = list(thermo.calculate_free_energy.__defaults__)
    names = thermo.calculate_free_energy.__code__.co_varnames
    index = names.index("tmm_table") - 1
    original = defaults[index]
    defaults[index] = watched
    thermo.calculate_free_energy.__defaults__ = tuple(defaults)
    try:
        from conftest import DATA

        lines = (DATA / "query.fasta").read_text().splitlines()
        query = "".join(line.strip() for line in lines[1:] if not line.startswith(">"))
        config = SifiConfig()
        sirnas = sirna.create_sirnas(query, config.sirna_size)
        for position, (_name, sequence) in enumerate(sirnas, start=1):
            neighbour = None if position in (1, 2) else sirnas[position - 3][1]
            efficiency.calculate_efficiency(config, sequence, neighbour, 0.5)
    finally:
        defaults[index] = original
        thermo.calculate_free_energy.__defaults__ = tuple(defaults)

    assert watched.matched == 0


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
