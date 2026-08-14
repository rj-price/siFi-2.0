"""The RNAi efficiency rules: strand selection, end stability, accessibility.

Ported from ``legacy/sifi_pipeline.py`` (``free_energy3``,
``free_energy_dangling_ends``, ``strand_selection``, ``end_stability``,
``pair_probability``, ``check_efficient``, ``calculate_efficiency``). Every
function is free of I/O and of Qt; each takes the :class:`~sifi2.config.SifiConfig`
that supplies the thresholds and the on/off flags.

The numbers here are pinned bit-for-bit by ``tests/golden/thermo_dangling.json``,
``check_efficient.json`` and ``calculate_efficiency.json``.
"""

from __future__ import annotations

from .config import SifiConfig
from .sirna import reverse_complement
from .thermo import calculate_free_energy

__all__ = [
    "calculate_efficiency",
    "check_efficient",
    "end_stability",
    "free_energy3",
    "free_energy_dangling_ends",
    "pair_probability",
    "strand_selection",
]


def free_energy3(config: SifiConfig, sirna_sequence: str) -> tuple[float, float]:
    """ΔG of the two 5' ends of a siRNA duplex, without dangling ends.

    For a 21mer ``GGGATGGCTCAAAGGCGTAGT``: the sense 5' end is positions 0-2
    (``GGG``) and the antisense 5' end is positions 16-18 (``GTA``), i.e. the
    three nucleotides inside the 2 nt 3' overhang.
    """
    sense_five_seq = sirna_sequence[config.sirna_start_position : config.end_nucleotides]
    antisense_five_seq = sirna_sequence[
        config.sirna_size - config.overhang - config.end_nucleotides : config.sirna_size - config.overhang
    ]
    return (
        calculate_free_energy(sense_five_seq),
        calculate_free_energy(antisense_five_seq),
    )


def free_energy_dangling_ends(config: SifiConfig, sirna_sequence: str, sirna_sequence_n2: str) -> tuple[float, float]:
    """ΔG of the two 5' ends, taking the neighbouring siRNA as the dangling end.

    ``sirna_sequence_n2`` is the siRNA two positions upstream, whose reverse
    complement supplies the antisense strand of this duplex.
    """
    rc = reverse_complement(sirna_sequence_n2).strip()

    sense_five_seq = sirna_sequence[config.sirna_start_position : config.end_nucleotides]
    sense_c_seq = rc[config.sirna_size - 5 : config.sirna_size - 1]

    antisense_five_seq = rc[config.sirna_start_position : config.end_nucleotides]
    antisense_c_seq = sirna_sequence[config.sirna_size - 5 : config.sirna_size - 1]

    return (
        calculate_free_energy(sense_five_seq, check=True, strict=True, c_seq=sense_c_seq[::-1], shift=1),
        calculate_free_energy(antisense_five_seq, check=True, strict=True, c_seq=antisense_c_seq[::-1], shift=1),
    )


def _end_energies(config: SifiConfig, sirna_sequence: str, sirna_sequence_n2: str | None) -> tuple[float, float]:
    """Sense and antisense 5' ΔG, with dangling ends where they are available.

    The first two siRNAs of a query have no ``n-2`` neighbour, so they fall back
    to the plain three-nucleotide calculation.
    """
    if sirna_sequence_n2 is not None:
        return free_energy_dangling_ends(config, sirna_sequence, sirna_sequence_n2)
    return free_energy3(config, sirna_sequence)


def strand_selection(config: SifiConfig, sirna_sequence: str, sirna_sequence_n2: str | None) -> bool:
    """True when the antisense strand is the one favoured for RISC loading."""
    sense, antisense = _end_energies(config, sirna_sequence, sirna_sequence_n2)
    return antisense >= sense


def end_stability(config: SifiConfig, sirna_sequence: str, sirna_sequence_n2: str | None) -> bool:
    """True when the 5' end asymmetry reaches ``end_stability_threshold``."""
    sense, antisense = _end_energies(config, sirna_sequence, sirna_sequence_n2)
    return (antisense - sense) >= config.end_stability_threshold


def pair_probability(config: SifiConfig, lunp_data_xmer: float) -> bool:
    """True when the target site is accessible enough at the chosen window."""
    return lunp_data_xmer >= config.accessibility_threshold


def check_efficient(
    config: SifiConfig,
    strand_selection: bool,
    end_stability: bool,
    target_site_accessibility: bool,
) -> bool:
    """Combine whichever of the three rules the user enabled.

    A siRNA passes when every *enabled* rule passes; with all three disabled the
    original hard-coded ``True``, and that is kept. The original spelled this out
    as eight mutually exclusive ``if`` blocks; the collapsed form here is
    verified against all 64 combinations in ``tests/golden/check_efficient.json``.
    """
    enabled = [
        result
        for enabled_flag, result in (
            (config.strand_check, strand_selection),
            (config.end_check, end_stability),
            (config.accessibility_check, target_site_accessibility),
        )
        if enabled_flag
    ]
    return all(enabled) if enabled else True


def calculate_efficiency(
    config: SifiConfig,
    sirna_sequence: str,
    sirna_sequence_n2: str | None,
    lunp_data_xmer: float,
) -> tuple[bool | None, bool | None, bool | None, bool | None, bool | None]:
    """Score one siRNA.

    Returns ``(is_efficient, strand_selection, end_stability,
    target_site_accessibility, thermo_efficient)``, where the last is the verdict
    of the three thermodynamic rules before the terminal-nucleotide rule is
    applied on top.
    """
    if config.no_efficience:
        return False, None, None, None, None

    strand = strand_selection(config, sirna_sequence, sirna_sequence_n2)
    end = end_stability(config, sirna_sequence, sirna_sequence_n2)
    accessibility = pair_probability(config, lunp_data_xmer)
    thermo_efficient = check_efficient(config, strand, end, accessibility)

    if not config.terminal_check:
        return thermo_efficient, strand, end, accessibility, thermo_efficient

    # Terminal-nucleotide rule, on the sense strand: A/T at position 19 (of 21)
    # favours antisense loading, and the identity of position 2 then decides.
    # Note the comparisons are against uppercase literals; lowercase FASTA input
    # therefore fails every branch (``PLAN.md`` Phase 6 defect 4).
    if sirna_sequence[config.sirna_size - 3] in ("A", "T"):
        if sirna_sequence[1] in ("A", "T"):
            is_efficient = bool(thermo_efficient)
        elif config.accessibility_check:
            is_efficient = bool(accessibility)
        else:
            is_efficient = True
    elif sirna_sequence[1] in ("G", "C"):
        is_efficient = bool(thermo_efficient)
    else:
        is_efficient = False

    return is_efficient, strand, end, accessibility, thermo_efficient
