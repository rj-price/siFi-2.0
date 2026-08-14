"""siRNA windowing and the FASTA/tab files handed to bowtie.

Ported from ``legacy/sifi_pipeline.py`` (``create_sirnas``,
``create_single_fasta_file``, ``create_multi_fasta_file``). The file formats are
part of the contract with the golden fixtures — in particular the multi-FASTA
headers are ``"> sirnaN"``, with a space after the ``>``.
"""

from __future__ import annotations

import os

__all__ = [
    "create_sirnas",
    "read_tab_file",
    "reverse_complement",
    "write_multi_fasta",
    "write_single_fasta",
]

_COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA sequence.

    Replaces ``Bio.Seq.Seq(...).reverse_complement()``; kept local so the core
    stays independent of Biopython's ``Seq``, which since 1.78 no longer
    subclasses ``str`` (``PLAN.md`` Risks).
    """
    return seq.translate(_COMPLEMENT)[::-1]


def create_sirnas(query_sequence: str, sirna_size: int, start_position: int = 0) -> list[tuple[str, str]]:
    """Slide a ``sirna_size`` window over ``query_sequence``.

    Returns ``[(name, sequence)]`` with names ``sirna1``, ``sirna2``, ... — the
    name carries the 1-based start position on the query.

    The original called ``sirna.upper()`` and discarded the result; that no-op is
    preserved here (``PLAN.md`` Phase 6 defect 4) because the terminal-nucleotide
    rule compares against uppercase literals, so lowercase input changes the
    result and fixing it is a deliberate, attributable change.
    """
    start = start_position
    end = start_position + sirna_size
    sirnas = []
    while len(query_sequence[start:end]) == sirna_size:
        sirnas.append((f"sirna{start + 1}", query_sequence[start:end]))
        start += 1
        end += 1
    return sirnas


def write_single_fasta(query_name: str, query_sequence: str, path: str | os.PathLike[str]) -> str:
    """Write one sequence as FASTA, as RNAplfold's input."""
    with open(path, "w") as handle:
        handle.write(f">{query_name}\n")
        handle.write(f"{query_sequence}\n")
    return str(path)


def write_multi_fasta(
    sirnas: list[tuple[str, str]],
    fasta_path: str | os.PathLike[str],
    tab_path: str | os.PathLike[str],
) -> tuple[str, str]:
    """Write the siRNAs as both multi-FASTA (for bowtie) and tab (for lookup).

    The ``"> name"`` header keeps the original's space after the ``>``: it ends
    up in bowtie's read names, and the fixtures were captured with it.
    """
    with open(fasta_path, "w") as fasta, open(tab_path, "w") as tab:
        for name, sequence in sirnas:
            fasta.write(f"> {name}\n")
            fasta.write(f"{sequence}\n")
            tab.write(f"{name}\t{sequence}\n")
    return str(fasta_path), str(tab_path)


def read_tab_file(path: str | os.PathLike[str]) -> list[list[str]]:
    """Read a tab file back as ``[[name, sequence], ...]``.

    The original split each raw line, so every sequence carried a trailing
    newline and each use site had to remember to ``strip()`` it; the newline is
    dropped here instead.
    """
    with open(path) as handle:
        return [line.split("\t") for line in handle.read().splitlines()]
