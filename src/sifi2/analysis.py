"""Aggregation of pipeline records into tables and per-target position sets.

Ported from ``legacy/general_helpers.py``. These functions take the records
themselves — the original re-read and re-parsed the JSON file on every call, via
a hand-rolled ``iterparse`` that worked around a limitation ``json.load`` does
not have. :func:`load_records` is the file-reading convenience wrapper.

The aggregation semantics are pinned by ``tests/golden/analysis.json``, except
where ``PLAN.md`` Phase 6 fixed them deliberately; see :func:`get_target_data`.
"""

from __future__ import annotations

import json
import os
import re
import time
from collections import Counter
from itertools import groupby

__all__ = [
    "create_gbk",
    "get_table_data",
    "get_target_data",
    "group_ranges",
    "load_records",
    "validate_fasta_seq",
    "validate_seq",
]


def load_records(f_in: str | os.PathLike[str]) -> list[dict]:
    """Read a pipeline JSON file into a list of records."""
    with open(f_in) as handle:
        return json.load(handle)


def get_table_data(records: list[dict]) -> list[list]:
    """Per-target summary rows: ``[hit_name, total hits, efficient hits]``."""
    hit_counter = Counter(record["hit_name"] for record in records)
    efficient_counter = Counter(record["hit_name"] for record in records if record["is_efficient"])
    table_data = []
    for hit_name, hits in hit_counter.most_common():
        for efficient_name, efficient_hits in efficient_counter.most_common():
            if hit_name == efficient_name:
                table_data.append([hit_name, hits, efficient_hits])
        if hit_name not in efficient_counter:
            table_data.append([hit_name, hits, 0])
    return table_data


def get_target_data(records: list[dict], sirna_size: int) -> tuple[dict, dict, dict, list]:
    """Positions covered per target, for the design plot.

    Returns ``(off_target_dict, main_target_dict, efficient_dict,
    main_hits_histo)``: covered query positions per off-target and per main
    target, the positions of efficient siRNAs per target, and the flat list of
    main-target positions the histogram is built from.

    Both dictionaries are per-target. The original assigned each key the running
    union of everything seen so far, so every target was credited with every
    other target's positions, and each key held a different snapshot of that
    union depending on where its last record fell (``PLAN.md`` Phase 6 defect 4).

    Main-target positions have the off-target positions subtracted, as the
    original intended: a position that is also hit off-target is shown as an
    off-target on the design plot. The subtraction now uses the *complete* set of
    off-target positions rather than however much of it had accumulated by then.
    """
    off_target_positions: set[int] = set()
    off_target_dict: dict[str, set[int]] = {}
    main_target_dict: dict[str, set[int]] = {}
    efficient_dict: dict[str, list[int]] = {}
    main_hits_histo: list[int] = []
    ready_sirnas: set[str] = set()
    main_target_ready: set[tuple] = set()

    for data in records:
        start = int(data["sirna_position"])
        plot_range_l = list(range(start, start + sirna_size))
        plot_range = set(plot_range_l)

        if data["is_efficient"] and data["sirna_name"] not in ready_sirnas:
            if data["hit_name"] in efficient_dict:
                efficient_dict[data["hit_name"]].extend(plot_range_l)
            else:
                efficient_dict[data["hit_name"]] = plot_range_l
            ready_sirnas.add(data["sirna_name"])

        if data["is_off_target"]:
            off_target_positions |= plot_range
            off_target_dict.setdefault(data["hit_name"], set()).update(plot_range)
        else:
            main_target_dict.setdefault(data["hit_name"], set()).update(plot_range)
            # The histogram counts each siRNA once however many main targets it
            # hits. The original extended it from the *set*, so the order was
            # CPython's set-iteration order and differs between Python 2 and 3.
            # Ascending order is used instead: same multiset, and the only
            # consumer is a histogram, where order cannot matter.
            if (data["sirna_name"], data["sirna_position"]) not in main_target_ready:
                main_target_ready.add((data["sirna_name"], data["sirna_position"]))
                main_hits_histo.extend(plot_range_l)

    main_target_dict = {name: positions - off_target_positions for name, positions in main_target_dict.items()}
    return off_target_dict, main_target_dict, efficient_dict, main_hits_histo


def group_ranges(data) -> list[tuple[int, int]]:
    """Collapse a sorted-ish sequence of positions into ``(first, last)`` runs."""
    ranges = []
    for _key, group in groupby(enumerate(data), lambda ix: ix[0] - ix[1]):
        positions = [item[1] for item in group]
        ranges.append((positions[0], positions[-1]))
    return ranges


def validate_seq(sequence: str) -> bool:
    """Validate a plain-text nucleotide sequence with no FASTA header."""
    sequence = sequence.strip().replace(" ", "").replace("\n", "")
    return re.compile("^[ACTGNRYSWKMBDHVEFILPQSXZ]*$", re.I).search(sequence) is not None


def validate_fasta_seq(sequence: str) -> int | bool:
    """Number of FASTA records in ``sequence``, or ``False`` if there are none."""
    sequence = sequence.replace(" ", "")
    matches = re.compile(">\\S*\n[ACTGNRYSWKMBDHVEFILPQSXZ]*", re.MULTILINE).findall(sequence)
    return len(matches) if matches else False


def create_gbk(
    main_target_dict: dict,
    off_target_dict: dict,
    query_name: str,
    query_sequence: str,
    out_file: str | os.PathLike[str],
) -> str:
    """Write the query with main- and off-target hits as GenBank features.

    Takes the sequence directly rather than copying a temp file to ``.fasta`` and
    re-reading it, and sets ``molecule_type``, which modern ``SeqIO.write`` needs
    for the GenBank format (``Bio.Alphabet`` was removed in Biopython 1.78).
    """
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqIO import write
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq(query_sequence), id=query_name, name=query_name, description="")
    record.annotations["molecule_type"] = "DNA"

    for prefix, targets in (("MT ", main_target_dict), ("OT ", off_target_dict)):
        for target, positions in targets.items():
            for first, last in group_ranges(sorted(positions)):
                record.features.append(SeqFeature(FeatureLocation(first, last), type=prefix + str(target)))

    record.annotations["date"] = time.strftime("%d") + "-" + time.strftime("%b").upper() + "-" + time.strftime("%Y")
    with open(out_file, "w") as handle:
        write(record, handle, "genbank")
    return str(out_file)
