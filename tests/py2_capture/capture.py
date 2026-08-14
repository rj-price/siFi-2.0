# -*- coding: utf-8 -*-
"""Capture golden fixtures from the Python 2 original (PLAN.md Phase 1).

Run under the throwaway py2.7 env, from anywhere:

    conda env create -f environment-py2.yml          # once
    conda run -n sifi2-py2 python tests/py2_capture/capture.py

Writes tests/golden/*.json. Re-running must reproduce the committed files
byte-for-byte -- that is the check that the goldens were never hand-edited.

Everything here runs against legacy/ unmodified. legacy/ uses Python 2 implicit
relative imports, so legacy/ goes on sys.path and its modules are imported as
top-level ones; the PyQt4 modules (popup, show_plot) are stubbed by placing this
directory *ahead* of legacy/ on sys.path.

Determinism rules, because these files are committed and diffed:
  * floats are stored as repr() strings, never rounded and never as JSON numbers
  * sets are serialised sorted
  * every dict is dumped with sort_keys=True
  * exceptions are recorded as {"error": "<type>: <message>"} so the port has to
    reproduce the failure, not just the successes
"""

import io
import json
import os
import sys
from contextlib import contextmanager

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(REPO, "tests", "data")
GOLDEN = os.path.join(REPO, "tests", "golden")

# Stubs first, then legacy -- order matters.
sys.path.insert(0, os.path.join(REPO, "legacy"))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

import free_energy  # noqa: E402
import general_helpers  # noqa: E402
from sifi_pipeline import SifiPipeline  # noqa: E402

SIRNA_SIZE = 21
ACCESSIBILITY_WINDOW = 8
QUERY_FASTA = os.path.join(DATA, "query.fasta")
BOWTIE_OUT = os.path.join(DATA, "bowtie_output_mm2.txt")
LUNP = os.path.join(DATA, "TUB3_fragment_lunp")


@contextmanager
def quiet():
    """Swallow the `print seq` in calculate_free_energy and `print 'ok'` in
    data_to_json, which otherwise bury the capture in noise."""
    saved = sys.stdout
    sys.stdout = io.BytesIO()
    try:
        yield
    finally:
        sys.stdout = saved


def dump(name, obj):
    path = os.path.join(GOLDEN, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print("wrote %-28s %8d bytes" % (name, os.path.getsize(path)))


def call(fn, *args, **kwargs):
    """Return repr() of the result, or a recorded error."""
    try:
        with quiet():
            return {"value": repr(fn(*args, **kwargs))}
    except Exception as exc:  # noqa: BLE001 - recording the failure is the point
        return {"error": "%s: %s" % (type(exc).__name__, exc)}


def make_pipeline(**overrides):
    """SifiPipeline with the GUI's own defaults; overrides by keyword."""
    kw = {
        "sirna_size": SIRNA_SIZE,
        "mismatches": 0,
        "accessibility_check": True,
        "accessibility_window": ACCESSIBILITY_WINDOW,
        "mode": 0,
        "strand_check": True,
        "end_check": True,
        "end_stability_treshold": 1.0,
        "target_site_accessibility_treshold": 0.1,
        "terminal_check": True,
        "no_efficience": False,
    }
    kw.update(overrides)
    return SifiPipeline(
        "testdb", DATA + "/", QUERY_FASTA, kw["sirna_size"], kw["mismatches"],
        kw["accessibility_check"], kw["accessibility_window"], "", "", kw["mode"],
        kw["strand_check"], kw["end_check"], kw["end_stability_treshold"],
        kw["target_site_accessibility_treshold"], "", kw["terminal_check"],
        kw["no_efficience"],
    )


def read_query():
    """The query sequence, read without Biopython so this cannot drift."""
    lines = open(QUERY_FASTA).read().splitlines()
    return "".join(line.strip() for line in lines[1:] if not line.startswith(">"))


def sirna_table(query_sequence):
    """[(name, sequence)] exactly as create_sirnas windows them."""
    out = []
    start, end = 0, SIRNA_SIZE
    for _ in range(len(query_sequence)):
        if len(query_sequence[start:end]) == SIRNA_SIZE:
            out.append(("sirna%d" % (start + 1), query_sequence[start:end]))
            start += 1
            end += 1
    return out


# --------------------------------------------------------------------------
# thermo_nn.json -- calculate_free_energy(seq), the default single-argument call
# --------------------------------------------------------------------------
def capture_thermo_nn(sirnas):
    bases = "ACGT"
    seqs = [a + b + c for a in bases for b in bases for c in bases]
    assert len(seqs) == 64
    # 50 real 21mers, evenly spread across the query rather than the first 50.
    step = max(1, len(sirnas) // 50)
    seqs += [s for _, s in sirnas[::step]][:50]
    # The 3mers the pipeline actually feeds it via free_energy3: the sense and
    # antisense 5' ends of every real siRNA.
    ends = set()
    for _, s in sirnas:
        ends.add(s[0:3])
        ends.add(s[SIRNA_SIZE - 2 - 3:SIRNA_SIZE - 2])
    seqs += sorted(ends)

    cases = []
    for seq in seqs:
        cases.append({"seq": seq, "result": call(free_energy.calculate_free_energy, seq)})
    dump("thermo_nn.json", {"description":
                            "free_energy.calculate_free_energy(seq) with all defaults "
                            "(RNA_NN3 / DNA_TMM1 / DNA_IMM1 / RNA_DE2, Na=20 K=50 saltcorr=5). "
                            "Floats are repr() strings.",
                            "cases": cases})


# --------------------------------------------------------------------------
# thermo_dangling.json -- the c_seq/shift=1 call shape, plus the Biopython Seq
# slicing that derives its arguments (PLAN.md Risks: Seq no longer subclasses str)
# --------------------------------------------------------------------------
def capture_thermo_dangling(sirnas):
    pipe = make_pipeline()
    cases = []
    # siRNA n uses siRNA n-2 as the dangling-end partner, so pairs start at index 2.
    step = max(1, (len(sirnas) - 2) // 60)
    for i in range(2, len(sirnas), step):
        sirna_sequence = sirnas[i][1]
        sirna_sequence_n2 = sirnas[i - 2][1]

        # Argument derivation copied from free_energy_dangling_ends, recorded as
        # strings so the port's replacement for Bio.Seq is pinned too.
        rc = Seq(sirna_sequence_n2).reverse_complement().strip()
        sense_five_seq = sirna_sequence[0:pipe.end_nucleotides]
        sense_c_seq = rc[SIRNA_SIZE - 5:SIRNA_SIZE - 1]
        antisense_five_seq = rc[0:pipe.end_nucleotides]
        antisense_c_seq = sirna_sequence[SIRNA_SIZE - 5:SIRNA_SIZE - 1]

        cases.append({
            "sirna_sequence": sirna_sequence,
            "sirna_sequence_n2": sirna_sequence_n2,
            "reverse_complement": str(rc),
            "sense": {
                "seq": str(sense_five_seq),
                "c_seq": str(sense_c_seq[::-1]),
                "shift": 1,
                "result": call(free_energy.calculate_free_energy, sense_five_seq,
                               check=True, strict=True, c_seq=sense_c_seq[::-1], shift=1),
            },
            "antisense": {
                "seq": str(antisense_five_seq),
                "c_seq": str(antisense_c_seq[::-1]),
                "shift": 1,
                "result": call(free_energy.calculate_free_energy, antisense_five_seq,
                               check=True, strict=True, c_seq=antisense_c_seq[::-1], shift=1),
            },
            # The two pipeline wrappers, so the pair is pinned as a unit.
            "free_energy_dangling_ends": call(
                pipe.free_energy_dangling_ends, sirna_sequence, sirna_sequence_n2),
            "free_energy3": call(pipe.free_energy3, sirna_sequence),
            "strand_selection": call(
                pipe.strand_selection, sirna_sequence, sirna_sequence_n2),
            "end_stability": call(
                pipe.end_stability, sirna_sequence, sirna_sequence_n2),
        })

    # The first two siRNAs take the no-dangling-ends path (sirna_sequence_n2 None).
    for i in (0, 1):
        cases.append({
            "sirna_sequence": sirnas[i][1],
            "sirna_sequence_n2": None,
            "free_energy3": call(pipe.free_energy3, sirnas[i][1]),
            "strand_selection": call(pipe.strand_selection, sirnas[i][1], None),
            "end_stability": call(pipe.end_stability, sirnas[i][1], None),
        })

    dump("thermo_dangling.json", {"description":
                                  "The pipeline's second call shape: "
                                  "calculate_free_energy(seq, c_seq=..., shift=1), plus the "
                                  "free_energy3 / free_energy_dangling_ends / strand_selection / "
                                  "end_stability wrappers. Floats are repr() strings.",
                                  "sirna_size": SIRNA_SIZE,
                                  "end_nucleotides": 3, "overhang": 2,
                                  "end_stability_treshold": repr(1.0),
                                  "cases": cases})


# --------------------------------------------------------------------------
# salt_correction.json -- methods 0-7 plus the out-of-range error
# --------------------------------------------------------------------------
def capture_salt_correction():
    seqs = [
        "GGGATGGCTCAAAGGCGTAGT",   # mixed
        "AAAAAAAAAAAAAAAAAAAAA",   # GC = 0
        "GCGCGCGCGCGCGCGCGCGCG",   # GC = 100
        "TTCTGGGAAGTGGTTTGCGCC",   # a real siRNA from the query
    ]
    ions = [
        {"Na": 20, "K": 50, "Tris": 0, "Mg": 0, "dNTPs": 0},   # siFi's own settings
        {"Na": 50, "K": 0, "Tris": 0, "Mg": 0, "dNTPs": 0},
        {"Na": 0, "K": 0, "Tris": 10, "Mg": 1.5, "dNTPs": 0.2},
        {"Na": 0, "K": 0, "Tris": 10, "Mg": 1.5, "dNTPs": 2.0},  # dNTPs >= Mg
        {"Na": 0, "K": 0, "Tris": 0, "Mg": 0, "dNTPs": 0},       # zero ions -> error
    ]
    cases = []
    for method in range(0, 9):
        for ion in ions:
            for seq in seqs:
                kw = dict(ion)
                kw["method"] = method
                kw["seq"] = seq
                kwargs = dict((k, repr(v)) for k, v in ion.items())
                kwargs["method"] = method
                kwargs["seq"] = seq
                cases.append({"kwargs": kwargs,
                              "result": call(free_energy.salt_correction, **kw)})
    # method 5/6/7 with no sequence at all
    for method in (5, 6, 7):
        cases.append({"kwargs": {"method": method, "seq": None, "Na": repr(20),
                                 "K": repr(50), "Tris": repr(0), "Mg": repr(0),
                                 "dNTPs": repr(0)},
                      "result": call(free_energy.salt_correction, Na=20, K=50,
                                     method=method, seq=None)})
    dump("salt_correction.json", {"description":
                                  "free_energy.salt_correction, methods 0-8 (8 is the error "
                                  "case) across four sequences and five ion combinations. "
                                  "Methods 6 and 7 are the SeqUtils.GC ones. "
                                  "Floats are repr() strings.",
                                  "cases": cases})


# --------------------------------------------------------------------------
# check_efficient.json -- exhaustive 2**3 inputs x 2**3 flags
# --------------------------------------------------------------------------
def capture_check_efficient():
    cases = []
    for sc in (False, True):
        for ec in (False, True):
            for ac in (False, True):
                pipe = make_pipeline(strand_check=sc, end_check=ec, accessibility_check=ac)
                for strand in (False, True):
                    for end in (False, True):
                        for access in (False, True):
                            cases.append({
                                "flags": {"strand_check": sc, "end_check": ec,
                                          "accessibility_check": ac},
                                "inputs": {"strand_selection": strand,
                                           "end_stability": end,
                                           "target_site_accessibility": access},
                                "result": call(pipe.check_efficient, strand, end, access),
                            })
    assert len(cases) == 64, len(cases)
    dump("check_efficient.json", {"description":
                                  "SifiPipeline.check_efficient, exhaustive: 2**3 inputs x "
                                  "2**3 check flags = 64 cases.",
                                  "cases": cases})


# --------------------------------------------------------------------------
# calculate_efficiency.json -- the full 5-tuple, including the terminal-rule branches
# --------------------------------------------------------------------------
def capture_calculate_efficiency(sirnas):
    # The terminal rule branches on sirna[sirna_size-3] (A/T vs not) and
    # sirna[1] (A/T, G/C, or neither). Pick real siRNAs covering every branch.
    wanted = {}
    for name, s in sirnas:
        s19 = "AT" if s[SIRNA_SIZE - 3] in "AT" else "other"
        s1 = "AT" if s[1] in "AT" else ("GC" if s[1] in "GC" else "other")
        wanted.setdefault((s19, s1), []).append((name, s))
    picks = []
    for key in sorted(wanted):
        picks.extend(wanted[key][:3])

    # Real accessibility values from the committed _lunp file, plus two synthetic
    # values that straddle the 0.1 threshold in both directions.
    lunp = load_lunp()
    lunp_values = [lunp(name) for name, _ in picks]

    cases = []
    for flags in [
        {},
        {"terminal_check": False},
        {"no_efficience": True},
        {"strand_check": False},
        {"end_check": False},
        {"accessibility_check": False},
        {"strand_check": False, "end_check": False, "accessibility_check": False},
        {"terminal_check": False, "accessibility_check": False},
    ]:
        pipe = make_pipeline(**flags)
        full = {"terminal_check": True, "no_efficience": False, "strand_check": True,
                "end_check": True, "accessibility_check": True}
        full.update(flags)
        for (name, sirna), real_lunp in zip(picks, lunp_values):
            idx = int(name[len("sirna"):]) - 1
            n2 = sirnas[idx - 2][1] if idx >= 2 else None
            for lunp_value in (real_lunp, 0.0, 1.0):
                cases.append({
                    "flags": full,
                    "sirna_name": name,
                    "sirna_sequence": sirna,
                    "sirna_sequence_n2": n2,
                    "lunp_data_xmer": repr(lunp_value),
                    "result": call(pipe.calculate_efficiency, sirna, n2, lunp_value),
                })
    dump("calculate_efficiency.json", {"description":
                                       "SifiPipeline.calculate_efficiency, returning the full "
                                       "5-tuple (is_efficient, strand_selection, end_stability, "
                                       "target_site_accessibility, thermo_effcicient). Real "
                                       "siRNAs chosen to cover every terminal-rule branch, each "
                                       "run at its real accessibility value and at 0.0/1.0. "
                                       "Floats are repr() strings.",
                                       "end_stability_treshold": repr(1.0),
                                       "target_site_accessibility_treshold": repr(0.1),
                                       "cases": cases})


# --------------------------------------------------------------------------
# sirna_windows.json -- create_sirnas on the real query
# --------------------------------------------------------------------------
def capture_sirna_windows(query_sequence):
    pipe = make_pipeline()
    sirna_file, tab_file = pipe.create_sirnas(query_sequence, SIRNA_SIZE)
    fasta = open(sirna_file).read()
    tab = open(tab_file).read()
    rows = [line.split("\t") for line in tab.splitlines()]
    dump("sirna_windows.json", {"description":
                                "SifiPipeline.create_sirnas on tests/data/query.fasta. "
                                "'fasta' and 'tab' are the two temp files' exact contents; note "
                                "the fasta headers are '> sirnaN' with a space after '>'.",
                                "query_length": len(query_sequence),
                                "sirna_size": SIRNA_SIZE,
                                "sirna_start_position": 0,
                                "count": len(rows),
                                "rows": rows,
                                "fasta": fasta,
                                "tab": tab})
    for path in (sirna_file, tab_file):
        os.remove(path)


# --------------------------------------------------------------------------
# bowtie_parse.json -- bowtie_to_lst on real bowtie 1.3.1 output
# --------------------------------------------------------------------------
def capture_bowtie_parse():
    pipe = make_pipeline()
    lines = open(BOWTIE_OUT).readlines()
    parsed = pipe.bowtie_to_lst(lines)
    shapes = {}
    for row in parsed:
        shapes[len(row)] = shapes.get(len(row), 0) + 1
    dump("bowtie_parse.json", {"description":
                               "SifiPipeline.bowtie_to_lst on tests/data/bowtie_output_mm2.txt "
                               "(real bowtie 1.3.1, -a -v 2 -y). bowtie always writes 8 columns "
                               "but leaves the mismatch descriptor empty on exact hits; the "
                               "leading .strip() then removes it, so those rows arrive as 7 "
                               "fields and get the integer 0 appended. Both shapes and both "
                               "strands are present.",
                               "input_file": "tests/data/bowtie_output_mm2.txt",
                               "input_lines": len(lines),
                               "field_count_histogram": dict(
                                   (str(k), v) for k, v in shapes.items()),
                               "rows": parsed})


# --------------------------------------------------------------------------
# analysis.json -- get_table_data / get_target_data
# --------------------------------------------------------------------------
def load_lunp():
    """Return name -> accessibility value, shaped exactly as data_to_json does."""
    lunp_data = np.loadtxt(LUNP, dtype="str")
    lunp_data = np.delete(lunp_data, np.r_[:SIRNA_SIZE - 1], 0)

    def get(sirna_name):
        idx = int(sirna_name.split("sirna")[1]) - 1
        return lunp_data[idx, :].astype(np.float).tolist()[ACCESSIBILITY_WINDOW]

    return get


def capture_analysis():
    """Build a real pipeline JSON with data_to_json, commit it, then pin the two
    aggregation functions against it."""
    pipe = make_pipeline(mode=0)
    query_sequence = read_query()
    pipe.len_seq = len(query_sequence)
    pipe.sirna_l = [[name, seq] for name, seq in sirna_table(query_sequence)]

    bowtie_lines = open(BOWTIE_OUT).readlines()
    pipe.bowtie_data_l = pipe.bowtie_to_lst(bowtie_lines)

    lunp_data = np.loadtxt(LUNP, dtype="str")
    lunp_data = np.delete(lunp_data, np.r_[:SIRNA_SIZE - 1], 0)

    # TUB3 is the query's own gene, so it is the main target; every paralogue hit
    # is an off-target. This is what the Qt dialog would have been used to pick.
    main_targets = ["NM_125665.4"]

    with quiet():
        json_lst = pipe.data_to_json(
            "TUB3_fragment", pipe.bowtie_data_l, False, lunp_data, main_targets)

    # One record per line: compact enough to commit, still diffable line-wise.
    pipeline_json = os.path.join(DATA, "pipeline_design.json")
    with open(pipeline_json, "w") as fh:
        fh.write("[\n")
        for i, rec in enumerate(json_lst):
            fh.write("  %s%s\n" % (json.dumps(rec, sort_keys=True),
                                   "," if i < len(json_lst) - 1 else ""))
        fh.write("]\n")

    table_data = general_helpers.get_table_data(pipeline_json)
    off_target_dict, main_target_dict, efficient_dict, main_hits_histo = \
        general_helpers.get_target_data(pipeline_json, SIRNA_SIZE)

    dump("analysis.json", {"description":
                           "general_helpers.get_table_data / get_target_data over "
                           "tests/data/pipeline_design.json, which is itself real "
                           "SifiPipeline.data_to_json output for design mode with "
                           "main_targets=['NM_125665.4']. Sets are serialised sorted; "
                           "table_data is sorted because Counter.most_common breaks ties "
                           "by dict order.",
                           "input_file": "tests/data/pipeline_design.json",
                           "sirna_size": SIRNA_SIZE,
                           "main_targets": main_targets,
                           "record_count": len(json_lst),
                           "table_data": sorted(table_data),
                           "off_target_dict": dict(
                               (k, sorted(v)) for k, v in off_target_dict.items()),
                           "main_target_dict": dict(
                               (k, sorted(v)) for k, v in main_target_dict.items()),
                           "efficient_dict": dict(
                               (k, list(v)) for k, v in efficient_dict.items()),
                           "main_hits_histo": main_hits_histo})
    print("wrote %-28s %8d bytes" % ("../data/pipeline_design.json",
                                     os.path.getsize(pipeline_json)))


def main():
    if not os.path.isdir(GOLDEN):
        os.makedirs(GOLDEN)
    query_sequence = read_query()
    sirnas = sirna_table(query_sequence)
    print("query %d bp -> %d siRNAs of %d nt" % (
        len(query_sequence), len(sirnas), SIRNA_SIZE))

    capture_thermo_nn(sirnas)
    capture_thermo_dangling(sirnas)
    capture_salt_correction()
    capture_check_efficient()
    capture_calculate_efficiency(sirnas)
    capture_sirna_windows(query_sequence)
    capture_bowtie_parse()
    capture_analysis()


if __name__ == "__main__":
    main()
