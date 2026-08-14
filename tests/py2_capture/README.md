# Golden-fixture capture (Phase 1)

`tests/golden/*.json` is the fidelity contract of this port: the ported Python 3
functions must reproduce these values exactly. They were captured by running
`legacy/` — unmodified — under Python 2.7, which is the only way to know what the
original actually computes.

## Regenerating

```bash
conda env create -f environment-py2.yml          # once; creates sifi2-py2
conda run -n sifi2-py2 python tests/py2_capture/capture.py
```

Re-running must leave `tests/golden/` **byte-for-byte identical**. That is the
check that the fixtures were never hand-edited to make a test pass:

```bash
md5sum tests/golden/*.json > /tmp/sums
conda run -n sifi2-py2 python tests/py2_capture/capture.py
md5sum -c /tmp/sums
```

`tests/test_golden_fixtures.py` checks the fixtures are well-formed and cover
what `PLAN.md` Phase 1 says they cover; it runs in the normal `sifi2` env.

## How it runs legacy/ headlessly

`legacy/sifi_pipeline.py` imports `popup` and `show_plot`, both PyQt4. This
directory contains stub modules of those names and puts itself *ahead* of
`legacy/` on `sys.path`, so the import succeeds. Only the import matters:
`capture.py` calls the pure methods directly and never touches `run_pipeline` or
`get_main_target`, the only users of the two Qt classes. Both stubs raise if
called, so that stays true.

`legacy/` itself goes on `sys.path` (rather than being imported as a package)
because its modules use Python 2 implicit relative imports.

## Conventions

These files are committed and diffed, so the capture is deterministic:

- **Floats are `repr()` strings**, never JSON numbers and never rounded. Python
  2.7 and Python 3 both use the shortest-round-trip repr, so the string is
  identical on both; `test_recorded_reprs_round_trip_under_python3` pins that.
- **Sets are serialised sorted**, and every dict is dumped with `sort_keys=True`.
  `table_data` is sorted because `Counter.most_common` breaks ties by dict order.
- **Exceptions are recorded** as `{"error": "<type>: <message>"}`. The port has
  to reproduce the failures too, not just the successes.

## Committed inputs (`tests/data/`)

| File | What it is |
|---|---|
| `query.fasta` | 500 bp of *Arabidopsis thaliana* TUB3 (`NM_125665.4`, bases 301–800) |
| `query_multi.fasta` | 3 records, for the Phase 4 batch-mode test |
| `reference.fasta` | the real β-tubulin paralogue family — TUB3, TUB2, TUB4, TUB8 ×2, plus TUB4 reverse-complemented |
| `bowtie_output_mm2.txt` | real bowtie 1.3.1 output, `-a -v 2 -y`, 1430 hits |
| `TUB3_fragment_lunp` | real RNAplfold 2.7.2 output, `-W 80 -L 40 -u 21 -T 22.00` |
| `pipeline_design.json` | real `data_to_json` output, design mode, `main_targets=['NM_125665.4']` |

The reference is a paralogue family so that a query drawn from one member
genuinely cross-hits the others — real main targets and real off-targets, which
is what siFi exists to tell apart. The reverse-complemented entry is there so the
fixtures contain minus-strand hits, which `PLAN.md` Phase 6 defect 1 is about.

`bowtie_output_mm2.txt` and `TUB3_fragment_lunp` are committed rather than
regenerated so `capture.py` needs neither bowtie nor ViennaRNA in the py2 env.
The scripts that produced them are recorded in `PLAN.md`'s Phase 1 notes; both
use the legacy argv verbatim.

## Fixtures

| File | Pins |
|---|---|
| `thermo_nn.json` | `calculate_free_energy(seq)` — all 64 trinucleotides, 50 real 21mers, and every 3mer `free_energy3` actually feeds it |
| `thermo_dangling.json` | the `c_seq=..., shift=1` call shape, plus `free_energy3` / `free_energy_dangling_ends` / `strand_selection` / `end_stability`, and the `Bio.Seq` slicing that derives the arguments |
| `salt_correction.json` | methods 0–8 × 4 sequences × 5 ion combinations, including the two `SeqUtils.GC` methods |
| `check_efficient.json` | exhaustive: 2³ inputs × 2³ flags = 64 |
| `calculate_efficiency.json` | the full 5-tuple, over real siRNAs covering all four terminal-rule leaves |
| `sirna_windows.json` | `create_sirnas` windowing, and the exact temp-file contents |
| `bowtie_parse.json` | `bowtie_to_lst`, both row shapes and both strands |
| `analysis.json` | `get_table_data` / `get_target_data` |

## Defects deliberately captured

The fixtures record the original's behaviour, **bugs included** — that is the
point. Phase 6 fixes them one attributable commit at a time, and each fix should
show up as a deliberate golden diff. Visible here:

- **Minus-strand hits dropped.** 1430 bowtie rows in, 1370 JSON records out; the
  60 minus-strand hits are discarded by the `if strand == '+':` indentation.
  `test_analysis_pins_the_minus_strand_defect` asserts the defect is present.
- **`salt_correction` method 7** raises `ValueError: math domain error` for some
  ion combinations rather than reporting anything useful.
