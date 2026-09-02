# siFi 2.0

Long double-stranded RNAi-target design and off-target prediction.

This is a **fork of [siFi21-](https://github.com/snowformatics/siFi21-)**, ported from Python 2 / PyQt4 to
Python 3 and exposed as a command-line tool. Upstream is Windows-only, GUI-only and unmaintained; it noted
that the software needed "a complete rewrite to meet current Python standards and library versions". This is
that rewrite.

The method is unchanged. The port preserves the original's numerical behaviour exactly — the thermodynamics,
efficiency rules and off-target logic all reproduce fixtures captured from the Python 2 original bit-for-bit.
See [Relationship to upstream](#relationship-to-upstream) for the handful of deliberate differences.

## What's different

- **Runs on Python 3.11+**, on Linux, macOS and Windows; no PyQt4, no `py2exe`, no bundled `.exe` binaries,
  **and no administrator rights** — nothing is installed or copied outside your own user directories.
- **A real CLI.** `sifi design`, `sifi offtarget` and `sifi db`, scriptable and usable on an HPC cluster.
- **Results are written to files you choose.** Upstream's advice was to hunt through the system temp
  directory for a `tmpXXXX.json` after each run; here `--outdir` gets you JSON, TSV and a plot per query
  record, named after the record.
- **Multi-FASTA batch mode.** Every record in the query file is processed (the original stopped after the
  first), optionally in parallel with `--threads`.
- **Parameters that were hard-coded are now flags** — RNAplfold's window, span and temperature, and the
  free-energy overhang and end-nucleotide counts.
- **A GUI that is optional.** The desktop interface is ported to PyQt5 (`sifi gui`), but it is a client of
  the same pipeline the CLI uses, and it installs separately — the core never imports Qt.
- **Tested.** 1314 tests, including golden fixtures captured from the Python 2 original and an end-to-end
  baseline for both modes.
- **Known bugs fixed**, one attributable commit each — minus-strand hits are no longer dropped, design mode
  with no hits no longer silently produces nothing, and three smaller defects. Each fix's commit message says
  which numbers moved.

## Installation

The pipeline shells out to `bowtie` (v1, *not* bowtie2), `bowtie-build` and `RNAplfold`, so conda is the
easiest route — the environment file pins all three alongside the Python dependencies:

```bash
conda env create -f environment.yml
conda activate sifi2
sifi --help
```

`environment.lock.yml` holds a fully pinned export if you need to reproduce an exact environment.

To install into an existing environment instead, make sure `bowtie`, `bowtie-build` and `RNAplfold` are on
`PATH` (or pass `--bowtie-path` / `--rnaplfold-path`), then:

```bash
pip install -e .          # CLI only
pip install -e '.[gui]'   # ...and the desktop GUI
```

## Quickstart

```bash
# 1. index your reference sequences
sifi db build --fasta reference.fasta --name mydb

# 2. design an RNAi trigger, naming the genes you actually want to silence
sifi design --query query.fasta --db mydb --main-target GENE1 --main-target GENE2 --outdir results/

# 3. or check an existing trigger for off-targets
sifi offtarget --query trigger.fasta --db mydb --outdir results/
```

Databases live in a platform-appropriate data directory by default (`~/.local/share/sifi2/databases` on
Linux); `--db-location` overrides it, and `sifi db list` / `sifi db remove` manage them.

## The GUI

```bash
sifi gui        # or: sifi-gui
```

The desktop interface is the original's, ported to PyQt5: paste or open a sequence, pick a database, set the
siRNA and efficiency options, press **Start**. The plot opens in its own window with matplotlib's zoom, and
exports to PNG, CSV and — in design mode — GenBank.

It needs the `gui` extra (`pip install 'sifi2[gui]'`, or `pyqt` in the conda environment), and it is strictly
optional: `import sifi2.pipeline` works in an environment with no Qt at all, which is what lets the CLI run on
a headless cluster node.

Three things differ from the upstream GUI:

- **No elevated privileges.** Upstream shipped Windows `bowtie`/`RNAplfold` executables inside its install
  directory and copied them into `%LOCALAPPDATA%` on every start-up; from a `Program Files` install that
  needs administrator rights. Here the binaries are found on `PATH`, databases live in your user data
  directory (the same one the CLI uses, so both see the same databases) and preferences live in `QSettings`.
- **Multi-record FASTA is accepted**, and each record gets its own plot in the results window. Upstream
  refused it, pointing at a batch mode that was never written.
- **Cancelling the "which mode?" question cancels**, rather than silently starting a design run.

## The two modes

**`design`** — you know which genes you want to silence. Every siRNA window of the query is scored for
efficiency against your chosen main targets, and hits to anything else are reported as off-targets. Exactly
one of `--main-target` (repeatable), `--main-target-file` (one name per line) or `--all-targets-main` is
required. A `--main-target` matching no hit warns and continues rather than failing, since it is nearly
always a typo and the consequence — everything reported as an off-target — is otherwise silent.

**`offtarget`** — you have a trigger sequence already and want to know what else it hits. No main targets,
and the end-stability and accessibility rules default off.

## Output

Per query record, in `--outdir`:

| file | contents |
| --- | --- |
| `<record>.json` | every siRNA/hit record, as a JSON list of objects |
| `<record>.tsv` | the same records as a tab-separated table |
| `<record>.png` | the design or off-target plot (`--no-plot` to skip) |

Both tabular formats carry the same fields:

```
query_name  sirna_name  sirna_position  sirna_sequence  is_efficient  strand_selection
end_stability  target_site_accessibility  accessibility_value  is_off_target  hit_name
reference_strand_pos  strand  mismatches
```

Record ids are sanitised for use as filenames. Duplicate ids in one query file are a hard error rather than a
silent overwrite.

## Options worth knowing

Run `sifi design --help` for the full list. The defaults reproduce the original GUI's:

- `--sirna-size` (21), `--sirna-start-position` (0), `--mismatches` (0).
- The four efficiency rules — `--strand`, `--terminal`, `--end-stability`, `--accessibility` — each with
  `--no-` counterparts, plus `--end-stability-threshold` (1.0), `--accessibility-threshold` (0.1) and
  `--accessibility-window` (8). `--end-stability` and `--accessibility` default on for `design` and off for
  `offtarget`.
- `--efficiency/--no-efficiency` turns scoring off entirely. It defaults on, but flips off automatically when
  `--mismatches > 0` — the GUI disabled all four efficiency widgets in that case, so this matches it.
- RNAplfold: `--winsize` (80), `--span` (40), `--temperature` (22). Free energy: `--overhang` (2),
  `--end-nucleotides` (3). None of these were reachable from the GUI.
- `--threads` processes several query records at once. Output stays in FASTA order.

Exit codes: `0` success, `1` for a handled error (reported as a message, never a traceback), `2` for a
command-line error.

## Development

```bash
conda run -n sifi2 pytest
conda run -n sifi2 ruff check src/ tests/
conda run -n sifi2 ruff format --check src/ tests/
```

The GUI tests run through Qt's `offscreen` platform plugin, so they need no display; they skip if PyQt5 is
not installed. The Qt `.ui` and `.qrc` sources live in `src/sifi2/gui/resources/` and are compiled with
`pyuic5 --from-imports` / `pyrcc5` — edit the sources, never the generated files.

`PLAN.md` is the source of truth for the port and records every phase and deliberate deviation. `CLAUDE.md`
carries the working conventions. Repository layout:

```
src/sifi2/        the Python 3 port
src/sifi2/gui/    the PyQt5 GUI — the only place Qt is imported
legacy/           the frozen Python 2 original — never edited, never imported from src/
tests/golden/     JSON fixtures captured from legacy/ under Python 2.7
tests/baseline/   the port's own end-to-end output for both modes
```

`tests/golden/` records the original's behaviour *including its bugs* and stays byte-for-byte identical;
where the port deliberately diverges, the test carries the divergence and names the defect.

## Relationship to upstream

Deliberate differences beyond the port itself, all covered by tests:

- Five behavioural bugs fixed (above), so results differ from upstream where those bugs applied — most
  visibly, minus-strand hits now appear in the output.
- One ordering change: a histogram that inherited Python 2's set-iteration order is now emitted ascending.
  Same values, stable across versions.
- Shaded regions in the design plot are drawn as contiguous patches rather than one rectangle per base,
  which removes the seams the overlapping patches produced.

## Licence and citation

Attribution-NonCommercial-ShareAlike 2.0 Generic (CC BY-NC-SA 2.0) — see `LICENSE`. Copyright © 2015–2019
Stefanie Lück and the siFi authors; Python 3 port, CLI and tests © 2026 Jordan Price.

If you use siFi, cite the original publication:

> Lück S, Kreszies T, Strickert M, Schweizer P, Douchkov D (2019). siRNA-Finder (si-Fi) Software for RNAi-Target
> Design and Off-Target Prediction. *Frontiers in Plant Science* 10:1023.
> [doi:10.3389/fpls.2019.01023](https://doi.org/10.3389/fpls.2019.01023)

Upstream's [quick help PDF](https://github.com/snowformatics/siFi21-/blob/7ef3e82003883156ecfa2cc73c5d36e93e02f9c2/Quick_help.pdf)
documents the method and the GUI it describes.
