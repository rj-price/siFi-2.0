# siFi 2.0 — Python 2 → Python 3 rewrite with CLI and tests

> **This file lives at the repo root as `PLAN.md` and is the source of truth for the port.**
> It is written to be executed one phase at a time, in a fresh session per phase.

## How to use this plan

Each phase is self-contained. To start a phase in a clean context:

1. Read `CLAUDE.md` and this file's **Context**, **Status** and the target phase section — nothing else.
2. Do the phase.
3. Tick the phase in **Status** below, add a one-line note on anything that surprised you, and commit
   (`PLAN.md` update included in the same commit).
4. Clear context before the next phase.

Phases are strictly ordered — **Phase 1 must be complete before Phase 2**, since the golden fixtures are what
make every later phase verifiable. Phases 4 and 5 can be swapped. Phase 7 is deferred indefinitely.

## Status

- [ ] **Phase 0** — repo restructure, `sifi2` conda env, `pyproject.toml`
- [ ] **Phase 1** — golden fixtures captured from the Python 2 original ← *gate for everything below*
- [ ] **Phase 2** — `thermo.py` ported, `test_thermo.py` green
- [ ] **Phase 3** — rest of core ported, headless, Qt coupling broken
- [ ] **Phase 4** — CLI with batch mode
- [ ] **Phase 5** — plots + end-to-end baseline
- [ ] **Phase 6** — known defects fixed, one attributable commit each
- [ ] **Phase 7** — GUI (deferred)

**Notes from completed phases:** _(append here as you go)_

## Context

siFi (Lück et al. 2019, doi:10.3389/fpls.2019.01023) designs long dsRNA RNAi constructs and predicts
off-targets. The upstream code is Python 2.7 / PyQt4, Windows-only, packaged with `py2exe`, and unmaintained —
upstream's own README says it needs "a complete rewrite to meet current Python standards and library versions".

The goal is a Python 3 rewrite that **retains all functionality**, adds a **CLI** so it can run headlessly on the
HPC cluster, and is **covered by tests** — with the numerical output of the RNAi-efficiency and off-target
predictions proven identical to the original.

Three exploration passes established the state of play:

- The hand-written code is small: **~2,900 lines** across 12 modules (excluding 58k lines of generated Qt resource data).
- The science core **cannot currently be imported headlessly** — `sifi_pipeline.py` imports `popup` (PyQt4) and
  `show_plot` (PyQt4 + matplotlib), and calls a **modal Qt dialog in the middle of the pipeline** to pick main targets.
- Several dependencies have **hard breaking changes**, verified against the installed toolchain: `Bio.Alphabet`
  removed (Biopython ≥1.78), `SeqUtils.GC` removed and replaced by `gc_fraction` **which returns a fraction, not a
  percentage**, and `np.float` removed (NumPy ≥1.24).
- The original has **real defects**: minus-strand bowtie hits are silently dropped from the JSON, only the first
  FASTA record is ever processed, and two thermodynamic table keys are unreachable due to a leading-space typo.

Decisions taken: dedicated `sifi2` conda env; **CLI-first, GUI deferred**; **bit-for-bit fidelity** proven against
golden fixtures captured from the running Python 2 original; **port faithfully first, fix defects in separate
attributable commits**; **finish batch mode** so the CLI accepts multi-FASTA.

## Target repository shape

Remote is now `git@github.com:rj-price/siFi-2.0.git` (currently `origin` still points at `snowformatics/siFi21-`).

```
sifi-2.0/
├── environment.yml            # pinned sifi2 env
├── environment-py2.yml        # throwaway py2.7 reference env (fixture capture only)
├── pyproject.toml             # package metadata + `sifi` console entry point
├── ruff.toml                  # already written
├── LICENSE                    # CC BY-NC-SA 2.0, per upstream README
├── legacy/                    # frozen Python 2 original — never edited, never imported by src/
├── src/sifi2/
│   ├── thermo.py              # ← free_energy.py
│   ├── sirna.py               # ← create_sirnas / windowing / FASTA writing
│   ├── bowtie.py              # ← database_helpers.py + run_bowtie
│   ├── rnaplfold.py           # ← run_rnaplfold
│   ├── efficiency.py          # ← strand_selection / end_stability / check_efficient / calculate_efficiency
│   ├── pipeline.py            # ← sifi_pipeline.py, headless
│   ├── analysis.py            # ← general_helpers.py (table/target aggregation)
│   ├── plots.py               # ← show_plot.py, Agg rendering only
│   └── cli.py
└── tests/
    ├── golden/                # JSON fixtures captured from legacy/ under Python 2.7
    └── test_*.py
```

`legacy/` keeps the original importable under Python 2.7 for fixture regeneration. Move it with `git mv` so
history follows the files.

---

## Phase 0 — Repo and environment setup

1. Repoint `origin` to `git@github.com:rj-price/siFi-2.0.git`; keep upstream as a second remote for reference.
2. `git mv` the 12 Python modules + `Resources/` into `legacy/`. Delete outright: root `sifi_2015_rc.py` (0 bytes —
   see hazard below), `Resources/wizard_ui.py` (empty generated husk, no matching `.ui`), `Thumbs.db` ×3, `.idea/`,
   `imageformats/` and `ToCopy/{Bowtie,RNAplfold}` (Windows `.exe`/`.dll`). Keep `ToCopy/Images` — the icons are
   needed for the eventual GUI, but convert the `.tif` files to PNG so Qt needs no TIFF plugin.
3. Write `environment.yml` for the `sifi2` env and create it at `/mnt/apps/users/jnprice/conda/envs/sifi2`:
   `python=3.12, biopython, numpy, matplotlib, seaborn, pandas, bowtie=1.3.1, viennarna=2.7.2, pytest,
   pytest-cov, ruff`. (PyQt5 is deliberately **not** included until the GUI phase.) Pin exact builds after
   creation via `conda env export --no-builds`.
4. Add `pyproject.toml` with a `sifi = "sifi2.cli:main"` console script.

**Verify:** `conda run -n sifi2 python -c "import Bio, numpy, matplotlib"`, `conda run -n sifi2 bowtie --version`,
`conda run -n sifi2 RNAplfold --help`.

---

## Phase 1 — Capture golden fixtures from the Python 2 original

This must happen **before any porting**, and is the reason the whole plan is trustworthy.

`biopython=1.76` for `py27` is available on conda-forge (verified), so the original code can genuinely be run:

```yaml
# environment-py2.yml
dependencies: [python=2.7, biopython=1.76, numpy]
```

`legacy/free_energy.py` imports only `math`, `warnings` and four trivial Biopython touchpoints, so it is
importable under this env with no stubbing. `legacy/sifi_pipeline.py` needs `popup` and `show_plot` stubbed —
add `tests/py2_capture/` containing dummy `popup.py`/`show_plot.py` modules placed earlier on `sys.path`, which
lets every pure method be exercised without Qt.

Write `tests/py2_capture/capture.py`, run it under the py2 env, and commit its output to `tests/golden/`:

| Fixture | Source | Coverage |
|---|---|---|
| `thermo_nn.json` | `calculate_free_energy(seq)` | all 64 trinucleotides + ~50 real 21mers, `repr()` precision |
| `thermo_dangling.json` | `calculate_free_energy(seq, c_seq=..., shift=1)` | the pipeline's second call shape |
| `salt_correction.json` | `salt_correction(...)` | methods 1–7, incl. the two that use `SeqUtils.GC` |
| `check_efficient.json` | `SifiPipeline.check_efficient` | **exhaustive: 2³ inputs × 2³ flags = 64 cases** |
| `calculate_efficiency.json` | `calculate_efficiency` | full 5-tuple, incl. terminal-rule branches |
| `sirna_windows.json` | `create_sirnas` | windowing on a real query sequence |
| `bowtie_parse.json` | `bowtie_to_lst` | a captured real bowtie output, 7- and 8-field rows |
| `analysis.json` | `get_table_data` / `get_target_data` | derived from a checked-in JSON, sets serialised sorted |

Record floats at full `repr()` precision, never rounded. Also commit a small real query FASTA and a small
reference FASTA so an end-to-end index can be rebuilt reproducibly.

**Note on end-to-end:** a full original-vs-port JSON comparison is *not* achievable — the original's design mode
blocks on a Qt dialog and the bundled binaries are Windows `.exe`s. End-to-end output is therefore baselined on
the **port** in Phase 5, with the per-function goldens above carrying the fidelity guarantee.

---

## Phase 2 — Port the pure numerical core (`thermo.py`)

`free_energy.py` is a vendored, modified copy of `Bio.SeqUtils.MeltingTemp` that returns ΔG instead of Tm.
Port it first, in isolation, with `tests/test_thermo.py` green against `thermo_*.json` before anything else moves.

Preserve verbatim, these are the numerical contract:
- The seven lookup tables. **`DNA_TMM1` has two keys with a leading space** (`' CC/GC'`, `' GG/CA'`) which makes
  them unreachable by the exact-string lookup. Keep the spaces — "tidying" the whitespace changes the output.
- The siFi-specific default table selection: RNA nearest-neighbours (`RNA_NN3`) with **DNA** mismatch tables
  (`DNA_TMM1`, `DNA_IMM1`) and `RNA_DE2` dangling ends; `Na=20, K=50, saltcorr=5`. These differ from Biopython's
  defaults and a well-meaning "modernisation" would silently change them.
- `tao = 273.15 + 22` and `deltaG = (deltaH * 1000 - tao * deltaS) / 1000` — the single most load-bearing
  expression in the codebase.

Changes to make:
- **Drop Biopython entirely from this module.** `back_transcribe`, `complement` and the `GC` check are a few lines
  each; inlining them removes the `gc_fraction` percentage/fraction trap and makes `thermo.py` dependency-free
  and trivially testable. Where `SeqUtils.GC(seq) / 100` appears, the correct replacement is `gc_fraction(seq)`
  **with the `/ 100` removed** — getting this wrong is a silent 100× error in the salt correction.
- Delete the dead `Tm` computation (lines 424–429): computed, never returned.
- Delete `print seq` at line 284 — it fires on *every* energy computation.
- Give `_check` a proper error for an unrecognised `method` (currently `UnboundLocalError`).

---

## Phase 3 — Port the rest of the core, headless

Ported module by module, each with tests against its golden fixture.

**Breaking the Qt coupling** — three changes make a headless pipeline possible:
1. `get_main_target()`'s modal `popup.ListSelection` becomes an **injected callable**
   (`main_target_selector: Callable[[list[str]], list[str]]`). The CLI passes a function backed by
   `--main-target` / `--main-target-file` / `--all-targets-main`; the future GUI passes the dialog.
2. Plot construction moves out of the pipeline entirely. `pipeline.py` returns data; `plots.py` renders it.
3. `run_pipeline` **stops being a `@property`** — running an entire multi-minute pipeline on attribute access is
   the reason the GUI freezes. It becomes `run()`.

**Constructor** — replace the 17 positional args with a `@dataclass SifiConfig`. Include the six constants the GUI
never exposed (`winsize=80`, `span=40`, `temperature=22`, `sirna_start_position=0`, `overhang=2`,
`end_nucleotides=3`) as fields so the CLI can surface them. Note `mode` is `0 = design`, anything else = off-target,
and the Qt checkbox flags arrive as `0/2`, not `0/1` — normalise to `bool` at the boundary.

**Mechanical Python 3 fixes** (each is a hard failure, not a warning):

| Site | Problem | Fix |
|---|---|---|
| `general_helpers.py:15,176` | `Bio.Alphabet` removed | Delete the import; set `record.annotations["molecule_type"] = "DNA"` — modern `SeqIO.write(..., "genbank")` *requires* it |
| `general_helpers.py:162` | `lambda (i,x): i-x` | Python 3 **syntax error** → `lambda ix: ix[0] - ix[1]` |
| `general_helpers.py:163` | `map(...)` then indexed | `list(map(...))` |
| `general_helpers.py:139` | `range(...)` then `.extend()` | `list(range(...))` |
| `general_helpers.py:177,185` | `.iteritems()` | `.items()` |
| `sifi_pipeline.py:252` | `np.float` | `float` |
| `sifi_pipeline.py:352` | `prc.stdin.write(str)` | `Popen(..., text=True)` + `communicate(input=seq + "\n")` |
| `sifi_pipeline.py:153,222` | `print x` statements | remove (they are debug noise) |
| `database_helpers.py:6,38-41` | `from types import *`, `StringType` asserts | delete the asserts |
| `database_helpers.py:70` | `WindowsError` | **`NameError` on Linux at catch time** → `except OSError` |
| `database_helpers.py:20` | `(file_size/1000000)*3` | `//` to preserve the Py2 floor |
| `popup.py:44` | `xrange` | `range` (GUI phase) |
| all modules | implicit relative imports | explicit package imports |

**Subprocess and filesystem** — the four `os.chdir()` calls exist because Windows resolves executables from the
CWD; **POSIX does not**, so they are simultaneously useless and dangerous (none is ever restored). Remove all four
and resolve `bowtie`/`bowtie-build`/`RNAplfold` via `shutil.which()` with an explicit `--bowtie-path` override.
Add return-code checking: `run_bowtie` currently cannot distinguish "bowtie crashed" from "no hits", because
`mkstemp` has already created the output file it checks for.

Replace all seven `mkstemp` sites with one helper (`fd, path = mkstemp(); os.close(fd); return path`) — the current
idiom leaks **6 file descriptors per query sequence** plus a temp dir. Use `tempfile.TemporaryDirectory` for the
run scratch space so it is cleaned up.

**Resource-import hazard for the GUI phase:** the generated Qt files end with a bare `import sifi_2015_rc`. Under
Python 2's implicit relative imports that found `Resources/sifi_2015_rc.py` (3.8 MB); under Python 3 it resolves
against `sys.path` and finds the **0-byte root-level file**, importing successfully and registering nothing — so
every icon silently renders blank. Delete the root file (Phase 0) and regenerate with `pyrcc5`.

---

## Phase 4 — CLI

`sifi2/cli.py` using `argparse`, two subcommands mirroring the two modes, plus DB management:

```
sifi design     --query in.fasta --db mydb [--main-target NAME ...| --all-targets-main] ...
sifi offtarget  --query in.fasta --db mydb ...
sifi db build   --fasta ref.fasta --name mydb
sifi db list / sifi db remove
```

Every GUI parameter becomes a flag, with the GUI's own defaults: `--sirna-size 21`, `--mismatches 0`,
`--terminal/--no-terminal` (on), `--strand/--no-strand` (on), `--end-stability` (design on, off-target off) with
`--end-stability-threshold 1.0`, `--accessibility` (design on, off-target off) with
`--accessibility-threshold 0.1` and `--accessibility-window 8`. Expose the previously-hidden
`--winsize/--span/--temperature`. Outputs: `--outdir` with `<query>.json`, `<query>.tsv` and `<query>.png`.

**Batch mode:** the original's `for seq_record in SeqIO.parse(...)` loop `return`s inside the first iteration, so
multi-FASTA is silently truncated — and the GUI refuses >1 record citing a batch mode that was never written.
The CLI processes every record, emitting one result set per record. This is the main reason a CLI is worth having
here. Add `--threads` to parallelise across records (safe only once `os.chdir` is gone).

Replace `main.py`'s Windows path derivation (`storageLocation(DataLocation).split('Local')[0] + '/Local/'`, which
produces nonsense on Linux) with `platformdirs`, overridable by `--db-location`.

---

## Phase 5 — Plots and end-to-end baseline

Port `show_plot.py`'s two rendering functions into `plots.py` as **pure functions taking a matplotlib `Figure`** —
Agg for the CLI, embeddable in a canvas later. Do **not** port `create_plots.py` (dead ancestor) but do take its
`xaxis.set_ticks(np.arange(0, len, 50))`, lost in `show_plot.py`.

Fix while porting: `figsize` is **unbound when there are exactly 5 hits** (`> 5` should be `>= 5`) → `NameError`;
and `filename += filename + '.png'` doubles the path instead of appending the extension.

Then build a small bowtie index from the committed reference FASTA and run both modes end to end, committing the
resulting JSON as the end-to-end baseline. Consolidate the five identical copies of `show_info_message` and the
duplicated export/print/menu code shared by `imageviewer.py` and `show_plot.py`.

---

## Phase 6 — Fix the known defects, one attributable commit each

Only now, with goldens green. Each commit changes the golden file **deliberately** and says which numbers moved.

1. **Minus-strand hits dropped.** In `data_to_json`, the JSON-record construction is indented inside
   `if strand == '+':`, so every minus-strand bowtie hit is discarded. For off-target *prediction* this is a
   substantive scientific defect. Fix, and re-baseline.
2. **Design mode with no bowtie hits yields nothing.** The `no_target` path sets `strand = None`, which the same
   indentation then discards, so the efficiency-only plot the code explicitly provides for is dead. Fix.
3. **`DNA_TMM1` leading-space keys.** Decide explicitly: fixing them makes two terminal mismatches contribute
   energy for the first time, changing ΔG for affected sequences. Recommend fixing, in its own commit, with the
   golden diff attached — but this is a judgement call worth a second opinion, since it changes published-tool output.
4. Minor: `sirna.upper()` is a no-op (return value discarded) while the terminal rule compares against uppercase
   literals, so lowercase FASTA input silently fails every terminal test; `main_target_dict` aliases one
   accumulating set across all keys (`general_helpers.py:154`); `mode` selection returns design mode if the
   dialog is closed with the window X.

---

## Phase 7 — GUI (deferred, separate effort)

Out of scope for the initial rewrite; recorded so it is not lost. Regenerate the `.ui`/`.qrc` files with
`pyuic5`/`pyrcc5` rather than hand-porting. The migration surface: ~40 `QtGui.*` widget classes move to
`QtWidgets`; `QPrinter`/`QPrintDialog` move to `QtPrintSupport`; three old-style `SIGNAL`/`SLOT` connections
(`main.py:108`, `threads.py:20`, `db_wizard.py:51`) become new-style with a declared `pyqtSignal`;
`QDesktopServices.storageLocation` is **removed in Qt5** (5 call sites) → `QStandardPaths.writableLocation`;
`QFileDialog.getOpenFileName` now returns a **tuple**, breaking `.isNull()` at `main.py:311` and `db_wizard.py:79`;
`backend_qt4agg` → `backend_qtagg`; the `"gtk"` style no longer exists. Also decide whether GenBank export
(`imageviewer.py` → `general_helpers.create_gbk`) returns — it is currently unreachable dead code, so "retain all
functionality" does not strictly require it.

---

## Verification

- **Unit/golden:** `conda run -n sifi2 pytest -q` — every ported function green against `tests/golden/`.
  `test_thermo.py` is the load-bearing one; `check_efficient` is exhaustively pinned at all 64 combinations.
- **Lint:** `conda run -n sifi2 ruff check src/ tests/` and `ruff format --check`. `legacy/` stays excluded —
  ruff cannot parse Python 2.
- **End-to-end:** build a bowtie index from the committed reference FASTA, then
  `sifi offtarget --query tests/data/query.fasta --db test_db --outdir /tmp/out` and
  `sifi design --query tests/data/query.fasta --db test_db --all-targets-main --outdir /tmp/out`; diff the JSON
  against the Phase 5 baseline and eyeball the PNGs.
- **Batch:** a 3-record multi-FASTA produces 3 result sets — the check that the early-`return` bug is gone.
- **Headless:** `python -c "import sifi2.pipeline"` in an env with **no PyQt installed** must succeed. This is the
  regression test that the Qt coupling stayed broken.
- **Fixture regeneration:** `tests/py2_capture/capture.py` re-runnable under the py2 env to confirm goldens were
  not hand-edited.

## Risks

- **Silent numerical drift** is the main risk, and the Phase 1 goldens exist specifically to catch it. The
  highest-risk single change is `SeqUtils.GC / 100` → `gc_fraction` (100× if done naively); it sits in salt-correction
  methods 6–7 which siFi never reaches at `saltcorr=5`, so it is pinned but not load-bearing.
- **Biopython `Seq` no longer subclasses `str`** (≥1.78), and `free_energy_dangling_ends` chains
  `reverse_complement().strip()[...][::-1]`. Pin this under Python 2 before touching it.
- Bowtie 1.3.1 accepts the original's exact argv (`-a -v N -y <index> -f <reads> <out>`) and RNAplfold 2.7.2
  still has `-W/-L/-u/-T` — both verified, so the external-tool interface should port unchanged.
- Upstream is **CC BY-NC-SA 2.0**: the fork must stay non-commercial and share-alike, and carry attribution.
  There is no `LICENSE` file in the repo today, only a line in the README — add one.
- **Context rot across phases.** The exploration behind this plan was expensive; the file:line detail in each
  phase is deliberately included so no phase needs to re-derive it. If a phase turns out to need context this
  file does not carry, add it to the phase section rather than re-exploring next time.

## Phase-0 deliverable: update `CLAUDE.md`

`CLAUDE.md` currently describes the repo as it was found. Once Phase 0 restructures it, update it to point at
`src/sifi2/` and `legacy/`, name the `sifi2` env, and reference this plan — otherwise every future session will
start from a stale picture of the tree.
