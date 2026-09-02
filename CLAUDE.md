# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

siFi — long dsRNA RNAi-target design and off-target prediction (Lück et al. 2019, doi:10.3389/fpls.2019.01023).
The original code is **Python 2 / PyQt4**, Windows-only, packaged with `py2exe`. Upstream is unmaintained.
This repo is a fork porting it to Python 3 with a CLI, an optional PyQt5 GUI, and tests.

## The port is planned — read `PLAN.md` first

`PLAN.md` is the source of truth for this rewrite. It is phased, and each phase is meant to be executed in a
fresh session: read `PLAN.md`'s Context, Status and the one target phase, do that phase, tick it off, commit,
then clear context. Do not start a phase before the one above it in Status is ticked.

## Goal of this fork

Port to **Python 3**, retain all existing functionality, and expose the analysis as a **CLI tool** (the GUI must
not be the only entry point; it is now ported too, as an optional extra, and is a client of the same pipeline). When choosing between preserving an odd upstream construct and writing idiomatic
Python 3, prefer Python 3 — but preserve *numerical behaviour* exactly; this is scientific code and the
thermodynamics must not silently change.

## Repository layout (post Phase 0)

```
legacy/            frozen Python 2 original — never edit, never import from src/
                   (12 modules + setup.py + Resources/ with the .ui/.qrc and generated Qt files)
src/sifi2/         the Python 3 port; see PLAN.md for the module map. Ported: thermo,
                   config (the SifiConfig dataclass), sirna, bowtie, rnaplfold, efficiency,
                   analysis, pipeline, cli, plots. The core is complete and the known defects are
                   fixed
src/sifi2/gui/     the PyQt5 GUI (Phase 7) — app, dialogs, workers, wizard, results, and
                   resources/ with the .ui/.qrc sources and their generated output. The ONLY
                   place in the package that may import Qt
tests/baseline/    the port's own end-to-end output for both modes, pinned by test_end_to_end.py;
                   see its README.md — these are not golden fixtures, they pin the port against itself
tests/golden/      JSON fixtures captured from legacy/ under Python 2.7 (Phase 1)
tests/data/        real inputs the fixtures are built from — query/reference FASTA, captured
                   bowtie and RNAplfold output, and a real data_to_json result
tests/py2_capture/ the capture script, the Qt stubs, and make_inputs.py; see its README.md
ToCopy/Images/     the two TIFF icons converted to PNG; the GUI's copies live in gui/resources/Images
environment.yml    curated sifi2 env spec; environment.lock.yml is the generated pin
environment-py2.yml throwaway sifi2-py2 env, for regenerating the fixtures only
```

`legacy/` is kept so the original stays importable under Python 2.7 for fixture regeneration. Do not lint,
format or "fix" anything in it.

## Python 2 constructs to fix, not imitate

Present throughout `legacy/`, so do not pattern-match on surrounding code when porting:
- `print x` statements (`free_energy.py`, and commented-out debug lines everywhere)
- `from types import *` + `assert type(x) is StringType` (`database_helpers.py`) — drop these asserts or use `isinstance`
- `dict.iteritems()` (`general_helpers.py`), `xrange` (`popup.py`)
- **Integer division**: `/` on ints changes meaning in Python 3. Check every division in `free_energy.py`,
  `sifi_pipeline.py` and `general_helpers.py` before assuming a port is behaviour-preserving.

## Known defects — fixed, and how to add to them

Every defect `PLAN.md` listed is fixed: the leaked `mkstemp` descriptors and the four unrestored `os.chdir()`
calls went in Phase 3, and Phase 6 did the five that move numbers (minus-strand hits dropped, design mode with no
hits producing nothing, the unreachable `DNA_TMM1` keys, the discarded `sirna.upper()`, and the per-target
position sets). The last one — `legacy/main.py:279`, where closing the mode dialog silently selected design
mode — went with the GUI in Phase 7, along with the Qt-only clean-ups Phase 5 handed over.

If another one turns up, fix it the same way: **one attributable commit per defect**, whose message says which
numbers moved and by how much. Never fold a behaviour change into an unrelated commit.

## External binaries

The pipeline shells out to `bowtie`, `bowtie-build` and `RNAplfold`. All three are in the `sifi2` env
(bowtie 1.3.1 — v1, `.ebwt` indices, *not* bowtie2; ViennaRNA 2.7.2). Resolve them with `shutil.which()`, never
by `os.chdir()`. The Windows `.exe` builds that shipped with upstream have been deleted from this fork.

## Qt and generated files

The GUI's Qt sources are `src/sifi2/gui/resources/{sifi2015.ui, db_wizard.ui, sifi_2015.qrc}`;
`ui_sifi2015.py`, `ui_db_wizard.py` and `sifi_2015_rc.py` beside them are **generated** — edit the sources and
regenerate, never hand-edit the output (the commands are in `resources/__init__.py`). `pyuic5` must be run with
**`--from-imports`**, or the generated file ends with a bare `import sifi_2015_rc` that resolves against
`sys.path`. `legacy/Resources/` holds the originals of all three and stays frozen.

**Qt is confined to `src/sifi2/gui/`.** No module in `src/sifi2/*.py` may import PyQt5, and none may import
`sifi2.gui` at module level — `sifi gui` imports it inside its handler. `import sifi2.pipeline` and
`import sifi2.cli` must keep working with no Qt installed; two tests in `tests/test_pipeline.py` enforce this.
PyQt5 is now in the `sifi2` env and in the `gui` extra, so it is present but optional.

## Environment

The project env is `sifi2` at `/mnt/apps/users/jnprice/conda/envs/sifi2` (Python 3.12). Run everything through
it: `conda run -n sifi2 <command>`. `sifi2` is installed into it editable, so `import sifi2` and the `sifi`
console script work without reinstalling after edits.

The throwaway `sifi2-py2` env (`environment-py2.yml`) runs `legacy/` to regenerate the golden fixtures. It
already exists; you only need it if `tests/golden/` has to be rebuilt. Note `repo.anaconda.com` is blocked on
this cluster, so any new env spec needs `- nodefaults` in its channel list.

## Linting

`conda run -n sifi2 ruff check src/ tests/` and `ruff format --check src/ tests/`. `ruff.toml` selects
`E,W,F,UP,B,SIM,I` and excludes `legacy/`, `ToCopy/`, `tests/py2_capture/` and the pyuic5/pyrcc5 output in
`src/sifi2/gui/resources/` — ruff cannot parse Python 2 syntax, `tests/py2_capture/` must keep running under
Python 2.7, and the generated Qt files are not edited by hand. Never point it at any of them.

## Testing

`conda run -n sifi2 pytest`. Tests live in `tests/`, with the golden fixtures in `tests/golden/`. The GUI
tests set `QT_QPA_PLATFORM=offscreen` themselves, so they need no display. Note `conda run` buffers all output
until the process exits, so a hanging test looks like silence — debug by calling the env's python directly. The fidelity
guarantee of this port is that the ported functions reproduce those fixtures exactly, so a failing golden test
means the port is wrong — never adjust a fixture to make a test pass unless that change is the deliberate
subject of the commit.

The fixtures record the original's behaviour **including its bugs**, and they keep doing so now that Phase 6 has
fixed those bugs in the port: `tests/golden/` is a capture of the original and stays byte-for-byte identical.
Where the port deliberately diverges, the *test* carries the divergence and says which defect it is; the
whole-pipeline diff lands in `tests/baseline/`, which is regenerated as part of the commit that moves it.
Regenerate the goldens with
`conda run -n sifi2-py2 python tests/py2_capture/capture.py`; it is deterministic, so a re-run must leave
`tests/golden/` byte-for-byte identical. See `tests/py2_capture/README.md`.
