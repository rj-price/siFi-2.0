# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

siFi — long dsRNA RNAi-target design and off-target prediction (Lück et al. 2019, doi:10.3389/fpls.2019.01023).
The original code is **Python 2 / PyQt4**, Windows-only, packaged with `py2exe`. Upstream is unmaintained.
This repo is a fork porting it to Python 3 with a CLI and tests.

## The port is planned — read `PLAN.md` first

`PLAN.md` is the source of truth for this rewrite. It is phased, and each phase is meant to be executed in a
fresh session: read `PLAN.md`'s Context, Status and the one target phase, do that phase, tick it off, commit,
then clear context. Do not start a phase before the one above it in Status is ticked.

## Goal of this fork

Port to **Python 3**, retain all existing functionality, and expose the analysis as a **CLI tool** (the GUI must
not be the only entry point). When choosing between preserving an odd upstream construct and writing idiomatic
Python 3, prefer Python 3 — but preserve *numerical behaviour* exactly; this is scientific code and the
thermodynamics must not silently change.

## Repository layout (post Phase 0)

```
legacy/            frozen Python 2 original — never edit, never import from src/
                   (12 modules + setup.py + Resources/ with the .ui/.qrc and generated Qt files)
src/sifi2/         the Python 3 port; see PLAN.md for the module map (thermo, sirna, bowtie,
                   rnaplfold, efficiency, pipeline, analysis, plots, cli)
tests/golden/      JSON fixtures captured from legacy/ under Python 2.7 (Phase 1)
tests/data/        real inputs the fixtures are built from — query/reference FASTA, captured
                   bowtie and RNAplfold output, and a real data_to_json result
tests/py2_capture/ the capture script, the Qt stubs, and make_inputs.py; see its README.md
ToCopy/Images/     icons kept for the eventual GUI, converted from TIFF to PNG
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

## Known defects worth fixing during the port

Fix these in **separate, attributable commits** (Phase 6), not silently while porting. `PLAN.md` Phase 6 has the
full list with file:line; the recurring ones are the leaked `mkstemp` file descriptors, the four unrestored
`os.chdir()` calls used to locate the binaries, `prc.stdin.write(seq)` needing a text-mode pipe, and the
minus-strand bowtie hits silently dropped from the JSON.

## External binaries

The pipeline shells out to `bowtie`, `bowtie-build` and `RNAplfold`. All three are in the `sifi2` env
(bowtie 1.3.1 — v1, `.ebwt` indices, *not* bowtie2; ViennaRNA 2.7.2). Resolve them with `shutil.which()`, never
by `os.chdir()`. The Windows `.exe` builds that shipped with upstream have been deleted from this fork.

## Qt and generated files

`legacy/Resources/ui_sifi2015.py`, `ui_db_wizard.py` and `sifi_2015_rc.py` are **generated** from `sifi2015.ui`,
`db_wizard.ui` and `sifi_2015.qrc` — edit the `.ui`/`.qrc` sources and regenerate, don't hand-edit the outputs.
Under PyQt5 the tools are `pyuic5` / `pyrcc5`. The GUI is deferred to Phase 7; PyQt5 is deliberately **not** in
the `sifi2` env, and `import sifi2.pipeline` must keep working without it.

## Environment

The project env is `sifi2` at `/mnt/apps/users/jnprice/conda/envs/sifi2` (Python 3.12). Run everything through
it: `conda run -n sifi2 <command>`. `sifi2` is installed into it editable, so `import sifi2` and the `sifi`
console script work without reinstalling after edits.

The throwaway `sifi2-py2` env (`environment-py2.yml`) runs `legacy/` to regenerate the golden fixtures. It
already exists; you only need it if `tests/golden/` has to be rebuilt. Note `repo.anaconda.com` is blocked on
this cluster, so any new env spec needs `- nodefaults` in its channel list.

## Linting

`conda run -n sifi2 ruff check src/ tests/` and `ruff format --check src/ tests/`. `ruff.toml` selects
`E,W,F,UP,B,SIM,I` and excludes `legacy/`, `ToCopy/` and `tests/py2_capture/` — ruff cannot parse Python 2
syntax, and `tests/py2_capture/` must keep running under Python 2.7, so never point it at either.

## Testing

`conda run -n sifi2 pytest`. Tests live in `tests/`, with the golden fixtures in `tests/golden/`. The fidelity
guarantee of this port is that the ported functions reproduce those fixtures exactly, so a failing golden test
means the port is wrong — never adjust a fixture to make a test pass unless that change is the deliberate
subject of the commit.

The fixtures record the original's behaviour **including its bugs** — that is deliberate. `PLAN.md` Phase 6
fixes them one attributable commit at a time, each showing up as a reviewed golden diff. Regenerate with
`conda run -n sifi2-py2 python tests/py2_capture/capture.py`; it is deterministic, so a re-run must leave
`tests/golden/` byte-for-byte identical. See `tests/py2_capture/README.md`.
