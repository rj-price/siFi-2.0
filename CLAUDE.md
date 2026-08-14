# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

siFi — long dsRNA RNAi-target design and off-target prediction (Lück et al. 2019, doi:10.3389/fpls.2019.01023).
Original code is **Python 2 / PyQt4**, Windows-only, packaged with `py2exe`. Upstream is unmaintained.

## The port is planned — read `PLAN.md` first

`PLAN.md` is the source of truth for this rewrite. It is phased, and each phase is meant to be executed in a
fresh session: read `PLAN.md`'s Context, Status and the one target phase, do that phase, tick it off, commit,
then clear context. Do not start a phase before the one above it in Status is ticked.

## Goal of this fork

Port to **Python 3**, retain all existing functionality, and expose the analysis as a **CLI tool** (the GUI must
not be the only entry point). When choosing between preserving an odd upstream construct and writing idiomatic
Python 3, prefer Python 3 — but preserve *numerical behaviour* exactly; this is scientific code and the
thermodynamics must not silently change.

## Python 2 constructs to fix, not imitate

Present throughout, so do not pattern-match on surrounding code:
- `print x` statements (`free_energy.py`, and commented-out debug lines everywhere)
- `from types import *` + `assert type(x) is StringType` (`database_helpers.py`) — drop these asserts or use `isinstance`
- `dict.iteritems()` (`general_helpers.py`), `xrange` (`popup.py`)
- **Integer division**: `/` on ints changes meaning in Python 3. Check every division in `free_energy.py`,
  `sifi_pipeline.py` and `general_helpers.py` before assuming a port is behaviour-preserving.

## Known defects worth fixing during the port

- `tempfile.mkstemp()` is called for its path only and the file descriptor is never closed
  (`sifi_pipeline.py`, several call sites) — leaks fds over a long run.
- `os.chdir()` is used to locate the bowtie / RNAplfold binaries before `subprocess.Popen`
  (`sifi_pipeline.run_bowtie`, `run_rnaplfold`, `database_helpers.create_bowtie_database`). Global process
  state; must go before anything runs concurrently or as a library. Use `cwd=` / absolute paths instead.
- `prc.stdin.write(seq)` in `run_rnaplfold` writes `str` — needs bytes or a text-mode pipe in Python 3.

## External binaries

The pipeline shells out to `bowtie`, `bowtie-build` and `RNAplfold`. `ToCopy/` contains **Windows `.exe`
builds** of these and is not usable on Linux — install Linux bowtie (v1, `.ebwt` indices — *not* bowtie2) and
ViennaRNA's RNAplfold instead. Never assume `ToCopy/` binaries can be executed here.

## Qt and generated files

`Resources/ui_sifi2015.py`, `Resources/ui_db_wizard.py`, `Resources/wizard_ui.py` and `Resources/sifi_2015_rc.py`
are **generated** from `sifi2015.ui`, `db_wizard.ui` and `sifi_2015.qrc` — edit the `.ui`/`.qrc` sources and
regenerate, don't hand-edit the outputs. Under PyQt5 the tools are `pyuic5` / `pyrcc5`.

## Environment

No conda env for this project exists yet. Envs live in `/mnt/apps/users/jnprice/conda/envs`; run tools with
`conda run -n <env> <command>`. Until an env is built, **nothing in this repo can be executed** — reason
statically and say so rather than claiming a change was verified by running it.

## Linting

`ruff.toml` is configured (`UP` pyupgrade rules on, generated Qt files excluded). Ruff **cannot parse Python 2
syntax**, so files that still contain `print x` will report a syntax error — that is expected and doubles as a
checklist of what is left to port. Ruff is not installed yet; add it to the project env.

## Testing

There are no tests. The port needs them: pin the current outputs of `free_energy.py` and the siRNA-generation /
scoring steps as fixtures before changing them, so the Python 3 version can be shown to match.
