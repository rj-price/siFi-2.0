"""The siFi 2.0 desktop GUI (``PLAN.md`` Phase 7).

Everything Qt lives under this subpackage, and nothing outside it imports Qt:
``sifi2.pipeline`` and friends must stay importable in an environment with no
PyQt installed, which is what made the CLI and the HPC batch runs possible in
the first place. The GUI is a *client* of that core — it builds a
:class:`sifi2.config.SifiConfig`, hands it to :class:`sifi2.pipeline.SifiPipeline`
on a worker thread, and renders the returned data with :mod:`sifi2.plots`.

PyQt5 is an optional dependency (``pip install sifi2[gui]``, or ``pyqt`` in the
conda env). Import this package only behind a try/except in code that must run
headless.

**No elevated privileges are needed, on any platform.** Upstream shipped
``bowtie`` and ``RNAplfold`` as Windows ``.exe``\\ s inside its install directory
and copied them, at every start-up, into ``%LOCALAPPDATA%``
(``general_helpers.copying_files``); combined with a py2exe installer that landed
in ``Program Files``, that is what made siFi an administrator-rights install.
Here the binaries are resolved on ``PATH`` with :func:`shutil.which` and nothing
is ever copied or written outside the user's own directories: databases go to
:func:`sifi2.cli.default_db_location` (a ``platformdirs`` user data directory)
and preferences to ``QSettings``. See :mod:`sifi2.gui.app`.
"""

from __future__ import annotations

__all__ = ["main"]


def main(argv: list[str] | None = None, db_location: str | None = None) -> int:
    """Entry point for the ``sifi-gui`` console script; see :mod:`sifi2.gui.app`."""
    from .app import main as _main

    return _main(argv, db_location=db_location)
