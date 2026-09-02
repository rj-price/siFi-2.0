"""``python -m sifi2.gui`` — the same entry point as the ``sifi-gui`` script."""

from .app import main

if __name__ == "__main__":
    raise SystemExit(main())
