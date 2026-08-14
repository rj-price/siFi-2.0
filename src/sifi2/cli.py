"""Command line entry point for siFi 2.0.

Placeholder until Phase 4 of ``PLAN.md``; the console script exists from Phase 0
so the packaging and entry point can be verified early.
"""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    print(
        "sifi: the command line interface is not implemented yet (PLAN.md, Phase 4).",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
