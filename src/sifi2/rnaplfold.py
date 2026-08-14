"""RNAplfold: local pair probabilities for target-site accessibility.

Ported from ``SifiPipeline.run_rnaplfold``. As with bowtie, the binary is
resolved with ``shutil.which()`` rather than by chdir'ing into a bundled
directory, the exit status is checked, and the pipe is opened in text mode —
``prc.stdin.write(seq)`` on a bytes pipe is a hard failure under Python 3.
"""

from __future__ import annotations

import os
import subprocess

import numpy as np

from .bowtie import resolve_binary

__all__ = ["RNAplfoldError", "accessibility_value", "load_lunp", "run_rnaplfold"]


class RNAplfoldError(RuntimeError):
    """RNAplfold exited non-zero or wrote no ``_lunp`` file."""


def run_rnaplfold(
    query_name: str,
    query_sequence: str,
    workdir: str | os.PathLike[str],
    sirna_size: int = 21,
    winsize: int = 80,
    span: int = 40,
    temperature: float = 22,
    rnaplfold_path: str | None = None,
) -> np.ndarray:
    """Run RNAplfold in ``workdir`` and return its ``_lunp`` table.

    RNAplfold names its output after the FASTA header, so the sequence is fed in
    as FASTA on stdin and the result read back from ``<query_name>_lunp``.
    """
    rnaplfold = resolve_binary("RNAplfold", rnaplfold_path)
    argv = [
        rnaplfold,
        "-W",
        f"{winsize:d}",
        "-L",
        f"{span:d}",
        "-u",
        f"{sirna_size:d}",
        "-T",
        f"{temperature:.2f}",
    ]
    process = subprocess.run(
        argv,
        input=f">{query_name}\n{query_sequence}\n",
        capture_output=True,
        text=True,
        cwd=str(workdir),
        check=False,
    )
    if process.returncode != 0:
        raise RNAplfoldError(f"RNAplfold failed (exit {process.returncode}): {process.stderr.strip()}")

    lunp_file = os.path.join(str(workdir), f"{query_name}_lunp")
    if not os.path.exists(lunp_file):
        raise RNAplfoldError(f"RNAplfold wrote no accessibility file for query {query_name!r}")
    return load_lunp(lunp_file, sirna_size)


def load_lunp(path: str | os.PathLike[str], sirna_size: int = 21) -> np.ndarray:
    """Load an RNAplfold ``_lunp`` table, as the pipeline consumes it.

    Values stay strings — the original loaded with ``dtype='str'`` and converted
    per row — and the first ``sirna_size - 1`` rows are dropped because their
    windows are incomplete. Row ``i`` then corresponds to ``sirnaI`` with
    ``i = I - 1``; column ``u`` is the probability that ``u`` bases ending at
    that position are unpaired (column 0 is the position itself).
    """
    lunp_data = np.loadtxt(path, dtype="str")
    return np.delete(lunp_data, np.r_[: sirna_size - 1], 0)


def accessibility_value(lunp_data: np.ndarray, sirna_name: str, accessibility_window: int) -> float:
    """Unpaired probability for one siRNA at the chosen window size.

    ``np.float`` was removed in NumPy 1.24; the original's ``.astype(np.float)``
    is plain ``float`` here, which is what it always meant.
    """
    index = int(sirna_name.split("sirna")[1]) - 1
    return lunp_data[index, :].astype(float).tolist()[accessibility_window]
