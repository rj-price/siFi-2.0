"""bowtie 1.x: index construction, alignment, output parsing and DB management.

Ported from ``legacy/database_helpers.py`` and ``SifiPipeline.run_bowtie``.

Two things changed structurally. The binaries are resolved with
``shutil.which()`` instead of ``os.chdir()`` into a bundled directory: the chdir
existed only because Windows resolves executables from the current working
directory, it was never restored, and on POSIX it is simply dangerous. And the
return code is now checked — the original could not tell "bowtie crashed" from
"no hits", because ``mkstemp`` had already created the output file it tested for.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time

__all__ = [
    "BowtieError",
    "all_databases",
    "bowtie_to_lst",
    "build_database",
    "database_exists",
    "delete_databases",
    "resolve_binary",
    "run_bowtie",
]

#: The six files ``bowtie-build`` writes for one index.
INDEX_EXTENSIONS = [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]


class BowtieError(RuntimeError):
    """bowtie or bowtie-build exited non-zero, or could not be found."""


def resolve_binary(name: str, path_hint: str | None = None) -> str:
    """Find an executable, honouring an explicit path or directory override."""
    if path_hint:
        candidate = os.path.join(path_hint, name) if os.path.isdir(path_hint) else path_hint
        resolved = shutil.which(candidate)
        if resolved is None:
            raise BowtieError(f"{name} not found at {path_hint!r}")
        return resolved
    resolved = shutil.which(name)
    if resolved is None:
        raise BowtieError(f"{name} not found on PATH; pass an explicit path")
    return resolved


def run_bowtie(
    sirna_fasta: str | os.PathLike[str],
    db_prefix: str | os.PathLike[str],
    mismatches: int,
    out_path: str | os.PathLike[str],
    bowtie_path: str | None = None,
) -> list[str]:
    """Align the siRNA multi-FASTA against a bowtie index and return its rows.

    The argv is the original's: ``-a`` (all alignments), ``-v N`` (end-to-end,
    at most N mismatches, no quality weighting), ``-y`` (try hard).
    """
    bowtie = resolve_binary("bowtie", bowtie_path)
    argv = [bowtie, "-a", "-v", str(mismatches), "-y", str(db_prefix), "-f", str(sirna_fasta), str(out_path)]
    process = subprocess.run(argv, capture_output=True, text=True, check=False)
    if process.returncode != 0:
        raise BowtieError(f"bowtie failed (exit {process.returncode}): {process.stderr.strip()}")
    if not os.path.exists(out_path):
        return []
    with open(out_path) as handle:
        return handle.readlines()


def bowtie_to_lst(bowtie_data: list[str]) -> list[list]:
    """Split bowtie's output rows into fields.

    bowtie 1.3.1 always writes eight columns but leaves the mismatch descriptor
    empty on an exact hit; the leading ``strip()`` then takes the trailing tab
    with it, so those rows arrive as seven fields and get the integer ``0``
    appended. Both shapes are in ``tests/golden/bowtie_parse.json``.
    """
    bowtie_lst = []
    for line in bowtie_data:
        fields = line.strip().split("\t")
        if len(fields) == 7:
            fields.append(0)
        bowtie_lst.append(fields)
    return bowtie_lst


def build_database(
    db_name: str,
    fasta_file: str | os.PathLike[str],
    db_location: str | os.PathLike[str],
    bowtie_build_path: str | None = None,
) -> str:
    """Build a bowtie index, returning its prefix. Raises on failure."""
    bowtie_build = resolve_binary("bowtie-build", bowtie_build_path)
    prefix = os.path.join(str(db_location), str(db_name))
    process = subprocess.run([bowtie_build, str(fasta_file), prefix], capture_output=True, text=True, check=False)
    if process.returncode != 0:
        raise BowtieError(f"bowtie-build failed (exit {process.returncode}): {process.stderr.strip()}")
    if not database_exists(db_name, db_location):
        raise BowtieError(f"bowtie-build reported success but wrote no index for {db_name!r}")
    return prefix


def database_exists(db_name: str, db_location: str | os.PathLike[str]) -> bool:
    """True when every file of the six-file index is present."""
    prefix = os.path.join(str(db_location), str(db_name))
    return all(os.path.exists(prefix + extension) for extension in INDEX_EXTENSIONS)


def all_databases(db_location: str | os.PathLike[str]) -> dict[str, list]:
    """Summarise the indices in ``db_location``: name, approximate size, date."""
    database_dict: dict[str, list] = {"Database name": [], "Database size (MB)": [], "Created": []}
    for db_file in sorted(os.listdir(db_location)):
        if not db_file.endswith(".ebwt"):
            continue
        name = db_file.split(".")[0]
        if name in database_dict["Database name"]:
            continue
        file_date, file_size = _size_and_date(os.path.join(str(db_location), db_file))
        database_dict["Database name"].append(name)
        database_dict["Created"].append(file_date)
        # The original's estimate: three times the size of the first index file.
        # `//` keeps the Python 2 floor division it was written against.
        database_dict["Database size (MB)"].append((file_size // 1000000) * 3)
    return database_dict


def _size_and_date(db_file: str) -> tuple[str, int]:
    """Creation date as ``DD.Mon.YYYY`` plus size in bytes."""
    parts = time.ctime(os.path.getctime(db_file)).split(" ")
    # time.ctime pads a single-digit day with a space, giving an empty field.
    parts = [part for part in parts if part]
    file_date = f"{parts[2]}.{parts[1]}.{parts[4]}"
    return file_date, os.stat(db_file).st_size


def delete_databases(db_list: list[str], db_location: str | os.PathLike[str]) -> list[str]:
    """Delete the named indices; returns the names actually removed.

    The original caught ``WindowsError``, which is not defined on Linux and so
    raised ``NameError`` at catch time — the ``except`` clause could never run.
    """
    deleted = []
    for db in db_list:
        prefix = os.path.join(str(db_location), str(db))
        removed_any = False
        for extension in INDEX_EXTENSIONS:
            try:
                os.remove(prefix + extension)
                removed_any = True
            except OSError:
                pass
        if removed_any:
            deleted.append(db)
    return deleted
