"""Command line entry point for siFi 2.0 (``PLAN.md`` Phase 4).

Two analysis subcommands mirroring the GUI's two modes, plus database
management::

    sifi design     --query in.fasta --db mydb --all-targets-main
    sifi offtarget  --query in.fasta --db mydb
    sifi db build   --fasta ref.fasta --name mydb
    sifi db list
    sifi db remove  mydb

Every parameter the GUI exposed is a flag here, with the GUI's own defaults —
including the six constants (``--winsize``, ``--span``, ``--temperature``,
``--sirna-start-position``, ``--overhang``, ``--end-nucleotides``) that were
buried in the pipeline body and reachable from neither the GUI nor a user.

Two mode-dependent defaults are taken from ``main.py``'s ``set_design_settings``
and ``set_offtarget_settings``: end-stability and accessibility checking are on
for design and off for off-target prediction. Both are overridable.

**Batch mode.** The original's ``for seq_record in SeqIO.parse(...)`` loop
``return``ed inside its first iteration, so every record after the first was
silently dropped — and the GUI refused multi-record input citing a batch mode
that was never written. Here every record is processed and gets its own set of
output files (JSON, TSV and PNG); ``--threads`` runs records concurrently, which is
safe now that the ``os.chdir()`` calls are gone.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from . import bowtie, rnaplfold
from .config import Mode, SifiConfig
from .pipeline import PipelineCancelled, QueryResult, SifiPipeline

__all__ = ["main"]

#: Columns of the per-query TSV, in the order ``data_to_json`` builds them.
TSV_COLUMNS = [
    "query_name",
    "sirna_name",
    "sirna_position",
    "sirna_sequence",
    "is_efficient",
    "strand_selection",
    "end_stability",
    "target_site_accessibility",
    "accessibility_value",
    "is_off_target",
    "hit_name",
    "reference_strand_pos",
    "strand",
    "mismatches",
]


class CliError(Exception):
    """A user-facing failure: reported as a message, not a traceback."""


# ----------------------------------------------------------------------
# Argument parsing
# ----------------------------------------------------------------------
def default_db_location() -> str:
    """Where databases live unless ``--db-location`` says otherwise.

    ``main.py`` derived this from ``QDesktopServices.storageLocation`` and then
    did ``.split('Local')[0] + '/Local/'`` to it, which is Windows-specific and
    produces nonsense on any path without a ``Local`` component.
    """
    from platformdirs import user_data_dir

    return os.path.join(user_data_dir("sifi2", appauthor=False), "databases")


def _add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Flags shared by ``design`` and ``offtarget``."""
    required = parser.add_argument_group("input and output")
    required.add_argument("--query", required=True, metavar="FASTA", help="query FASTA; every record is processed")
    required.add_argument("--db", required=True, metavar="NAME", help="name of the bowtie database to search")
    required.add_argument("--db-location", metavar="DIR", help=f"database directory (default: {default_db_location()})")
    required.add_argument("--outdir", default=".", metavar="DIR", help="output directory (default: the current one)")

    sirnas = parser.add_argument_group("siRNA windowing")
    sirnas.add_argument("--sirna-size", type=int, default=21, metavar="N", help="siRNA length (default: 21)")
    sirnas.add_argument(
        "--sirna-start-position", type=int, default=0, metavar="N", help="offset of the first window (default: 0)"
    )
    sirnas.add_argument(
        "--mismatches", type=int, default=0, metavar="N", help="mismatches allowed by bowtie (default: 0)"
    )

    rules = parser.add_argument_group("efficiency rules")
    rules.add_argument(
        "--strand", action=argparse.BooleanOptionalAction, default=True, help="strand-selection rule (default: on)"
    )
    rules.add_argument(
        "--terminal",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="terminal-nucleotide rule (default: on)",
    )
    rules.add_argument(
        "--end-stability",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="end-stability rule (default: on for design, off for off-target)",
    )
    rules.add_argument(
        "--end-stability-threshold", type=float, default=1.0, metavar="X", help="ΔΔG threshold (default: 1.0)"
    )
    rules.add_argument(
        "--accessibility",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="target-site accessibility rule (default: on for design, off for off-target)",
    )
    rules.add_argument(
        "--accessibility-threshold",
        type=float,
        default=0.1,
        metavar="X",
        help="unpaired-probability threshold (default: 0.1)",
    )
    rules.add_argument(
        "--accessibility-window", type=int, default=8, metavar="N", help="accessibility x-mer (default: 8)"
    )
    rules.add_argument(
        "--efficiency",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="score efficiency at all (default: on, but off when --mismatches > 0, as in the GUI)",
    )

    folding = parser.add_argument_group("RNAplfold (previously not user-visible)")
    folding.add_argument("--winsize", type=int, default=80, metavar="N", help="folding window (default: 80)")
    folding.add_argument("--span", type=int, default=40, metavar="N", help="maximum base-pair span (default: 40)")
    folding.add_argument("--temperature", type=float, default=22, metavar="C", help="folding temperature (default: 22)")

    energy = parser.add_argument_group("free energy (previously not user-visible)")
    energy.add_argument("--overhang", type=int, default=2, metavar="N", help="3' overhang length (default: 2)")
    energy.add_argument(
        "--end-nucleotides", type=int, default=3, metavar="N", help="nucleotides scored per end (default: 3)"
    )

    binaries = parser.add_argument_group("external binaries (default: found on PATH)")
    binaries.add_argument("--bowtie-path", metavar="PATH", help="bowtie binary, or the directory holding it")
    binaries.add_argument("--rnaplfold-path", metavar="PATH", help="RNAplfold binary, or the directory holding it")

    parser.add_argument(
        "--plot", action=argparse.BooleanOptionalAction, default=True, help="write <query>.png (default: on)"
    )
    parser.add_argument("--threads", type=int, default=1, metavar="N", help="query records to process at once")
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress progress messages")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sifi",
        description="siFi 2.0 — long dsRNA RNAi-target design and off-target prediction.",
        epilog="See https://doi.org/10.3389/fpls.2019.01023 for the method.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    design = subparsers.add_parser(
        "design",
        help="design an RNAi trigger against chosen main targets",
        description="RNAi design: score siRNAs against the main targets you choose, reporting the rest as off-targets.",
    )
    _add_common_arguments(design)
    targets = design.add_argument_group("main targets (exactly one is required)")
    exclusive = targets.add_mutually_exclusive_group(required=True)
    exclusive.add_argument(
        "--main-target", action="append", default=[], metavar="NAME", help="a main target; repeat for several"
    )
    exclusive.add_argument("--main-target-file", metavar="FILE", help="file of main target names, one per line")
    exclusive.add_argument("--all-targets-main", action="store_true", help="treat every hit as a main target")

    offtarget = subparsers.add_parser(
        "offtarget",
        help="predict off-targets of an existing RNAi trigger",
        description="Off-target prediction: report every database hit of every siRNA, with no main-target split.",
    )
    _add_common_arguments(offtarget)

    db = subparsers.add_parser("db", help="manage bowtie databases")
    db_subparsers = db.add_subparsers(dest="db_command", required=True, metavar="SUBCOMMAND")

    db_build = db_subparsers.add_parser("build", help="build a bowtie index from a FASTA file")
    db_build.add_argument("--fasta", required=True, metavar="FASTA", help="reference sequences to index")
    db_build.add_argument("--name", required=True, metavar="NAME", help="name to give the database")
    db_build.add_argument("--db-location", metavar="DIR", help=f"default: {default_db_location()}")
    db_build.add_argument("--bowtie-build-path", metavar="PATH", help="bowtie-build binary, or its directory")

    db_list = db_subparsers.add_parser("list", help="list the available databases")
    db_list.add_argument("--db-location", metavar="DIR", help=f"default: {default_db_location()}")

    db_remove = db_subparsers.add_parser("remove", help="delete databases")
    db_remove.add_argument("names", nargs="+", metavar="NAME", help="databases to delete")
    db_remove.add_argument("--db-location", metavar="DIR", help=f"default: {default_db_location()}")

    return parser


# ----------------------------------------------------------------------
# Wiring arguments to the pipeline
# ----------------------------------------------------------------------
def config_from_args(args: argparse.Namespace, mode: Mode) -> SifiConfig:
    """Build a :class:`SifiConfig`, resolving the mode-dependent defaults."""
    is_design = mode is Mode.DESIGN
    end_check = is_design if args.end_stability is None else args.end_stability
    accessibility_check = is_design if args.accessibility is None else args.accessibility
    # The GUI disabled every efficiency widget once mismatches > 0 and derived
    # no_efficience from that; the same rule applies here unless overridden.
    no_efficience = (args.mismatches > 0) if args.efficiency is None else not args.efficiency

    return SifiConfig(
        bowtie_db=args.db,
        db_location=args.db_location or default_db_location(),
        mode=mode,
        sirna_size=args.sirna_size,
        sirna_start_position=args.sirna_start_position,
        mismatches=args.mismatches,
        strand_check=args.strand,
        end_check=end_check,
        accessibility_check=accessibility_check,
        terminal_check=args.terminal,
        no_efficience=no_efficience,
        end_stability_threshold=args.end_stability_threshold,
        accessibility_threshold=args.accessibility_threshold,
        accessibility_window=args.accessibility_window,
        winsize=args.winsize,
        span=args.span,
        temperature=args.temperature,
        overhang=args.overhang,
        end_nucleotides=args.end_nucleotides,
        bowtie_path=args.bowtie_path,
        rnaplfold_path=args.rnaplfold_path,
    )


def read_main_target_file(path: str | os.PathLike[str]) -> list[str]:
    """Read main target names, one per line; blanks and ``#`` comments ignored."""
    try:
        lines = Path(path).read_text().splitlines()
    except OSError as error:
        raise CliError(f"cannot read --main-target-file {path}: {error}") from error
    names = [line.strip() for line in lines if line.strip() and not line.startswith("#")]
    if not names:
        raise CliError(f"--main-target-file {path} lists no target names")
    return names


def make_main_target_selector(args: argparse.Namespace):
    """A non-interactive stand-in for the GUI's modal ``popup.ListSelection``.

    The pipeline calls it with ``[(hit_name, hit_count)]`` and expects the names
    to treat as main targets. Names that matched nothing in the database are
    worth a warning: they are far more likely to be a typo than a deliberate
    request to designate every hit an off-target.
    """
    if args.all_targets_main:
        return lambda choices: [name for name, _count in choices]

    wanted = args.main_target or read_main_target_file(args.main_target_file)

    def selector(choices: list[tuple[str, int]]) -> list[str]:
        available = {name for name, _count in choices}
        for name in wanted:
            if name not in available:
                print(f"sifi: warning: main target {name!r} matched no hit in the database", file=sys.stderr)
        return list(wanted)

    return selector


# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------
def safe_stem(query_name: str) -> str:
    """A filename-safe stem for a FASTA record id.

    FASTA ids routinely carry ``/``, ``|`` and whitespace, none of which can go
    into a path unescaped.
    """
    stem = "".join(character if character.isalnum() or character in "-._" else "_" for character in query_name)
    return stem.strip("._") or "query"


def write_results(result: QueryResult, outdir: Path, stem: str, config: SifiConfig, plot: bool = True) -> list[Path]:
    """Write ``<query>.json``, ``<query>.tsv`` and ``<query>.png``; returns what was written."""
    json_path = outdir / f"{stem}.json"
    tsv_path = outdir / f"{stem}.tsv"

    with open(json_path, "w") as handle:
        json.dump(result.records, handle, indent=1)

    with open(tsv_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_COLUMNS)
        for record in result.records:
            writer.writerow(["" if record[column] is None else record[column] for column in TSV_COLUMNS])

    written = [json_path, tsv_path]
    if plot:
        # Imported here so `sifi db ...` and a failed run never pay for
        # matplotlib and seaborn, which together are most of the import cost.
        from . import plots

        png_path = outdir / f"{stem}.png"
        plots.save_plot(result, png_path, config.sirna_size, config.is_design)
        written.append(png_path)
    return written


# ----------------------------------------------------------------------
# Commands
# ----------------------------------------------------------------------
def run_analysis(args: argparse.Namespace, mode: Mode) -> int:
    config = config_from_args(args, mode)
    log = (lambda *a, **k: None) if args.quiet else print

    if args.threads < 1:
        raise CliError("--threads must be at least 1")
    if not os.path.exists(args.query):
        raise CliError(f"query file not found: {args.query}")
    if not bowtie.database_exists(config.bowtie_db, config.db_location):
        raise CliError(
            f"no bowtie database named {config.bowtie_db!r} in {config.db_location} "
            f"(build one with: sifi db build --fasta ref.fasta --name {config.bowtie_db})"
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    selector = make_main_target_selector(args) if mode is Mode.DESIGN else None
    pipeline = SifiPipeline(config, main_target_selector=selector)

    queries = read_queries(args.query)
    log(f"sifi: {len(queries)} query record(s) from {args.query}")

    if args.threads == 1:
        results = [(name, pipeline.run_query(name, sequence)) for name, sequence in queries]
    else:
        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            # `map` preserves input order, so the reports below stay in FASTA order.
            futures = pool.map(lambda query: pipeline.run_query(*query), queries)
            results = list(zip([name for name, _ in queries], futures, strict=True))

    return report_results(results, outdir, config, log, plot=args.plot)


def read_queries(query_file: str) -> list[tuple[str, str]]:
    """Every record of the query FASTA as ``(id, sequence)``."""
    from Bio import SeqIO

    queries = [(record.id, str(record.seq)) for record in SeqIO.parse(query_file, "fasta")]
    if not queries:
        raise CliError(f"no FASTA records found in {query_file}")

    # Records share an output filename if they share an id, so the second would
    # overwrite the first.
    duplicates = sorted(name for name, count in Counter(name for name, _ in queries).items() if count > 1)
    if duplicates:
        raise CliError(f"duplicate record ids in {query_file}: {', '.join(duplicates)}")
    return queries


def report_results(
    results: list[tuple[str, QueryResult]], outdir: Path, config: SifiConfig, log, plot: bool = True
) -> int:
    """Write each result and summarise; non-zero only if nothing was produced."""
    written_any = False
    for name, result in results:
        if result.message:
            # Design mode with no database hits still has efficiency records to
            # write, so the message is not by itself a reason to skip the query.
            log(f"sifi: {name}: {result.message}")
        if not result.has_hits:
            if not result.message:
                log(f"sifi: {name}: no records produced")
            continue
        paths = write_results(result, outdir, safe_stem(name), config, plot=plot)
        written_any = True
        efficient = sum(1 for record in result.records if record["is_efficient"])
        log(
            f"sifi: {name}: {len(result.records)} record(s), {efficient} efficient, "
            f"{len(result.table_data)} target(s) → {', '.join(path.name for path in paths)}"
        )
    if not written_any:
        log("sifi: no results were produced for any query record")
        return 1
    return 0


def run_db(args: argparse.Namespace) -> int:
    db_location = args.db_location or default_db_location()

    if args.db_command == "build":
        if not os.path.exists(args.fasta):
            raise CliError(f"reference FASTA not found: {args.fasta}")
        os.makedirs(db_location, exist_ok=True)
        prefix = bowtie.build_database(args.name, args.fasta, db_location, args.bowtie_build_path)
        print(f"sifi: built database {args.name!r} at {prefix}")
        return 0

    if args.db_command == "list":
        if not os.path.isdir(db_location):
            print(f"sifi: no database directory at {db_location}")
            return 1
        databases = bowtie.all_databases(db_location)
        if not databases["Database name"]:
            print(f"sifi: no databases in {db_location}")
            return 1
        print(f"{'Database name':<30} {'Size (MB)':>10}  Created")
        for name, size, created in zip(
            databases["Database name"], databases["Database size (MB)"], databases["Created"], strict=True
        ):
            print(f"{name:<30} {size:>10}  {created}")
        return 0

    deleted = bowtie.delete_databases(args.names, db_location)
    for name in args.names:
        print(f"sifi: {'deleted' if name in deleted else 'no such database:'} {name}")
    return 0 if len(deleted) == len(args.names) else 1


# ----------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(sys.argv[1:] if argv is None else argv)
    try:
        if args.command == "db":
            return run_db(args)
        return run_analysis(args, Mode.DESIGN if args.command == "design" else Mode.OFFTARGET)
    except (CliError, PipelineCancelled, bowtie.BowtieError, rnaplfold.RNAplfoldError) as error:
        print(f"sifi: error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
