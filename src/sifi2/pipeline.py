"""The siFi pipeline, headless.

Ported from ``legacy/sifi_pipeline.py``. Three changes break the Qt coupling
(``PLAN.md`` Phase 3):

1. Picking the main targets was a modal ``popup.ListSelection`` called from the
   middle of the pipeline. It is now an injected callable — the CLI backs it with
   ``--main-target``/``--all-targets-main``, a GUI can still back it with the
   dialog.
2. Plotting has moved out entirely. :meth:`SifiPipeline.run` returns data;
   ``plots.py`` renders it (Phase 5).
3. ``run_pipeline`` was a ``@property``, so a multi-minute pipeline ran on
   attribute access — which is exactly why the GUI froze. It is now ``run()``.

Design mode:
    1. Split the query into siRNAs and write them as multi-FASTA.
    2. bowtie against the database; the user picks which hits are main targets.
    3. RNAplfold over the query for target-site accessibility.
    4. Score every siRNA: strand selection, end stability, accessibility.

Off-target mode is the same without step 2 — every hit is reported and no
efficiency verdict is attributed to a main target.

``import sifi2.pipeline`` must keep working in an environment with no PyQt
installed; that is a Phase 3 acceptance test.
"""

from __future__ import annotations

import os
import tempfile
from collections import Counter
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field

import numpy as np

from . import analysis, bowtie, rnaplfold, sirna
from .config import SifiConfig
from .efficiency import calculate_efficiency

__all__ = ["MainTargetSelector", "PipelineCancelled", "QueryResult", "SifiPipeline"]

#: Called with ``[(hit_name, hit_count)]`` sorted by count, descending; returns
#: the names to treat as main targets, or ``None`` to cancel the run.
MainTargetSelector = Callable[[list[tuple[str, int]]], "list[str] | None"]

NO_TARGETS_MESSAGE = (
    "No targets found. Please make sure that the query and/or database sequences are in correct orientation."
)


class PipelineCancelled(Exception):
    """The main-target selector cancelled the run."""


@dataclass
class QueryResult:
    """Everything one query sequence produced."""

    query_name: str
    query_sequence: str
    records: list[dict] = field(default_factory=list)
    table_data: list[list] = field(default_factory=list)
    main_targets: list[str] | None = None
    message: str | None = None

    @property
    def has_hits(self) -> bool:
        return bool(self.records)


class SifiPipeline:
    """Run the siFi analysis over one or many query sequences."""

    def __init__(self, config: SifiConfig, main_target_selector: MainTargetSelector | None = None):
        self.config = config
        self.main_target_selector = main_target_selector

    # ------------------------------------------------------------------
    # Driving
    # ------------------------------------------------------------------
    def run(self, query_file: str | os.PathLike[str]) -> list[QueryResult]:
        """Process every record in a (multi-)FASTA query file.

        The original ``return``ed inside the first iteration of this loop, so
        only the first record of a multi-FASTA was ever analysed.
        """
        return list(self.iter_run(query_file))

    def iter_run(self, query_file: str | os.PathLike[str]) -> Iterator[QueryResult]:
        """As :meth:`run`, yielding each result as it is finished."""
        from Bio import SeqIO  # local: keeps the module importable without Biopython

        for seq_record in SeqIO.parse(str(query_file), "fasta"):
            yield self.run_query(seq_record.id, str(seq_record.seq))

    def run_query(self, query_name: str, query_sequence: str) -> QueryResult:
        """Run the pipeline for a single query sequence."""
        config = self.config
        sirnas = sirna.create_sirnas(query_sequence, config.sirna_size, config.sirna_start_position)

        with tempfile.TemporaryDirectory(prefix="sifi2_") as workdir:
            sirna_fasta = os.path.join(workdir, "sirnas.fasta")
            sirna_tab = os.path.join(workdir, "sirnas.tab")
            sirna.write_multi_fasta(sirnas, sirna_fasta, sirna_tab)

            bowtie_rows = bowtie.run_bowtie(
                sirna_fasta,
                os.path.join(config.db_location, config.bowtie_db),
                config.mismatches,
                os.path.join(workdir, "bowtie.out"),
                config.bowtie_path,
            )
            bowtie_data = bowtie.bowtie_to_lst(bowtie_rows)

            lunp_data = rnaplfold.run_rnaplfold(
                query_name,
                query_sequence,
                workdir,
                sirna_size=config.sirna_size,
                winsize=config.winsize,
                span=config.span,
                temperature=config.temperature,
                rnaplfold_path=config.rnaplfold_path,
            )

        main_targets = None
        if config.is_design and bowtie_data:
            main_targets = self.get_main_targets(bowtie_data)

        if bowtie_data:
            input_data: Iterable = bowtie_data
            no_target = False
        elif config.is_design:
            # No hits: the efficiency plot is still worth having, so the siRNAs
            # themselves stand in for the bowtie rows.
            input_data = [[name, sequence] for name, sequence in sirnas]
            no_target = True
        else:
            # Off-target mode with no hits has nothing to report.
            return QueryResult(query_name, query_sequence, message=NO_TARGETS_MESSAGE)

        records = self.data_to_json(query_name, input_data, no_target, lunp_data, main_targets, sirnas)
        if not records:
            return QueryResult(query_name, query_sequence, main_targets=main_targets, message=NO_TARGETS_MESSAGE)
        return QueryResult(
            query_name,
            query_sequence,
            records=records,
            table_data=analysis.get_table_data(records),
            main_targets=main_targets,
        )

    # ------------------------------------------------------------------
    # Main targets
    # ------------------------------------------------------------------
    def get_main_targets(self, bowtie_data: list[list]) -> list[str]:
        """Ask the injected selector which hits are the intended targets.

        Everything it does not choose is reported as an off-target.
        """
        all_targets = Counter(row[2] for row in bowtie_data)
        choices = sorted(all_targets.items(), key=lambda item: item[1], reverse=True)
        if self.main_target_selector is None:
            raise PipelineCancelled("design mode needs a main_target_selector; pass one to SifiPipeline()")
        main_targets = self.main_target_selector(choices)
        if main_targets is None:
            raise PipelineCancelled("main target selection was cancelled")
        return list(main_targets)

    # ------------------------------------------------------------------
    # Record construction
    # ------------------------------------------------------------------
    def data_to_json(
        self,
        query_name: str,
        input_data: Iterable,
        no_target: bool,
        lunp_data: np.ndarray,
        main_targets: list[str] | None,
        sirnas: list[tuple[str, str]],
    ) -> list[dict]:
        """Turn bowtie rows (or bare siRNAs) into scored records.

        One record per bowtie hit, carrying the siRNA, its position on the query,
        the target it hit and its efficiency verdict.

        Note the ``strand == "+"`` guard: minus-strand hits produce no record at
        all, and neither does the no-hits design path, since it sets
        ``strand = None``. Both are ``PLAN.md`` Phase 6 defects 1 and 2, kept
        here so that fixing them is a reviewed golden diff rather than a silent
        change of published numbers.
        """
        config = self.config
        json_lst = []

        for data_split in input_data:
            if not no_target:
                sirna_name = data_split[0]
                strand = data_split[1]
                hit_name = data_split[2]
                # Target position from bowtie, 0-based offset.
                reference_strand_pos = int(data_split[3])
                # Position on the query sequence, starting at 1.
                query_position = int(sirna_name.split("sirna")[1])
                sirna_sequence = data_split[4]
                missmatches = data_split[7]
                off_target = (hit_name not in main_targets) if config.is_design else None
            else:
                # No bowtie hits: the siRNA file carries the efficiency data.
                sirna_name = data_split[0]
                query_position = int(sirna_name.split("sirna")[1])
                sirna_sequence = data_split[1]
                off_target = False
                strand = None
                hit_name = False
                reference_strand_pos = None
                missmatches = None

            if strand != "+":
                continue

            # The dangling-end partner is the siRNA two positions upstream, so
            # the first two siRNAs of a query have none.
            sirna_sequence_n2 = None if query_position in (1, 2) else sirnas[query_position - 3][1].strip()

            lunp_data_xmer = rnaplfold.accessibility_value(lunp_data, sirna_name, config.accessibility_window)

            (
                is_efficient,
                strand_selection,
                end_stability,
                target_site_accessibility,
                _thermo_efficient,
            ) = calculate_efficiency(config, sirna_sequence, sirna_sequence_n2, lunp_data_xmer)

            json_lst.append(
                {
                    "query_name": query_name,
                    "sirna_name": sirna_name,
                    "sirna_position": query_position,
                    "sirna_sequence": sirna_sequence,
                    "is_efficient": is_efficient,
                    "strand_selection": strand_selection,
                    "end_stability": end_stability,
                    "target_site_accessibility": target_site_accessibility,
                    "accessibility_value": lunp_data_xmer,
                    "is_off_target": off_target,
                    "hit_name": hit_name,
                    "reference_strand_pos": reference_strand_pos,
                    "strand": strand,
                    "mismatches": missmatches,
                }
            )

        return json_lst
