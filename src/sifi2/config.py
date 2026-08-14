"""Run configuration for the siFi pipeline.

The original passed 17 positional arguments to ``SifiPipeline.__init__`` and hid
six further constants in the body where neither the GUI nor a user could reach
them (``winsize``, ``span``, ``temperature``, ``sirna_start_position``,
``overhang``, ``end_nucleotides``). All of them are fields here, so the CLI can
surface any of them (``PLAN.md`` Phase 3/4).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

__all__ = ["Mode", "SifiConfig", "qt_flag"]


class Mode(StrEnum):
    """The original's ``mode``: ``0`` meant design, anything else off-target."""

    DESIGN = "design"
    OFFTARGET = "offtarget"

    @classmethod
    def from_legacy(cls, mode: int) -> Mode:
        return cls.DESIGN if mode == 0 else cls.OFFTARGET


def qt_flag(value: object) -> bool:
    """Normalise a Qt checkbox state to ``bool``.

    ``QCheckBox.checkState()`` yields ``0``/``2``, never ``0``/``1``, and the
    original stored that straight into ``strand_check`` and friends. Everything
    inside the port is a real ``bool``; conversion happens at the boundary.
    """
    return bool(value)


@dataclass
class SifiConfig:
    """Everything a pipeline run needs, apart from the query sequences."""

    # Database
    bowtie_db: str = ""
    db_location: str = ""

    # What to run
    mode: Mode = Mode.DESIGN

    # siRNA windowing
    sirna_size: int = 21
    sirna_start_position: int = 0
    mismatches: int = 0

    # Efficiency rules
    strand_check: bool = True
    end_check: bool = True
    accessibility_check: bool = True
    terminal_check: bool = True
    no_efficience: bool = False
    end_stability_threshold: float = 1.0
    accessibility_threshold: float = 0.1
    accessibility_window: int = 8

    # RNAplfold
    winsize: int = 80
    span: int = 40
    temperature: float = 22

    # Free energy end rules
    overhang: int = 2
    end_nucleotides: int = 3

    # Binaries; None means "resolve on PATH with shutil.which"
    bowtie_path: str | None = None
    bowtie_build_path: str | None = None
    rnaplfold_path: str | None = None

    @property
    def is_design(self) -> bool:
        return self.mode is Mode.DESIGN
