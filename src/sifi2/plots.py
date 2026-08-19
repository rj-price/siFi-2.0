"""Rendering of the two siFi plots (``PLAN.md`` Phase 5).

Ported from ``legacy/show_plot.py``, which was a ``QMainWindow`` that built its
own figure, drew into it, saved a PNG to a leaked ``mkstemp`` path and carried a
print dialog, an export menu and a Qt navigation toolbar along with it. Here the
two rendering functions are **pure functions taking a matplotlib ``Figure``** —
the CLI hands them an Agg figure, and a future GUI can hand them the figure
behind a ``FigureCanvasQTAgg`` and get the same drawing.

Nothing in this module imports ``pyplot``: pyplot owns global figure state and
picks an interactive backend, neither of which a library should do. Figures are
built with :func:`design_figure` / :func:`offtarget_figure`, or supplied by the
caller.

Three fixes made while porting:

* ``figsize`` was **unbound for exactly 5 targets** — the chain of ``elif``\\ s
  ended ``len(hit_overview) > 5``, so five targets raised ``UnboundLocalError``
  before anything was drawn. Five now sizes with the many-targets case.
* ``ax1.xaxis.set_ticks(np.arange(0, query_length, 50))`` is taken back from
  ``legacy/create_plots.py``, the dead ancestor, where it was lost in the move
  to ``show_plot.py``. Without it matplotlib picks its own tick spacing and the
  x axis no longer reads as sequence position.
* ``the_table.set_fontsize(12)`` did nothing while automatic font sizing was
  still on, so the target table was drawn at whatever size fitted.

``show_plot.py``'s ``export_image``/``export_table`` are not ported: they are
file dialogs, and their ``filename += filename + '.png'`` bug (which doubles the
path rather than appending the extension) has no analogue here — the CLI names
its own output files.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np

from . import analysis

if TYPE_CHECKING:  # pragma: no cover - typing only
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from .pipeline import QueryResult

__all__ = [
    "apply_style",
    "design_figure",
    "offtarget_figsize",
    "offtarget_figure",
    "plot_design",
    "plot_offtarget",
    "save_plot",
]

#: Targets drawn in the off-target plot; the rest are summarised in the table.
MAX_PLOTTED_TARGETS = 5

#: X tick spacing, in bases.
TICK_SPACING = 50


def apply_style() -> None:
    """Apply the original's seaborn styling.

    ``sns.set(style="white", palette="muted", color_codes=True)`` is global rc
    state, which is why it is a separate call rather than something the plotting
    functions do to their caller behind its back. ``color_codes=True`` is not
    cosmetic here: it rebinds the single-letter ``'r'`` and ``'b'`` colours that
    both plots use to seaborn's muted versions.
    """
    import seaborn as sns

    sns.set_theme(style="white", palette="muted", color_codes=True)


def _histogram(positions, length: int) -> np.ndarray:
    """Counts per query position, padded to at least ``length`` entries."""
    return np.bincount(np.asarray(list(positions), dtype=np.intp), minlength=length)


def _format_coord(x, y) -> str:
    """The GUI's status-line readout; harmless and useful in an embedded canvas."""
    return f"Sequence position = {int(x)} , Nr. of siRNAs = {int(y)} "


# ----------------------------------------------------------------------
# Figure construction
# ----------------------------------------------------------------------
def _figure(figsize: tuple[float, float], dpi: float) -> Figure:
    from matplotlib.figure import Figure

    return Figure(figsize=figsize, dpi=dpi, facecolor="#FFFFFF")


def design_figure() -> Figure:
    """A figure sized for the design plot, which is always one panel."""
    return _figure((13, 4), 80)


def offtarget_figsize(target_count: int) -> tuple[tuple[float, float], float]:
    """``(figsize, dpi)`` for an off-target plot of ``target_count`` targets.

    One panel per target is stacked under the summary table, so the figure grows
    with the number of targets and then stops: past
    :data:`MAX_PLOTTED_TARGETS` only the table gets longer.

    The original left ``figsize`` **unbound at exactly 5** targets, raising
    ``UnboundLocalError``; five is handled by the many-targets branch here.
    """
    heights = {1: 5, 2: 6, 3: 8, 4: 10}
    if target_count in heights:
        return (13, heights[target_count]), 80
    if target_count < 1:
        return (13, 5), 80
    return (13, 20), 65


def offtarget_figure(target_count: int) -> Figure:
    """A figure sized for an off-target plot of ``target_count`` targets."""
    figsize, dpi = offtarget_figsize(target_count)
    return _figure(figsize, dpi)


# ----------------------------------------------------------------------
# The design plot
# ----------------------------------------------------------------------
def plot_design(figure: Figure, result: QueryResult, sirna_size: int) -> Axes:
    """Draw the RNAi design plot into ``figure`` and return its axes.

    Two histograms over the query — efficient siRNAs in red, main-target hits in
    blue — with main-target regions shaded green and off-target regions red.

    With no main targets chosen, only the efficiency histogram is drawn: that is
    the "design a construct against a sequence with no database hits" case the
    original explicitly provided for but could never reach, since ``data_to_json``
    discarded the records it needs (``PLAN.md`` Phase 6 defect 2, now fixed).
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    query_length = len(result.query_sequence)
    off_target_dict, main_target_dict, efficient_dict, main_hits_histo = analysis.get_target_data(
        result.records, sirna_size
    )

    if not result.main_targets:
        off_target_dict, main_target_dict, main_hits_histo = {}, {}, []

    off_target_positions: set[int] = set().union(*off_target_dict.values()) if off_target_dict else set()
    main_target_positions: set[int] = set().union(*main_target_dict.values()) if main_target_dict else set()

    efficient_positions: list[int] = []
    for positions in efficient_dict.values():
        efficient_positions.extend(positions)
    efficient_positions.sort()

    axes = figure.add_subplot(111)
    axes.format_coord = _format_coord

    efficient_histogram = _histogram(efficient_positions, query_length)
    axes.plot(efficient_histogram, "r-", label="Efficient siRNA hits")
    axes.plot(_histogram(main_hits_histo, query_length), "b-", label="Main target siRNA hits")

    # Shade the covered regions. Upstream drew one 1 bp patch per position;
    # contiguous runs are collapsed into a single patch instead, because
    # abutting translucent patches seam visibly where they overlap — the shading
    # came out striped, and a 500 bp query is 500 artists.
    for positions, colour in ((main_target_positions, "g"), (off_target_positions, "r")):
        for first, last in analysis.group_ranges(sorted(positions)):
            patch = Rectangle((first, 0.0), last - first + 1, sirna_size + 3, color=colour, alpha=0.2, linewidth=0)
            axes.add_patch(patch)
            patch.set_clip_on(False)

    max_y = max(int(efficient_histogram.max()), sirna_size) + 3
    axes.set_ylim([0, max_y])
    axes.set_xlim([0, query_length])
    axes.xaxis.set_ticks(np.arange(0, query_length, TICK_SPACING))
    axes.set_xlabel("mRNA sequence position\n\n", fontsize=12)
    axes.set_ylabel("Nr of siRNAs", color="b", fontsize=12)
    axes.set_title("RNAi design plot\n\n", fontsize=14)
    figure.tight_layout()

    handles = [
        Rectangle((0, 0), 1, 1, fc="g", alpha=0.2),
        Rectangle((0, 0), 1, 1, fc="r", alpha=0.2),
        Line2D([], [], color="r"),
        Line2D([], [], color="b"),
    ]
    axes.legend(
        handles,
        ["Main target", "Off target", "Efficient siRNAs", "All siRNAs"],
        bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
        loc=4,
        ncol=4,
        borderaxespad=0.5,
        frameon=True,
    )
    return axes


# ----------------------------------------------------------------------
# The off-target plot
# ----------------------------------------------------------------------
def offtarget_table_rows(table_data: list[list]) -> list[list]:
    """The table as drawn: the top targets, then a "N More targets" row."""
    if len(table_data) > MAX_PLOTTED_TARGETS:
        rows = [list(row) for row in table_data[:MAX_PLOTTED_TARGETS]]
        rows.append([f"{len(table_data) - MAX_PLOTTED_TARGETS} More targets", "...", "..."])
        return rows
    return [list(row) for row in table_data]


def plot_offtarget(figure: Figure, result: QueryResult, sirna_size: int) -> list[Axes]:
    """Draw the off-target plot into ``figure`` and return its axes, top first.

    A summary table of every target, then one histogram panel per target for the
    :data:`MAX_PLOTTED_TARGETS` most-hit ones: total siRNA hits in blue,
    efficient ones in red.
    """
    from matplotlib.lines import Line2D

    query_length = len(result.query_sequence)
    table_data = result.table_data or analysis.get_table_data(result.records)
    plotted = table_data[:MAX_PLOTTED_TARGETS]
    rows = len(plotted) + 1

    table_axes = figure.add_subplot(rows, 1, 1)
    table_axes.axis("off")
    table = table_axes.table(
        cellText=[[str(cell) for cell in row] for row in offtarget_table_rows(table_data)],
        colLabels=("Targets", "Total siRNA hits", "Efficient siRNA hits"),
        loc="center",
        cellLoc="left",
    )
    for cell in table.get_celld().values():
        cell.set_width(0.08)
    # Upstream set the font size with automatic sizing still on, so it was
    # recomputed at draw time and the call had no effect.
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(2.1, 2.1)

    all_axes = [table_axes]
    max_sirna_count = sirna_size

    for index, (hit_name, _total, _efficient) in enumerate(plotted, start=2):
        all_positions: list[int] = []
        efficient_positions: list[int] = []
        for record in result.records:
            if record["hit_name"] != hit_name:
                continue
            # 1-based, so the panels line up with the design plot's x axis.
            start = int(record["sirna_position"]) + 1
            positions = range(start, start + sirna_size)
            all_positions.extend(positions)
            if record["is_efficient"]:
                efficient_positions.extend(positions)

        axes = figure.add_subplot(rows, 1, index)
        axes.set_title(str(hit_name), loc="left")

        all_histogram = _histogram(all_positions, query_length)
        axes.plot(all_histogram, color="b", label="Total siRNA hits")
        axes.plot(_histogram(efficient_positions, query_length), color="r", label="Strand selected siRNAs")

        max_sirna_count = max(max_sirna_count, int(all_histogram.max()))
        axes.set_ylim([0, max_sirna_count + 3])
        axes.set_xlim([0, query_length + 5])
        axes.xaxis.set_ticks(np.arange(0, query_length, TICK_SPACING))
        axes.format_coord = _format_coord
        all_axes.append(axes)

    bottom = all_axes[-1]
    bottom.set_xlabel("\nRNAi trigger sequence position", fontsize=12)
    bottom.set_ylabel("siRNA counts per position", fontsize=12)

    figure.legend(
        [Line2D([], [], color="b"), Line2D([], [], color="r")],
        ["All siRNAs", "Efficient siRNAs"],
        loc="lower right",
        ncol=2,
        frameon=True,
    )
    # Upstream called ``subplots_adjust`` here too, but ``tight_layout`` runs
    # after it and recomputes every one of those values.
    figure.tight_layout()
    return all_axes


# ----------------------------------------------------------------------
# Saving
# ----------------------------------------------------------------------
def plot_result(figure: Figure, result: QueryResult, sirna_size: int, is_design: bool):
    """Draw whichever plot ``is_design`` calls for into ``figure``."""
    if is_design:
        return plot_design(figure, result, sirna_size)
    return plot_offtarget(figure, result, sirna_size)


def save_plot(
    result: QueryResult,
    path: str | os.PathLike[str],
    sirna_size: int,
    is_design: bool,
    style: bool = True,
) -> str:
    """Render ``result`` to a PNG at ``path`` and return the path.

    Uses the Agg canvas explicitly, so this works with no display and with no Qt
    installed — which the CLI needs and the headless import test guards.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    if style:
        apply_style()

    if is_design:
        figure = design_figure()
    else:
        table_data = result.table_data or analysis.get_table_data(result.records)
        figure = offtarget_figure(len(table_data))

    FigureCanvasAgg(figure)
    plot_result(figure, result, sirna_size, is_design)
    figure.savefig(os.fspath(path))
    return os.fspath(path)
