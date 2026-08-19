"""Phase 5: ``sifi2.plots`` — the two figures, rendered without Qt.

There is no golden fixture for a plot: the original could not be run to produce
one (it needed PyQt4 and a display), and a PNG is not a stable artefact to
compare against anyway. What is testable is the *structure* the original
prescribed — how many panels, what is on each axis, which regions are shaded —
plus the three defects fixed while porting, of which the unbound ``figsize`` at
exactly five targets is the one that made the plot raise instead of draw.
"""

import json

import pytest
from conftest import DATA
from matplotlib.figure import Figure

from sifi2 import plots
from sifi2.pipeline import QueryResult


@pytest.fixture(scope="module")
def records():
    """Real pipeline records: the captured Python 2 design run, 1370 of them."""
    with open(DATA / "pipeline_design.json") as handle:
        return json.load(handle)


@pytest.fixture(scope="module")
def query_sequence():
    from Bio import SeqIO

    return str(next(SeqIO.parse(str(DATA / "query.fasta"), "fasta")).seq)


@pytest.fixture
def result(records, query_sequence):
    from sifi2 import analysis

    return QueryResult(
        query_name="TUB3_fragment",
        query_sequence=query_sequence,
        records=records,
        table_data=analysis.get_table_data(records),
        main_targets=["NM_125665.4"],
    )


def make_result(records, query_sequence, main_targets=None):
    from sifi2 import analysis

    return QueryResult(
        query_name="q",
        query_sequence=query_sequence,
        records=records,
        table_data=analysis.get_table_data(records),
        main_targets=main_targets,
    )


# ----------------------------------------------------------------------
# Figure sizing — including the bug that made five targets raise
# ----------------------------------------------------------------------
@pytest.mark.parametrize(
    ("targets", "expected"),
    [(1, ((13, 5), 80)), (2, ((13, 6), 80)), (3, ((13, 8), 80)), (4, ((13, 10), 80)), (6, ((13, 20), 65))],
)
def test_offtarget_figsize_matches_the_original_for_the_sizes_it_handled(targets, expected):
    assert plots.offtarget_figsize(targets) == expected


def test_offtarget_figsize_is_defined_for_exactly_five_targets():
    """Upstream's ``elif`` chain jumped from ``== 4`` to ``> 5``, leaving
    ``figsize`` unbound at exactly five targets — an ``UnboundLocalError``
    before a single artist was drawn."""
    assert plots.offtarget_figsize(5) == ((13, 20), 65)


def test_design_figure_is_a_single_wide_panel():
    figure = plots.design_figure()
    assert figure.get_size_inches().tolist() == [13, 4]
    assert figure.get_dpi() == 80


# ----------------------------------------------------------------------
# The design plot
# ----------------------------------------------------------------------
def test_design_plot_draws_both_histograms_over_the_query(result):
    figure = Figure()
    axes = plots.plot_design(figure, result, 21)

    assert [line.get_label() for line in axes.get_lines()] == ["Efficient siRNA hits", "Main target siRNA hits"]
    assert axes.get_xlim() == (0, len(result.query_sequence))
    assert axes.get_title() == "RNAi design plot\n\n"


def test_design_plot_ticks_every_fifty_bases(result):
    """``xaxis.set_ticks(np.arange(0, length, 50))`` came from
    ``legacy/create_plots.py`` and was lost in the move to ``show_plot.py``."""
    axes = plots.plot_design(Figure(), result, 21)
    ticks = axes.get_xticks().tolist()
    assert ticks == list(range(0, len(result.query_sequence), 50))


def test_design_plot_y_axis_leaves_room_above_the_taller_of_the_two(result):
    axes = plots.plot_design(Figure(), result, 21)
    tallest = max(line.get_ydata().max() for line in axes.get_lines())
    assert axes.get_ylim() == (0, max(tallest, 21) + 3)


def test_design_plot_shades_main_targets_green_and_off_targets_red(result):
    axes = plots.plot_design(Figure(), result, 21)
    colours = {patch.get_facecolor()[:3] for patch in axes.patches}
    # Green and red, both at alpha 0.2, and nothing else.
    assert len(colours) == 2
    assert all(patch.get_alpha() == 0.2 for patch in axes.patches)


def test_design_plot_shades_contiguous_runs_as_one_patch(query_sequence):
    """Upstream drew one 1 bp rectangle per covered position; a run of 21
    consecutive positions is one patch here."""
    records = [
        {"sirna_position": 10, "sirna_name": "sirna10", "hit_name": "T", "is_efficient": False, "is_off_target": False}
    ]
    axes = plots.plot_design(Figure(), make_result(records, query_sequence, ["T"]), 21)
    assert len(axes.patches) == 1
    patch = axes.patches[0]
    assert (patch.get_x(), patch.get_width()) == (10, 21)


def test_design_plot_with_no_main_targets_shows_only_the_efficiency_histogram(result):
    """The original blanked all three target structures when nothing was chosen
    as a main target, so the plot falls back to "which siRNAs are efficient"."""
    result.main_targets = []
    axes = plots.plot_design(Figure(), result, 21)
    assert len(axes.patches) == 0
    main_histogram = axes.get_lines()[1].get_ydata()
    assert main_histogram.max() == 0


def test_design_plot_legend_names_all_four_series(result):
    axes = plots.plot_design(Figure(), result, 21)
    labels = [text.get_text() for text in axes.get_legend().get_texts()]
    assert labels == ["Main target", "Off target", "Efficient siRNAs", "All siRNAs"]


# ----------------------------------------------------------------------
# The off-target plot
# ----------------------------------------------------------------------
def test_offtarget_plot_is_a_table_plus_one_panel_per_target(result):
    figure = plots.offtarget_figure(len(result.table_data))
    axes = plots.plot_offtarget(figure, result, 21)

    targets = [row[0] for row in result.table_data]
    assert 1 <= len(targets) <= plots.MAX_PLOTTED_TARGETS
    assert len(axes) == len(targets) + 1
    assert axes[0].tables  # the summary table
    assert [ax.get_title(loc="left") for ax in axes[1:]] == targets


def test_offtarget_panels_plot_total_and_efficient_hits(result):
    axes = plots.plot_offtarget(plots.offtarget_figure(len(result.table_data)), result, 21)
    panel = axes[1]
    assert [line.get_label() for line in panel.get_lines()] == ["Total siRNA hits", "Strand selected siRNAs"]
    total, efficient = (line.get_ydata() for line in panel.get_lines())
    # Every efficient siRNA is also a hit, so one curve is bounded by the other.
    assert (efficient <= total).all()


def test_offtarget_axes_label_only_the_bottom_panel(result):
    axes = plots.plot_offtarget(plots.offtarget_figure(len(result.table_data)), result, 21)
    assert axes[-1].get_xlabel() == "\nRNAi trigger sequence position"
    assert axes[-1].get_ylabel() == "siRNA counts per position"


def test_offtarget_table_lists_every_target_when_there_are_five_or_fewer():
    table_data = [[f"target{index}", 10, 5] for index in range(5)]
    assert plots.offtarget_table_rows(table_data) == table_data


def test_offtarget_table_summarises_the_targets_it_cannot_draw():
    table_data = [[f"target{index}", 10, 5] for index in range(9)]
    rows = plots.offtarget_table_rows(table_data)
    assert len(rows) == plots.MAX_PLOTTED_TARGETS + 1
    assert rows[-1] == ["4 More targets", "...", "..."]


def test_offtarget_plot_draws_at_most_five_panels(query_sequence):
    records = [
        {
            "sirna_position": position,
            "sirna_name": f"sirna{position}",
            "hit_name": f"target{position}",
            "is_efficient": False,
            "is_off_target": True,
        }
        for position in range(1, 9)
    ]
    figure = plots.offtarget_figure(8)
    axes = plots.plot_offtarget(figure, make_result(records, query_sequence), 21)
    assert len(axes) == plots.MAX_PLOTTED_TARGETS + 1


# ----------------------------------------------------------------------
# Saving
# ----------------------------------------------------------------------
@pytest.mark.parametrize("is_design", [True, False])
def test_save_plot_writes_a_png_without_a_display(tmp_path, result, is_design):
    path = tmp_path / "plot.png"
    assert plots.save_plot(result, path, 21, is_design=is_design) == str(path)
    assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
