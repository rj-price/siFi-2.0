"""The results window: an embedded matplotlib canvas, exports and printing.

This is upstream's ``show_plot.DrawPlot`` and ``imageviewer.ImageViewer`` merged
into one window. There were two because the pipeline saved its figure to a PNG
and a separate viewer re-opened that PNG — so the plot was rasterised before it
was ever shown, the "zoom" was image scaling, and both classes carried their own
near-identical export menu, print handler and ``show_info_message``. Here the
figure is drawn straight onto a :class:`FigureCanvasQTAgg`, so zooming is
matplotlib's own and the exports come from the live figure. ``PLAN.md`` Phase 5
listed that de-duplication and handed it here, since all of it is Qt.

Fixed while porting, from ``PLAN.md`` Phase 5 and Phase 7:

* ``filename += filename + '.png'`` — in ``export_image`` and ``export_table``
  in *both* upstream classes — doubles the path instead of appending the
  extension, so "results" was saved as "resultsresults.png". See
  :func:`_with_suffix`.
* ``QFileDialog.getSaveFileName``/``getOpenFileName`` return a ``(path, filter)``
  tuple in Qt5, not a ``QString``, so upstream's ``.isNull()`` check would raise
  ``AttributeError`` on every file dialog.
* ``QPrinter``/``QPrintDialog`` moved to ``QtPrintSupport`` in Qt5.
* The GenBank export was unreachable dead code upstream — ``imageviewer`` was
  never shown, because the call that would have shown it is commented out in
  ``main.py``. It is kept and reachable here, so "retain all functionality"
  holds in the direction that matters.

Batch runs put every query record in one window, chosen from a combo box: the
original refused more than one record, citing a batch mode nobody had written.
"""

from __future__ import annotations

import csv
import os

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PyQt5 import QtCore, QtPrintSupport, QtWidgets

from .. import analysis, plots
from ..pipeline import QueryResult
from .dialogs import show_info_message

__all__ = ["ResultsWindow", "SifiNavigationToolbar"]


class SifiNavigationToolbar(NavigationToolbar2QT):
    """The matplotlib toolbar, trimmed to the three tools upstream kept."""

    toolitems = [item for item in NavigationToolbar2QT.toolitems if item[0] in ("Home", "Zoom", "Save")]


def _with_suffix(path: str, suffix: str) -> str:
    """``path`` with ``suffix`` appended unless it is already there.

    Upstream wrote ``filename += filename + '.png'``, which is ``filename`` twice
    over followed by the extension.
    """
    return path if path.lower().endswith(suffix.lower()) else path + suffix


class ResultsWindow(QtWidgets.QMainWindow):
    """One window per run, showing each query record's plot in turn."""

    closed = QtCore.pyqtSignal()

    def __init__(self, results: list[QueryResult], sirna_size: int, is_design: bool, parent=None):
        super().__init__(parent)
        self.results = [result for result in results if result.has_hits]
        self.sirna_size = sirna_size
        self.is_design = is_design
        self._canvas: FigureCanvasQTAgg | None = None
        self._toolbar: SifiNavigationToolbar | None = None

        self.setWindowTitle("RNAi design plot" if is_design else "Off-target plot")
        self.resize(1100, 700)

        central = QtWidgets.QWidget(self)
        self._layout = QtWidgets.QVBoxLayout(central)
        self.setCentralWidget(central)

        self.selector = QtWidgets.QComboBox()
        self.selector.addItems([result.query_name for result in self.results])
        self.selector.currentIndexChanged.connect(self.show_result)
        # A single query needs no chooser; a batch does.
        chooser = QtWidgets.QHBoxLayout()
        chooser.addWidget(QtWidgets.QLabel("Query sequence:"))
        chooser.addWidget(self.selector, stretch=1)
        self._chooser_widget = QtWidgets.QWidget()
        self._chooser_widget.setLayout(chooser)
        self._chooser_widget.setVisible(len(self.results) > 1)
        self._layout.addWidget(self._chooser_widget)

        self._create_menus()
        if self.results:
            self.show_result(0)

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------
    @property
    def current(self) -> QueryResult:
        return self.results[max(self.selector.currentIndex(), 0)]

    def show_result(self, index: int) -> None:
        """Draw the ``index``-th query record, replacing whatever is shown."""
        result = self.results[index]

        if self._toolbar is not None:
            self.removeToolBar(self._toolbar)
            self._toolbar.deleteLater()
        if self._canvas is not None:
            self._layout.removeWidget(self._canvas)
            self._canvas.deleteLater()

        if self.is_design:
            figure = plots.design_figure()
        else:
            table_data = result.table_data or analysis.get_table_data(result.records)
            figure = plots.offtarget_figure(len(table_data))

        self._canvas = FigureCanvasQTAgg(figure)
        plots.plot_result(figure, result, self.sirna_size, self.is_design)
        self._toolbar = SifiNavigationToolbar(self._canvas, self)
        self.addToolBar(self._toolbar)
        self._layout.addWidget(self._canvas)
        self._canvas.draw_idle()
        self.statusBar().showMessage(
            f"{result.query_name}: {len(result.records)} siRNA record(s), "
            f"{sum(1 for record in result.records if record['is_efficient'])} efficient, "
            f"{len(result.table_data)} target(s)"
        )

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt's spelling
        self.closed.emit()
        event.accept()

    # ------------------------------------------------------------------
    # Menus
    # ------------------------------------------------------------------
    def _create_menus(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        file_menu.addAction("&Print…", self.print_result)
        file_menu.addSeparator()
        file_menu.addAction("E&xit", self.close)

        export_menu = self.menuBar().addMenu("&Save as")
        export_menu.addAction("Image file", self.export_image)
        export_menu.addAction("Table (CSV)", self.export_table)
        if self.is_design:
            # Main/off-target features only mean something once a main target has
            # been chosen, which is design mode only — as upstream had it.
            export_menu.addAction("GenBank file", self.export_gbk)

    def _save_path(self, caption: str, suffix: str) -> str | None:
        """Ask for a save location and give it ``suffix`` if it has none."""
        path, _selected_filter = QtWidgets.QFileDialog.getSaveFileName(
            self, caption, self._suggested_name(suffix), f"*{suffix}"
        )
        return _with_suffix(path, suffix) if path else None

    def _suggested_name(self, suffix: str) -> str:
        from ..cli import safe_stem

        return os.path.join(
            QtCore.QStandardPaths.writableLocation(QtCore.QStandardPaths.HomeLocation),
            safe_stem(self.current.query_name) + suffix,
        )

    # ------------------------------------------------------------------
    # Exports
    # ------------------------------------------------------------------
    def export_image(self) -> None:
        """Save the figure as it is drawn, at the figure's own resolution."""
        path = self._save_path("Save image", ".png")
        if not path:
            return
        try:
            self._canvas.figure.savefig(path)
        except OSError as error:
            show_info_message(self, f"Could not save image: {error}")
            return
        show_info_message(self, "Image successfully saved!")

    def export_table(self) -> None:
        """Save the per-target summary as CSV."""
        path = self._save_path("Export table", ".csv")
        if not path:
            return
        result = self.current
        table_data = result.table_data or analysis.get_table_data(result.records)
        try:
            with open(path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter=";", lineterminator="\n")
                writer.writerow(["Targets", "Total siRNA hits", "Efficient siRNA hits"])
                writer.writerows(table_data)
        except OSError as error:
            show_info_message(self, f"Could not save table: {error}")
            return
        show_info_message(self, "Table successfully saved!")

    def export_gbk(self) -> None:
        """Save the query with its main- and off-target hits as GenBank features."""
        path = self._save_path("Export GenBank", ".gbk")
        if not path:
            return
        result = self.current
        off_target_dict, main_target_dict, _efficient, _histogram = analysis.get_target_data(
            result.records, self.sirna_size
        )
        try:
            analysis.create_gbk(main_target_dict, off_target_dict, result.query_name, result.query_sequence, path)
        except OSError as error:
            show_info_message(self, f"Could not save GenBank file: {error}")
            return
        show_info_message(self, "GenBank file successfully saved!")

    # ------------------------------------------------------------------
    # Printing
    # ------------------------------------------------------------------
    def print_result(self) -> None:
        """Print the canvas, scaled to the page and keeping its aspect ratio."""
        printer = QtPrintSupport.QPrinter()
        if QtPrintSupport.QPrintDialog(printer, self).exec_() != QtWidgets.QDialog.Accepted:
            return
        from PyQt5 import QtGui

        pixmap = self._canvas.grab()
        painter = QtGui.QPainter(printer)
        rect = painter.viewport()
        size = pixmap.size()
        size.scale(rect.size(), QtCore.Qt.KeepAspectRatio)
        painter.setViewport(rect.x(), rect.y(), size.width(), size.height())
        painter.setWindow(pixmap.rect())
        painter.drawPixmap(0, 0, pixmap)
        painter.end()
