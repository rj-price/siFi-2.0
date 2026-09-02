"""The GUI's small dialogs, and the one copy of ``show_info_message``.

Upstream had **five** identical copies of ``show_info_message`` — ``main.py:365``,
``show_plot.py:336``, ``imageviewer.py:123``, ``popup.py:57`` and
``db_wizard.py:84`` — one per widget class. ``PLAN.md`` Phase 5 asked for them to
be consolidated and could not, because every one of them is Qt; this is where
they land (:func:`show_info_message`).

Two upstream defects are fixed here rather than reproduced:

* :class:`ModeDialog` replaces ``main.py``'s ``self.mode = msg_box.exec_()``.
  ``exec_()`` returns ``0`` when a dialog is dismissed with the window's close
  button, and ``0`` was design mode, so cancelling the "which mode?" question
  silently started a design run (``PLAN.md`` Phase 6, handed to Phase 7). The
  dialog returns a :class:`~sifi2.config.Mode` or ``None``, and ``None`` means
  cancelled — a state the caller has to handle.
* :class:`TargetSelectionDialog` replaces ``popup.ListSelection``, whose header
  row (``"siRNA Hits\\tHit name"``) was a list item like any other, so the
  selection was read back with an ``index - 1`` offset that silently returned the
  wrong target if the header were ever removed. The header is a real header here.
"""

from __future__ import annotations

from PyQt5 import QtCore, QtWidgets

from ..config import Mode

__all__ = ["ModeDialog", "TargetSelectionDialog", "ask_main_targets", "ask_mode", "show_info_message"]


def show_info_message(parent, message: str, title: str = "Information") -> None:
    """Pop up an information box. The single copy of upstream's five."""
    QtWidgets.QMessageBox.information(parent, title, message)


class ModeDialog(QtWidgets.QDialog):
    """ "Which mode would you like to use?", with a cancel that means cancel."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Choose mode")
        self._mode: Mode | None = None

        label = QtWidgets.QLabel("Which mode would you like to use?")
        font = label.font()
        font.setPointSize(12)
        label.setFont(font)

        design = QtWidgets.QPushButton(" RNAi design ")
        design.setDefault(True)
        offtarget = QtWidgets.QPushButton(" Off-target prediction ")
        cancel = QtWidgets.QPushButton("Cancel")

        design.clicked.connect(lambda: self._choose(Mode.DESIGN))
        offtarget.clicked.connect(lambda: self._choose(Mode.OFFTARGET))
        cancel.clicked.connect(self.reject)

        buttons = QtWidgets.QHBoxLayout()
        for button in (design, offtarget, cancel):
            button.setFont(font)
            buttons.addWidget(button)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(label)
        layout.addLayout(buttons)

    def _choose(self, mode: Mode) -> None:
        self._mode = mode
        self.accept()

    def selected_mode(self) -> Mode | None:
        """The chosen mode, or ``None`` if the dialog was cancelled or closed."""
        return self._mode


def ask_mode(parent=None) -> Mode | None:
    """Show :class:`ModeDialog` modally; ``None`` when the user cancelled."""
    dialog = ModeDialog(parent)
    dialog.exec_()
    return dialog.selected_mode()


class TargetSelectionDialog(QtWidgets.QDialog):
    """Pick the main targets of a design run — the pipeline's injected selector.

    The pipeline calls its ``main_target_selector`` with ``[(hit_name, count)]``
    sorted by count; whatever is not chosen is reported as an off-target. With no
    hits at all the run still scores the siRNAs themselves, so "no targets found"
    is an outcome rather than an error (``PLAN.md`` Phase 6 defect 2).
    """

    def __init__(self, choices: list[tuple[str, int]], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select true target")
        self._choices = list(choices)

        layout = QtWidgets.QVBoxLayout(self)

        if self._choices:
            layout.addWidget(QtWidgets.QLabel("Tick the intended target(s) of the construct:"))
            self.tree = QtWidgets.QTreeWidget()
            self.tree.setHeaderLabels(["Hit name", "siRNA hits"])
            self.tree.setRootIsDecorated(False)
            for name, count in self._choices:
                item = QtWidgets.QTreeWidgetItem([str(name), str(count)])
                item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
                item.setCheckState(0, QtCore.Qt.Unchecked)
                self.tree.addTopLevelItem(item)
            self.tree.resizeColumnToContents(0)
            layout.addWidget(self.tree)

            select_all = QtWidgets.QPushButton("Select all")
            select_all.clicked.connect(self.select_all)
            layout.addWidget(select_all, alignment=QtCore.Qt.AlignLeft)
        else:
            self.tree = None
            label = QtWidgets.QLabel("Sorry, no targets found!\nPress OK to plot the efficiency of the query sequence.")
            label.setAlignment(QtCore.Qt.AlignCenter)
            layout.addWidget(label)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel, QtCore.Qt.Horizontal, self
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.resize(460, 350)

    def select_all(self) -> None:
        for index in range(self.tree.topLevelItemCount()):
            self.tree.topLevelItem(index).setCheckState(0, QtCore.Qt.Checked)

    def selected_targets(self) -> list[str]:
        """The ticked hit names, in the order they were offered."""
        if self.tree is None:
            return []
        return [
            self._choices[index][0]
            for index in range(self.tree.topLevelItemCount())
            if self.tree.topLevelItem(index).checkState(0) == QtCore.Qt.Checked
        ]


def ask_main_targets(choices: list[tuple[str, int]], parent=None) -> list[str] | None:
    """Run :class:`TargetSelectionDialog`; ``None`` cancels the pipeline run."""
    dialog = TargetSelectionDialog(choices, parent)
    if dialog.exec_() != QtWidgets.QDialog.Accepted:
        return None
    return dialog.selected_targets()
