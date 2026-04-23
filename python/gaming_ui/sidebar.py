from __future__ import annotations

from collections.abc import Callable
from typing import Any

from PySide6.QtCore import QAbstractListModel, QModelIndex, Qt
from PySide6.QtWidgets import QGridLayout, QGroupBox, QLabel, QLineEdit, QListView, QVBoxLayout, QWidget

from rxgaming_core import RxUnit, StructureSummary
from .units import format_value, label_for, from_display, UnitSystem


class UnitListModel(QAbstractListModel):
    def __init__(self, rx_units: list[RxUnit], unit_system: UnitSystem) -> None:
        super().__init__()
        self._rx_units = rx_units
        self._unit_system = unit_system

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._rx_units)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid():
            return None
        unit = self._rx_units[index.row()]
        if role == Qt.ItemDataRole.DisplayRole:
            return unit.name
        if role == Qt.ItemDataRole.ToolTipRole:
            return (
                f"{unit.name}\n"
                f"{label_for('tph', self._unit_system)}: {format_value('tph', unit.currentStructure.tph, self._unit_system)}\n"
                f"{label_for('ba', self._unit_system)}: {format_value('ba', unit.currentStructure.ba, self._unit_system)}\n"
                f"MCS: {format_value('mcs', unit.currentStructure.mcs, self._unit_system)}\n"
                f"CC: {format_value('cc', unit.currentStructure.cc, self._unit_system)}"
            )
        return None

    def refresh(self) -> None:
        self.layoutChanged.emit()


class StructureInfo(QWidget):
    def __init__(self, unit_system: UnitSystem) -> None:
        super().__init__()
        self._unit: RxUnit | None = None
        self._updating = False
        self._targets_changed_callback: Callable[[], None] | None = None
        self._unit_system = unit_system

        self.name_label = QLabel("No unit selected")
        self.result_label = QLabel("Result: -")

        self.current_header = QLabel("Current")
        self.target_header = QLabel("Target")
        self.treated_header = QLabel("Treated")

        self.current_tph = QLabel("-")
        self.current_ba = QLabel("-")
        self.current_mcs = QLabel("-")
        self.current_cc = QLabel("-")

        self.target_tph = QLineEdit()
        self.target_ba = QLineEdit()
        self.target_mcs = QLineEdit()
        self.target_cc = QLineEdit()
        self.target_edits = [self.target_tph, self.target_ba, self.target_mcs, self.target_cc]
        for edit in self.target_edits:
            edit.setMinimumWidth(84)

        self.treated_tph = QLabel("-")
        self.treated_ba = QLabel("-")
        self.treated_mcs = QLabel("-")
        self.treated_cc = QLabel("-")
        self.treated_values = [self.treated_tph, self.treated_ba, self.treated_mcs, self.treated_cc]

        self.treated_hint = QLabel("Run a treatment to populate treated values.")
        self.treated_hint.setWordWrap(True)

        layout = QGridLayout()
        layout.addWidget(self.name_label, 0, 0, 1, 4)
        layout.addWidget(QLabel("Measure"), 1, 0)
        layout.addWidget(self.current_header, 1, 1)
        layout.addWidget(self.target_header, 1, 2)
        layout.addWidget(self.treated_header, 1, 3)

        tph_label = label_for("tph", unit_system)
        metrics = [
            (tph_label, self.current_tph, self.target_tph, self.treated_tph),
            ("BA", self.current_ba, self.target_ba, self.treated_ba),
            ("MCS", self.current_mcs, self.target_mcs, self.treated_mcs),
            ("CC", self.current_cc, self.target_cc, self.treated_cc),
        ]
        for row, (label, current, target, treated) in enumerate(metrics, start=2):
            layout.addWidget(QLabel(label), row, 0)
            layout.addWidget(current, row, 1)
            layout.addWidget(target, row, 2)
            layout.addWidget(treated, row, 3)

        outer_layout = QVBoxLayout()
        outer_layout.addLayout(layout)
        outer_layout.addWidget(self.treated_hint)
        outer_layout.addWidget(self.result_label)
        outer_layout.addStretch(1)
        self.setLayout(outer_layout)

        for index, edit in enumerate(self.target_edits):
            edit.editingFinished.connect(lambda idx=index: self._apply_target_change(idx))

        self._set_treated_visual_state(False)

    def set_targets_changed_callback(self, callback: Callable[[], None]) -> None:
        self._targets_changed_callback = callback

    def update_for_unit(self, unit: RxUnit) -> None:
        self._unit = unit
        self._updating = True
        self.name_label.setText(f"Unit: {unit.name}")
        self.result_label.setText(f"Result: {unit.result.name}")

        self.current_tph.setText(format_value("tph", unit.currentStructure.tph, self._unit_system))
        self.current_ba.setText(format_value("ba", unit.currentStructure.ba, self._unit_system))
        self.current_mcs.setText(format_value("mcs", unit.currentStructure.mcs, self._unit_system))
        self.current_cc.setText(format_value("cc", unit.currentStructure.cc, self._unit_system))

        self.target_tph.setText(format_value("tph", unit.targetStructure.tph, self._unit_system))
        self.target_ba.setText(format_value("ba", unit.targetStructure.ba, self._unit_system))
        self.target_mcs.setText(format_value("mcs", unit.targetStructure.mcs, self._unit_system))
        self.target_cc.setText(format_value("cc", unit.targetStructure.cc, self._unit_system))

        self._update_treated_column(unit.treatedStructure)
        self._updating = False

    def _apply_target_change(self, metric_index: int) -> None:
        if self._updating or self._unit is None:
            return
        edit = self.target_edits[metric_index]
        try:
            value = float(edit.text())
        except ValueError:
            self._restore_target_value(metric_index)
            return

        attributes = ["tph", "ba", "mcs", "cc"]
        converted_value = from_display(attributes[metric_index], value, self._unit_system)
        setattr(self._unit.targetStructure, attributes[metric_index], converted_value)
        self._restore_target_value(metric_index)
        if self._targets_changed_callback is not None:
            self._targets_changed_callback()

    def _restore_target_value(self, metric_index: int) -> None:
        if self._unit is None:
            return
        values = [
            self._unit.targetStructure.tph,
            self._unit.targetStructure.ba,
            self._unit.targetStructure.mcs,
            self._unit.targetStructure.cc,
        ]
        metric_kinds = ["tph", "ba", "mcs", "cc"]
        self._updating = True
        self.target_edits[metric_index].setText(
            format_value(metric_kinds[metric_index], values[metric_index], self._unit_system)
        )
        self._updating = False

    def _update_treated_column(self, treated: StructureSummary | None) -> None:
        if treated is None:
            for label in self.treated_values:
                label.setText("-")
            self.treated_hint.setText("Run a treatment to populate treated values.")
            self._set_treated_visual_state(False)
            return

        for label_widget, metric_kind, value in zip(
            self.treated_values,
            ["tph", "ba", "mcs", "cc"],
            [treated.tph, treated.ba, treated.mcs, treated.cc],
        ):
            label_widget.setText(format_value(metric_kind, value, self._unit_system))
        self.treated_hint.setText("Treated values reflect the most recent stand treatment.")
        self._set_treated_visual_state(True)

    def _set_treated_visual_state(self, active: bool) -> None:
        if active:
            header_style = (
                "color: #1e4f2b;"
                "background-color: #e6f4ea;"
                "font-weight: 700;"
                "padding: 2px 6px;"
                "border-radius: 4px;"
            )
            value_style = (
                "color: #1e4f2b;"
                "background-color: #e6f4ea;"
                "padding: 2px 6px;"
                "border-radius: 4px;"
            )
            hint_style = "color: #1e4f2b; font-style: italic;"
        else:
            header_style = (
                "color: #7a7a7a;"
                "background-color: #f3f3f3;"
                "font-weight: 600;"
                "padding: 2px 6px;"
                "border-radius: 4px;"
            )
            value_style = (
                "color: #7a7a7a;"
                "background-color: #f3f3f3;"
                "padding: 2px 6px;"
                "border-radius: 4px;"
            )
            hint_style = "color: #7a7a7a; font-style: italic;"

        self.treated_header.setStyleSheet(header_style)
        self.treated_hint.setStyleSheet(hint_style)
        for label in self.treated_values:
            label.setStyleSheet(value_style)


class UnitSidebar(QWidget):
    def __init__(self, rx_units: list[RxUnit], unit_system: UnitSystem) -> None:
        super().__init__()
        self.model = UnitListModel(rx_units, unit_system)
        self.unit_list_view = QListView()
        self.unit_list_view.setModel(self.model)
        self.structure_info = StructureInfo(unit_system)

        units_group = QGroupBox("UNITS")
        units_layout = QVBoxLayout()
        units_layout.addWidget(self.unit_list_view)
        units_group.setLayout(units_layout)

        structure_group = QGroupBox("Structure")
        structure_layout = QVBoxLayout()
        structure_layout.addWidget(self.structure_info)
        structure_group.setLayout(structure_layout)

        layout = QVBoxLayout()
        layout.addWidget(units_group)
        layout.addWidget(structure_group)
        layout.addStretch(1)
        self.setMaximumWidth(290)
        self.setLayout(layout)
