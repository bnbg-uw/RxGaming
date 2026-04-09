from __future__ import annotations

from typing import Any

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QListView,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from rxgaming_core import ProjectSettings, ProjectArea, RxUnit, StructureSummary, TreatmentEngine, TreatmentResult
from widgets import SliderWithValue


class UnitTableModel(QAbstractTableModel):
    headers = ["Unit", "Area (ha)", "BA", "TPH", "MCS", "CC", "Result"]

    def __init__(self, rx_units: list[RxUnit]) -> None:
        super().__init__()
        self.rx_units = rx_units

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self.rx_units)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self.headers)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid() or role != Qt.ItemDataRole.DisplayRole:
            return None
        unit = self.rx_units[index.row()]
        current = unit.currentStructure
        if index.column() == 0:
            return unit.name
        if index.column() == 1:
            return f"{unit.areaHa:.2f}"
        if index.column() == 2:
            return f"{current.ba:.2f}"
        if index.column() == 3:
            return f"{current.tph:.2f}"
        if index.column() == 4:
            return f"{current.mcs:.2f}"
        if index.column() == 5:
            return f"{current.cc:.2f}"
        if index.column() == 6:
            return unit.result.name
        return None

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.DisplayRole:
            return self.headers[section]
        return super().headerData(section, orientation, role)

    def refresh(self) -> None:
        self.layoutChanged.emit()


class StructureInfo(QWidget):
    def __init__(self) -> None:
        super().__init__()
        layout = QVBoxLayout()
        self.current_label = QLabel()
        self.target_label = QLabel()
        self.treated_label = QLabel()
        self.result_label = QLabel()
        layout.addWidget(self.current_label)
        layout.addWidget(self.target_label)
        layout.addWidget(self.treated_label)
        layout.addWidget(self.result_label)
        layout.addStretch(1)
        self.setLayout(layout)

    def update_for_unit(self, unit: RxUnit) -> None:
        self.current_label.setText(self._format_structure("Current", unit.currentStructure))
        self.target_label.setText(self._format_structure("Target", unit.targetStructure))
        treated = unit.treatedStructure if unit.result == TreatmentResult.success else None
        self.treated_label.setText(self._format_structure("Treated", treated))
        self.result_label.setText(f"Result: {unit.result.name}")

    @staticmethod
    def _format_structure(title: str, structure: StructureSummary | None) -> str:
        if structure is None:
            return f"{title}\n-"
        return (
            f"{title}\n"
            f"BA: {structure.ba:.2f}\n"
            f"TPH: {structure.tph:.2f}\n"
            f"MCS: {structure.mcs:.2f}\n"
            f"CC: {structure.cc:.2f}"
        )


class GamingTabs(QTabWidget):
    def __init__(
        self,
        project_settings: ProjectSettings,
        project_area: ProjectArea,
        saved_state: dict[str, Any],
    ) -> None:
        super().__init__(None)
        self.project_settings = project_settings
        self.project_area = project_area
        self.rx_units = list(project_area.rxUnits)
        self.treatment_engine = TreatmentEngine()

        self.stand_tab = QWidget()
        self.land_tab = QWidget()
        self.addTab(self.stand_tab, "Stand View")
        self.addTab(self.land_tab, "Landscape View")

        self._build_stand_tab(saved_state)
        self._build_landscape_tab()
        self._restore_state(saved_state)
        self._update_all_views()

    def _build_stand_tab(self, saved_state: dict[str, Any]) -> None:
        outer_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        center_layout = QVBoxLayout()
        right_layout = QVBoxLayout()

        self.unit_label = QLabel("UNITS")
        self.unit_list_view = QListView()
        self.model = UnitTableModel(self.rx_units)
        self.unit_list_view.setModel(self.model)
        self.structure_info = StructureInfo()

        left_layout.addWidget(self.unit_label)
        left_layout.addWidget(self.unit_list_view)
        left_layout.addWidget(self.structure_info)

        self.raster_mode = QComboBox()
        self.raster_mode.addItems(["CHM", "Basin", "Clump Map", "Hillshade", "TAOs"])

        self.show_treatment_button = QPushButton("Show Treatment")
        self.show_treatment_button.setCheckable(True)

        self.dbh_min_spin = QDoubleSpinBox()
        self.dbh_min_spin.setRange(0.0, 1000.0)
        self.dbh_min_spin.setDecimals(2)
        self.dbh_max_spin = QDoubleSpinBox()
        self.dbh_max_spin.setRange(0.0, 1000.0)
        self.dbh_max_spin.setDecimals(2)
        self.dbh_max_spin.setValue(999.0)

        self.preview_dbh_slider = SliderWithValue(Qt.Orientation.Horizontal)
        self.preview_dbh_slider.setRange(1, 40)
        self.preview_dbh_slider.setValue(10)

        self.run_treatment_button = QPushButton("Run Treatment")

        self.raster_figure = Figure(figsize=(8, 8))
        self.raster_canvas = FigureCanvas(self.raster_figure)
        self.raster_axes = self.raster_figure.add_subplot(111)

        center_layout.addWidget(self.raster_mode)
        center_layout.addWidget(self.show_treatment_button)
        center_layout.addWidget(self.raster_canvas)

        form_layout = QFormLayout()
        form_layout.addRow("DBH Min", self.dbh_min_spin)
        form_layout.addRow("DBH Max", self.dbh_max_spin)
        form_layout.addRow("Preview DBH", self.preview_dbh_slider)

        right_layout.addLayout(form_layout)
        right_layout.addWidget(self.run_treatment_button)
        right_layout.addStretch(1)

        outer_layout.addLayout(left_layout, 1)
        outer_layout.addLayout(center_layout, 3)
        outer_layout.addLayout(right_layout, 1)
        self.stand_tab.setLayout(outer_layout)

        self.unit_list_view.selectionModel().currentChanged.connect(self._on_unit_changed)
        self.raster_mode.currentIndexChanged.connect(lambda _index: self.update_raster_canvas())
        self.show_treatment_button.toggled.connect(lambda _checked: self.update_raster_canvas())
        self.run_treatment_button.clicked.connect(self.run_treatment)

    def _build_landscape_tab(self) -> None:
        layout = QVBoxLayout()
        self.landscape_figure = Figure(figsize=(8, 6))
        self.landscape_canvas = FigureCanvas(self.landscape_figure)
        self.landscape_axes = self.landscape_figure.add_subplot(111)
        layout.addWidget(self.landscape_canvas)
        self.land_tab.setLayout(layout)

    def _restore_state(self, saved_state: dict[str, Any]) -> None:
        row = int(saved_state.get("GamingActivity.selected_unit", 0))
        row = max(0, min(row, len(self.rx_units) - 1)) if self.rx_units else 0
        if self.rx_units:
            self.unit_list_view.setCurrentIndex(self.model.index(row, 0))
        self.raster_mode.setCurrentIndex(int(saved_state.get("GamingActivity.raster_mode", 0)))
        self.show_treatment_button.setChecked(bool(saved_state.get("GamingActivity.show_treatment", False)))
        self.dbh_min_spin.setValue(float(saved_state.get("GamingActivity.dbh_min", 0.0)))
        self.dbh_max_spin.setValue(float(saved_state.get("GamingActivity.dbh_max", 999.0)))

    def _update_all_views(self) -> None:
        self.update_structure_info()
        self.update_raster_canvas()
        self.update_landscape_plot()

    def _on_unit_changed(self, current: QModelIndex, previous: QModelIndex) -> None:
        del previous
        if current.isValid():
            self._update_all_views()

    def current_unit_index(self) -> int:
        index = self.unit_list_view.currentIndex()
        return index.row() if index.isValid() else 0

    def current_unit(self) -> RxUnit:
        return self.rx_units[self.current_unit_index()]

    def current_raster_mode(self) -> int:
        return self.raster_mode.currentIndex()

    def showing_treatment_view(self) -> bool:
        return self.show_treatment_button.isChecked()

    def dbh_min(self) -> float:
        return self.dbh_min_spin.value()

    def dbh_max(self) -> float:
        return self.dbh_max_spin.value()

    def current_points(self) -> np.ndarray:
        unit = self.current_unit()
        if self.showing_treatment_view() and unit.result == TreatmentResult.success:
            return np.asarray(unit.get_treat_taos())
        return np.asarray(unit.get_taos())

    def current_raster_array(self) -> np.ndarray:
        unit = self.current_unit()
        use_treated = self.showing_treatment_view() and unit.result == TreatmentResult.success
        mode = self.current_raster_mode()
        if mode == 0:
            return np.asarray(unit.get_treat_chm() if use_treated else unit.get_chm())
        if mode == 1:
            return np.asarray(unit.get_treat_basin() if use_treated else unit.get_basin())
        if mode == 2:
            return np.asarray(unit.get_treat_clump_map() if use_treated else unit.get_clump_map())
        if mode == 3:
            return np.asarray(unit.get_treat_hillshade() if use_treated else unit.get_hillshade())
        points = self.current_points()
        if points.size == 0:
            return np.zeros((1, 1))
        return points[:, :2]

    def run_treatment(self, checked: bool = False) -> None:
        del checked
        unit = self.current_unit()
        self.treatment_engine.do_treatment(unit, self.dbh_min(), self.dbh_max())
        self.show_treatment_button.setChecked(True)
        self.model.refresh()
        self._update_all_views()

    def update_structure_info(self) -> None:
        self.structure_info.update_for_unit(self.current_unit())

    def update_raster_canvas(self) -> None:
        self.raster_axes.clear()
        data = self.current_raster_array()
        if self.current_raster_mode() == 4:
            if data.ndim == 2 and data.shape[1] >= 2:
                self.raster_axes.scatter(data[:, 0], data[:, 1], s=4)
                self.raster_axes.set_title("TAO Points")
            else:
                self.raster_axes.text(0.5, 0.5, "No point data", ha="center", va="center")
        else:
            self.raster_axes.imshow(data, cmap="viridis")
            self.raster_axes.set_title(self.raster_mode.currentText())
        self.raster_axes.set_xticks([])
        self.raster_axes.set_yticks([])
        self.raster_figure.tight_layout()
        self.raster_canvas.draw_idle()

    def update_landscape_plot(self) -> None:
        self.landscape_axes.clear()
        x = [unit.currentStructure.tph for unit in self.rx_units]
        y = [unit.currentStructure.ba for unit in self.rx_units]
        self.landscape_axes.scatter(x, y, c="#2f6fed")
        self.landscape_axes.set_xlabel("TPH")
        self.landscape_axes.set_ylabel("BA")
        self.landscape_axes.set_title("Landscape Overview")
        self.landscape_figure.tight_layout()
        self.landscape_canvas.draw_idle()
