from __future__ import annotations

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QGridLayout,
    QLabel,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)
from .units import UnitSystem


class VisualizeTab(QWidget):
    def __init__(self, unit_system: UnitSystem) -> None:
        super().__init__()
        self._unit_system = unit_system

        self.raster_figure = Figure()
        self.raster_canvas = FigureCanvas(self.raster_figure)
        self.raster_axes = self.raster_figure.add_subplot(111)

        layout = QVBoxLayout()
        layout.addWidget(self.raster_canvas, stretch=1)
        self.setLayout(layout)


class TreatmentReportTab(QWidget):
    def __init__(self, unit_system: UnitSystem) -> None:
        super().__init__()
        self._unit_system = unit_system

        self.report_label = QLabel("Treatment Report")
        self.report_label.setStyleSheet("font: 24pt;")
        self.current_label = QLabel("Current")
        self.current_label.setStyleSheet("font: 16pt;")
        self.current_label.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.current_label.setWordWrap(True)
        self.displayed_label = QLabel("Post-Treatment")
        self.displayed_label.setStyleSheet("font: 16pt;")
        self.displayed_label.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.displayed_label.setWordWrap(True)
        self.target_label = QLabel("Target")
        self.target_label.setStyleSheet("font: 16pt;")
        self.target_label.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.target_label.setWordWrap(True)

        self.current_ba_figure = Figure()
        self.current_ba_canvas = FigureCanvas(self.current_ba_figure)
        self.current_ba_axes = self.current_ba_figure.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.current_ba_figure.subplots_adjust(left=0.1, right=0.97, top=0.85, bottom=0.25)

        self.current_mcs_figure = Figure()
        self.current_mcs_canvas = FigureCanvas(self.current_mcs_figure)
        self.current_mcs_axes = self.current_mcs_figure.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.current_mcs_figure.subplots_adjust(left=0.1, right=0.97, top=0.85, bottom=0.25)

        self.displayed_ba_figure = Figure()
        self.displayed_ba_canvas = FigureCanvas(self.displayed_ba_figure)
        self.displayed_ba_axes = self.displayed_ba_figure.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.displayed_ba_figure.subplots_adjust(left=0.1, right=0.97, top=0.85, bottom=0.25)

        self.displayed_mcs_figure = Figure()
        self.displayed_mcs_canvas = FigureCanvas(self.displayed_mcs_figure)
        self.displayed_mcs_axes = self.displayed_mcs_figure.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.displayed_mcs_figure.subplots_adjust(left=0.1, right=0.97, top=0.85, bottom=0.25)

        self.displayed_mcs_prop = QLabel("")
        self.displayed_mcs_prop.setStyleSheet("font: 16pt;")
        self.displayed_mcs_prop.setWordWrap(True)
        self.displayed_mcs_prop.setMinimumHeight(10)
        self.displayed_mcs_canvas.setMinimumHeight(10)

        layout_grid = QGridLayout()
        layout_grid.addWidget(self.report_label, 0, 1)
        layout_grid.addWidget(self.current_label, 1, 0)
        layout_grid.addWidget(self.displayed_label, 2, 0)
        layout_grid.addWidget(self.target_label, 3, 0)
        layout_grid.addWidget(self.current_ba_canvas, 1, 1)
        layout_grid.addWidget(self.current_mcs_canvas, 1, 2)
        layout_grid.addWidget(self.displayed_ba_canvas, 2, 1)
        layout_grid.addWidget(self.displayed_mcs_canvas, 2, 2)
        layout_grid.addWidget(self.displayed_mcs_prop, 3, 2)
        layout_grid.setRowStretch(0, 0)
        layout_grid.setRowStretch(1, 2)
        layout_grid.setRowStretch(2, 2)
        layout_grid.setRowStretch(3, 1)

        layout = QVBoxLayout()
        layout.addLayout(layout_grid)
        self.setLayout(layout)


class CutReportTab(QWidget):
    def __init__(self, unit_system: UnitSystem) -> None:
        super().__init__()
        self._unit_system = unit_system

        self.page_label = QLabel("Cut trees info")
        self.page_label.setStyleSheet("font: 16pt;")
        self.cut_summary = QLabel("Cut Trees:")
        self.cut_summary.setStyleSheet("font: 16pt;")
        self.cut_summary.setWordWrap(True)

        self.cut_figure = Figure()
        self.cut_canvas = FigureCanvas(self.cut_figure)
        self.cut_axes = self.cut_figure.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.cut_figure.subplots_adjust(left=0.1, right=0.97, top=0.85, bottom=0.25)

        layout = QGridLayout()
        layout.addWidget(self.page_label, 0, 1)
        layout.addWidget(self.cut_summary, 1, 0)
        layout.addWidget(self.cut_canvas, 1, 1)
        layout.setColumnStretch(1, 1)
        layout.setRowStretch(1, 1)
        self.setLayout(layout)


class StandPages(QTabWidget):
    def __init__(self, visualize: VisualizeTab, treatment_report: TreatmentReportTab, cut_report: CutReportTab) -> None:
        super().__init__()
        self.addTab(visualize, "Visualize")
        self.addTab(treatment_report, "Treatment Report")
        self.addTab(cut_report, "Cut Report")
