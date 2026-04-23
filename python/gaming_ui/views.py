from __future__ import annotations

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from widgets import SliderWithValue
from .units import UnitSystem, dbh_to_display, display_name_for


class VisualizeTab(QWidget):
    def __init__(self, unit_system: UnitSystem) -> None:
        super().__init__()
        self._unit_system = unit_system

        self.raster_figure = Figure()
        self.raster_canvas = FigureCanvas(self.raster_figure)
        self.raster_axes = self.raster_figure.add_subplot(111)

        self.dbh_cutoff = SliderWithValue(Qt.Orientation.Horizontal)
        self.dbh_cutoff.setMinimum(0)
        self.dbh_cutoff.setMaximum(120)
        self.dbh_cutoff.setValue(int(round(dbh_to_display(76.2, unit_system))))

        self.raster_mode = QComboBox()
        self.raster_mode.addItems(["Canopy Model", "Basins", "Clumps"])

        self.show_treatment_button = QPushButton("Show Treatment")
        self.show_treatment_button.setCheckable(True)
        self.show_treatment_button.setStyleSheet(
            "QPushButton:checked { background-color: rgb(80, 80, 80); color: white; border: none; }"
        )

        controls = QGroupBox("Stand Controls")
        controls_layout = QGridLayout()
        controls_layout.addWidget(QLabel(display_name_for("dbh", unit_system)), 0, 0)
        controls_layout.addWidget(self.dbh_cutoff, 0, 1, 1, 3)
        controls_layout.addWidget(QLabel("View Mode:"), 1, 0)
        controls_layout.addWidget(self.raster_mode, 1, 1, 1, 3)
        controls_layout.addWidget(self.show_treatment_button, 2, 1, 1, 2)
        controls.setLayout(controls_layout)

        layout = QVBoxLayout()
        layout.addWidget(self.raster_canvas, stretch=1)
        layout.addWidget(controls)
        self.setLayout(layout)


class TreatmentReportTab(QWidget):
    def __init__(self, unit_system: UnitSystem) -> None:
        super().__init__()
        self._unit_system = unit_system

        self.report_label = QLabel("Treatment Report")
        self.report_label.setStyleSheet("font: 24pt;")
        self.current_label = QLabel("Current")
        self.current_label.setStyleSheet("font: 16pt;")
        self.displayed_label = QLabel("Post-Treatment")
        self.displayed_label.setStyleSheet("font: 16pt;")
        self.target_label = QLabel("Target")
        self.target_label.setStyleSheet("font: 16pt;")

        self.preview_slider = SliderWithValue(Qt.Orientation.Horizontal)
        self.preview_slider.setMinimum(0)
        self.preview_slider.setMaximum(120)
        self.preview_slider.setValue(int(round(dbh_to_display(76.2, unit_system))))

        self.current_ba_figure = Figure()
        self.current_ba_canvas = FigureCanvas(self.current_ba_figure)
        self.current_ba_axes = self.current_ba_figure.add_subplot(111)

        self.current_mcs_figure = Figure()
        self.current_mcs_canvas = FigureCanvas(self.current_mcs_figure)
        self.current_mcs_axes = self.current_mcs_figure.add_subplot(111)

        self.displayed_ba_figure = Figure()
        self.displayed_ba_canvas = FigureCanvas(self.displayed_ba_figure)
        self.displayed_ba_axes = self.displayed_ba_figure.add_subplot(111)

        self.displayed_mcs_figure = Figure()
        self.displayed_mcs_canvas = FigureCanvas(self.displayed_mcs_figure)
        self.displayed_mcs_axes = self.displayed_mcs_figure.add_subplot(111)

        self.target_summary = QLabel("-")
        self.target_summary.setStyleSheet("font: 14pt;")
        self.target_summary.setWordWrap(True)
        self.report_status = QLabel("")
        self.report_status.setWordWrap(True)

        header_layout = QHBoxLayout()
        header_layout.addWidget(QLabel(display_name_for("dbh", unit_system)))
        header_layout.addWidget(self.preview_slider)

        grid = QGridLayout()
        grid.addWidget(self.report_label, 0, 1)
        grid.addWidget(self.current_label, 1, 0)
        grid.addWidget(self.displayed_label, 2, 0)
        grid.addWidget(self.target_label, 3, 0)
        grid.addWidget(self.current_ba_canvas, 1, 1)
        grid.addWidget(self.current_mcs_canvas, 1, 2)
        grid.addWidget(self.displayed_ba_canvas, 2, 1)
        grid.addWidget(self.displayed_mcs_canvas, 2, 2)
        grid.addWidget(self.target_summary, 3, 1, 1, 2)
        grid.setRowStretch(1, 2)
        grid.setRowStretch(2, 2)

        layout = QVBoxLayout()
        layout.addLayout(header_layout)
        layout.addWidget(self.report_status)
        layout.addLayout(grid)
        self.setLayout(layout)


class CutReportTab(QWidget):
    def __init__(self, unit_system: UnitSystem) -> None:
        super().__init__()
        self._unit_system = unit_system

        self.page_label = QLabel("Cut Trees Info")
        self.page_label.setStyleSheet("font: 16pt;")
        self.cut_summary = QLabel("Cut Trees:")
        self.cut_summary.setStyleSheet("font: 16pt;")
        self.cut_summary.setWordWrap(True)

        self.cut_figure = Figure()
        self.cut_canvas = FigureCanvas(self.cut_figure)
        self.cut_axes = self.cut_figure.add_subplot(111)

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
