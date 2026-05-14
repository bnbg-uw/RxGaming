from __future__ import annotations

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide6.QtCore import QTimer, Qt
from PySide6.QtGui import QResizeEvent, QShowEvent
from PySide6.QtWidgets import (
    QGridLayout,
    QLabel,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)
from .units import UnitSystem


class VisualizeTab(QWidget):
    _RASTER_COLORBAR_PAD = 0.02
    _RASTER_COLORBAR_WIDTH = 0.03

    def __init__(self, unit_system: UnitSystem) -> None:
        super().__init__()
        self._unit_system = unit_system

        self.raster_figure = Figure()
        self.raster_canvas = FigureCanvas(self.raster_figure)
        self.raster_axes = self.raster_figure.add_subplot(111)
        # Keep both the raster and colorbar in fixed axes so redraws do not
        # resize the image area as colorbars are rebuilt or updated.
        self.raster_figure.subplots_adjust(left=0.03, right=0.84, top=0.94, bottom=0.04)
        self.raster_colorbar_axes = self.raster_figure.add_axes([0.86, 0.04, 0.03, 0.90])
        self.raster_colorbar_axes.set_visible(False)
        self.sync_raster_colorbar_axes()

        layout = QVBoxLayout()
        layout.addWidget(self.raster_canvas, stretch=1)
        self.setLayout(layout)
        QTimer.singleShot(0, self._sync_raster_colorbar_after_layout)

    def sync_raster_colorbar_axes(self) -> None:
        raster_position = self.raster_axes.get_position()
        left = min(
            raster_position.x1 + self._RASTER_COLORBAR_PAD,
            1.0 - self._RASTER_COLORBAR_WIDTH,
        )
        self.raster_colorbar_axes.set_position(
            [
                left,
                raster_position.y0,
                self._RASTER_COLORBAR_WIDTH,
                raster_position.height,
            ]
        )

    def _sync_raster_colorbar_after_layout(self) -> None:
        self.sync_raster_colorbar_axes()
        self.raster_canvas.draw_idle()

    def showEvent(self, event: QShowEvent) -> None:
        super().showEvent(event)
        QTimer.singleShot(0, self._sync_raster_colorbar_after_layout)

    def resizeEvent(self, event: QResizeEvent) -> None:
        super().resizeEvent(event)
        self.sync_raster_colorbar_axes()


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
        self.stand_area_label = QLabel("")
        self.stand_area_label.setStyleSheet("font: 16pt;")
        self.stand_area_label.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.stand_area_label.setWordWrap(True)

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
        layout_grid.addWidget(self.stand_area_label, 4, 0, 1, 3)
        layout_grid.setRowStretch(0, 0)
        layout_grid.setRowStretch(1, 2)
        layout_grid.setRowStretch(2, 2)
        layout_grid.setRowStretch(3, 1)
        layout_grid.setRowStretch(4, 0)

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
