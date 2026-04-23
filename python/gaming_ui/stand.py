from __future__ import annotations

from collections.abc import Callable

import numpy as np
from matplotlib.figure import Figure
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from PySide6.QtWidgets import QHBoxLayout, QWidget

from rxgaming_core import RxUnit, StructureSummary, TreatmentEngine, TreatmentResult
from .sidebar import UnitSidebar
from .state import StandViewState
from .units import array_to_display, dbh_from_display, dbh_to_display, display_name_for, format_value
from .views import CutReportTab, StandPages, TreatmentReportTab, VisualizeTab


class StandViewCoordinator(QWidget):
    def __init__(self, rx_units: list[RxUnit], state: StandViewState) -> None:
        super().__init__()
        self.rx_units = rx_units
        self.state = state
        self.unit_system = state.unit_system
        self._state_changed_callback: Callable[[], None] | None = None
        self.treatment_engine = TreatmentEngine()

        self.sidebar = UnitSidebar(rx_units, self.unit_system)
        self.visualize_tab = VisualizeTab(self.unit_system)
        self.treatment_report_tab = TreatmentReportTab(self.unit_system)
        self.cut_report_tab = CutReportTab(self.unit_system)
        self.pages = StandPages(self.visualize_tab, self.treatment_report_tab, self.cut_report_tab)

        layout = QHBoxLayout()
        layout.addWidget(self.sidebar)
        layout.addWidget(self.pages, stretch=1)
        self.setLayout(layout)

        self._connect_signals()
        self._restore_state()

    def set_state_changed_callback(self, callback: Callable[[], None]) -> None:
        self._state_changed_callback = callback

    @property
    def current_figure(self) -> Figure:
        if self.state.active_page == 1:
            return self.treatment_report_tab.displayed_ba_figure
        if self.state.active_page == 2:
            return self.cut_report_tab.cut_figure
        return self.visualize_tab.raster_figure

    def current_export_widget(self) -> QWidget:
        if self.state.active_page == 1:
            return self.treatment_report_tab
        if self.state.active_page == 2:
            return self.cut_report_tab
        return self.visualize_tab

    def current_unit_index(self) -> int:
        return self.state.selected_unit_index

    def current_unit(self) -> RxUnit:
        if not self.rx_units:
            raise IndexError("No units are available")
        index = max(0, min(self.state.selected_unit_index, len(self.rx_units) - 1))
        return self.rx_units[index]

    def current_raster_mode(self) -> int:
        return self.state.raster_mode

    def showing_treatment_view(self) -> bool:
        return self.state.show_treatment

    def dbh_min(self) -> float:
        return self.state.dbh_min

    def dbh_max(self) -> float:
        return self.state.dbh_max

    def preview_dbh(self) -> float:
        return self.state.preview_dbh

    def current_page(self) -> int:
        return self.state.active_page

    def current_points(self) -> np.ndarray:
        unit = self.current_unit()
        if self.state.show_treatment:
            return unit.get_treat_taos()
        return unit.get_taos()

    def current_raster_array(self) -> np.ndarray:
        unit = self.current_unit()
        if self.state.raster_mode == 0:
            return unit.get_treat_chm() if self.state.show_treatment else unit.get_chm()
        if self.state.raster_mode == 1:
            return unit.get_treat_basin() if self.state.show_treatment else unit.get_basin()
        if self.state.raster_mode == 2:
            return unit.get_treat_clump_map() if self.state.show_treatment else unit.get_clump_map()
        if self.state.raster_mode == 3:
            return unit.get_treat_hillshade() if self.state.show_treatment else unit.get_hillshade()

        points = self.current_points()
        if points.ndim != 2 or points.shape[1] < 2:
            return np.zeros((0, 2))
        return points[:, :2]

    def select_unit(self, index: int) -> None:
        if not self.rx_units:
            return
        row = max(0, min(index, len(self.rx_units) - 1))
        self.sidebar.unit_list_view.setCurrentIndex(self.sidebar.model.index(row, 0))

    def _connect_signals(self) -> None:
        selection_model = self.sidebar.unit_list_view.selectionModel()
        selection_model.currentChanged.connect(self._on_unit_changed)

        self.sidebar.structure_info.set_targets_changed_callback(self._on_targets_changed)
        self.visualize_tab.raster_mode.currentIndexChanged.connect(self._on_raster_mode_changed)
        self.visualize_tab.dbh_cutoff.valueChanged.connect(self._on_cutoff_changed)
        self.visualize_tab.show_treatment_button.toggled.connect(self._on_show_treatment_toggled)
        self.treatment_report_tab.preview_slider.valueChanged.connect(self._on_preview_changed)
        self.pages.currentChanged.connect(self._on_page_changed)

    def _restore_state(self) -> None:
        self.visualize_tab.raster_mode.setCurrentIndex(self.state.raster_mode)
        self.visualize_tab.dbh_cutoff.setValue(int(round(dbh_to_display(self.state.dbh_cutoff, self.unit_system))))
        self.visualize_tab.show_treatment_button.setChecked(self.state.show_treatment)
        self.treatment_report_tab.preview_slider.setValue(
            int(round(dbh_to_display(self.state.preview_dbh, self.unit_system)))
        )
        self.pages.setCurrentIndex(self.state.active_page)

        if self.rx_units:
            self.select_unit(self.state.selected_unit_index)
        self.refresh_all()

    def _on_unit_changed(self, current: object, previous: object) -> None:
        del previous
        row = current.row() if getattr(current, "isValid", lambda: False)() else 0
        self.state.selected_unit_index = max(0, row)
        self.refresh_all()
        self._notify_state_changed()

    def _on_targets_changed(self) -> None:
        self.sidebar.model.refresh()
        self.refresh_all()
        self._notify_state_changed()

    def _on_raster_mode_changed(self, index: int) -> None:
        self.state.raster_mode = index
        self.update_raster_canvas()
        self._notify_state_changed()

    def _on_cutoff_changed(self, value: int) -> None:
        self.state.dbh_cutoff = dbh_from_display(float(value), self.unit_system)
        self.update_treatment_report()
        self._notify_state_changed()

    def _on_show_treatment_toggled(self, checked: bool) -> None:
        if checked and self.rx_units:
            self.treatment_engine.do_treatment(self.current_unit(), self.dbh_min(), self.dbh_max())
            self.sidebar.model.refresh()
        self.state.show_treatment = checked
        self.refresh_all()
        self._notify_state_changed()

    def _on_preview_changed(self, value: int) -> None:
        self.state.preview_dbh = dbh_from_display(float(value), self.unit_system)
        self.update_treatment_report()
        self._notify_state_changed()

    def _on_page_changed(self, index: int) -> None:
        self.state.active_page = index
        self.refresh_current_page()
        self._notify_state_changed()

    def refresh_all(self) -> None:
        if not self.rx_units:
            self.update_raster_canvas()
            self.update_treatment_report()
            self.update_cut_report()
            return
        self.sidebar.structure_info.update_for_unit(self.current_unit())
        self.update_raster_canvas()
        self.update_treatment_report()
        self.update_cut_report()

    def refresh_current_page(self) -> None:
        if self.state.active_page == 1:
            self.update_treatment_report()
        elif self.state.active_page == 2:
            self.update_cut_report()
        else:
            self.update_raster_canvas()

    def update_raster_canvas(self) -> None:
        figure = self.visualize_tab.raster_figure
        figure.clear()
        axes = figure.add_subplot(111)
        self.visualize_tab.raster_axes = axes

        if not self.rx_units:
            axes.text(0.5, 0.5, "No units available", ha="center", va="center")
        else:
            unit = self.current_unit()
            data = self.current_raster_array()
            hillshade = unit.get_treat_hillshade() if self.state.show_treatment else unit.get_hillshade()
            colorbar = None

            if data.size == 0:
                axes.text(0.5, 0.5, "No raster data", ha="center", va="center")

            if self.state.raster_mode == 0:
                display_data = array_to_display("canopy_height", data, self.unit_system)
                img = axes.imshow(display_data, cmap="coolwarm", vmin=0)
                axes.imshow(hillshade, cmap="Greys", alpha=0.5)
                colorbar = figure.colorbar(
                    img,
                    orientation="vertical",
                    label=display_name_for("canopy_height", self.unit_system),
                )
            elif self.state.raster_mode == 1:
                img = axes.imshow(data, cmap="nipy_spectral")
                axes.imshow(hillshade, cmap="Greys", alpha=0.5)
                colorbar = figure.colorbar(img, orientation="vertical", label="Basin ID (unique values)")
            elif self.state.raster_mode == 2:
                colors = ("white", "#7bc043", "#fdf498", "#f37736", "#ee4035")
                cmap = LinearSegmentedColormap.from_list("Clump Colors", colors, 5)
                boundaries = (-0.5, 0.5, 1.5, 4.5, 9.5, 99)
                norm = BoundaryNorm(boundaries, len(boundaries))
                img = axes.imshow(data, cmap=cmap, norm=norm)
                if hillshade.size:
                    axes.imshow(hillshade, cmap="Greys", alpha=0.5)
                colorbar = figure.colorbar(
                    img,
                    orientation="vertical",
                    label="Clump Map (Clump bins)",
                    ticks=[0, 1, 3, 7, 55],
                )
                colorbar.ax.set_yticklabels(["0", "1", "2-4", "4-9", "10+"])
            else:
                if data.size == 0:
                    axes.text(0.5, 0.5, "No raster data", ha="center", va="center")
                else:
                    img = axes.imshow(data, cmap="Greys")
                    colorbar = figure.colorbar(img, orientation="vertical", label="Hillshade")

            suffix = " (Treated)" if self.state.show_treatment else ""
            axes.set_title(f"{unit.name} {self.visualize_tab.raster_mode.currentText()}{suffix}")
            if colorbar is not None:
                colorbar.ax.tick_params(labelsize=8)

        axes.set_xticks([])
        axes.set_yticks([])
        figure.tight_layout()
        self.visualize_tab.raster_canvas.draw_idle()

    def update_treatment_report(self) -> None:
        report = self.treatment_report_tab
        report.current_ba_axes.clear()
        report.current_mcs_axes.clear()
        report.displayed_ba_axes.clear()
        report.displayed_mcs_axes.clear()

        if not self.rx_units:
            report.report_status.setText("No units available.")
            report.target_summary.setText("-")
            for canvas in (
                report.current_ba_canvas,
                report.current_mcs_canvas,
                report.displayed_ba_canvas,
                report.displayed_mcs_canvas,
            ):
                canvas.draw_idle()
            return

        unit = self.current_unit()
        #simulated = unit.get_simulated_structures(self.state.preview_dbh)
        preview = None
        
        report.report_status.setText(
            f"{display_name_for('dbh', self.unit_system)}: {format_value('dbh', self.state.preview_dbh, self.unit_system, precision=0)}\n"
            f"Unit: {unit.name}\n"
            f"Treatment Result: {unit.result.name}"
        )
        report.target_summary.setText(
            self._format_target_summary(
                unit.currentStructure,
                unit.targetStructure,
                unit.treatedStructure,
                preview,
                self.unit_system,
            )
        )

        self._draw_comparison_bar(
            report.current_ba_axes,
            f"Current vs Target {display_name_for('ba', self.unit_system)}",
            [
                array_to_display("ba", np.asarray([unit.currentStructure.ba], dtype=float), self.unit_system)[0],
                array_to_display("ba", np.asarray([unit.targetStructure.ba], dtype=float), self.unit_system)[0],
            ],
            ["Current", "Target"],
            ylabel=display_name_for("ba", self.unit_system),
        )
        self._draw_comparison_bar(
            report.current_mcs_axes,
            "Current vs Target MCS",
            [unit.currentStructure.mcs, unit.targetStructure.mcs],
            ["Current", "Target"],
        )
        self._draw_comparison_bar(
            report.displayed_ba_axes,
            f"Preview vs Treated {display_name_for('ba', self.unit_system)}",
            [
                array_to_display(
                    "ba",
                    np.asarray([preview.ba if preview is not None else unit.currentStructure.ba], dtype=float),
                    self.unit_system,
                )[0],
                array_to_display("ba", np.asarray([unit.treatedStructure.ba], dtype=float), self.unit_system)[0],
            ],
            ["Preview", "Treated"],
            ylabel=display_name_for("ba", self.unit_system),
        )
        self._draw_comparison_bar(
            report.displayed_mcs_axes,
            "Preview vs Treated MCS",
            [preview.mcs if preview is not None else unit.currentStructure.mcs, unit.treatedStructure.mcs],
            ["Preview", "Treated"],
        )

        report.current_ba_figure.tight_layout()
        report.current_mcs_figure.tight_layout()
        report.displayed_ba_figure.tight_layout()
        report.displayed_mcs_figure.tight_layout()
        report.current_ba_canvas.draw_idle()
        report.current_mcs_canvas.draw_idle()
        report.displayed_ba_canvas.draw_idle()
        report.displayed_mcs_canvas.draw_idle()

    def update_cut_report(self) -> None:
        report = self.cut_report_tab
        report.cut_axes.clear()

        if not self.rx_units:
            report.cut_summary.setText("Cut Trees:\nNo units available.")
            report.cut_canvas.draw_idle()
            return

        unit = self.current_unit()
        cut_points = unit.get_cut_taos()
        if cut_points.ndim != 2 or cut_points.shape[0] == 0:
            report.cut_summary.setText(f"Cut Trees:\n{unit.name}\nNo cut trees are available yet.")
            report.cut_axes.text(0.5, 0.5, "No cut-tree data yet", ha="center", va="center")
            report.cut_axes.set_xticks([])
            report.cut_axes.set_yticks([])
        else:
            heights = cut_points[:, 2] if cut_points.shape[1] > 2 else cut_points[:, -1]
            display_heights = array_to_display("height", heights.astype(float), self.unit_system)
            report.cut_summary.setText(
                f"Cut Trees:\n{unit.name}\nCount: {cut_points.shape[0]}\n"
                f"Mean Height: {format_value('height', float(np.mean(heights)), self.unit_system)} "
                f"{display_name_for('height', self.unit_system).split()[-1][1:-1]}"
            )
            report.cut_axes.hist(display_heights, bins=20, color="#4a6fa5", edgecolor="white")
            report.cut_axes.set_title("Cut Tree Height Distribution")
            report.cut_axes.set_xlabel(display_name_for("height", self.unit_system))
            report.cut_axes.set_ylabel("Count")

        report.cut_figure.tight_layout()
        report.cut_canvas.draw_idle()

    def _notify_state_changed(self) -> None:
        if self._state_changed_callback is not None:
            self._state_changed_callback()

    @staticmethod
    def _draw_comparison_bar(
        axes: object,
        title: str,
        values: list[float],
        labels: list[str],
        ylabel: str | None = None,
    ) -> None:
        axes.bar(labels, values, color=["#6d8dad", "#c88c5a"])
        axes.set_title(title)
        if ylabel is not None:
            axes.set_ylabel(ylabel)
        axes.tick_params(axis="x", rotation=15)

    @staticmethod
    def _format_target_summary(
        current: StructureSummary,
        target: StructureSummary,
        treated: StructureSummary | None,
        preview: StructureSummary | None,
        unit_system: object,
    ) -> str:
        def structure_text(structure: StructureSummary | None) -> str:
            if structure is None:
                return "-"
            return ", ".join(
                [
                    f"{display_name_for('tph', unit_system)} {format_value('tph', structure.tph, unit_system)}",
                    f"{display_name_for('ba', unit_system)} {format_value('ba', structure.ba, unit_system)}",
                    f"MCS {format_value('mcs', structure.mcs, unit_system)}",
                    f"CC {format_value('cc', structure.cc, unit_system)}",
                ]
            )

        return (
            "Current\n"
            f"{structure_text(current)}\n\n"
            "Target\n"
            f"{structure_text(target)}\n\n"
            "Preview\n"
            f"{structure_text(preview)}\n\n"
            "Treated\n"
            f"{structure_text(treated)}"
        )
