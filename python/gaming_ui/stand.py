from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap, hsv_to_rgb
from scipy.stats import gaussian_kde
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QHBoxLayout, QWidget

from rxgaming_core import RxUnit, StructureSummary, TreatmentEngine
from .sidebar import UnitSidebar
from .state import StandViewState
from .units import FEET_PER_METER, array_to_display, dbh_from_display, dbh_to_display, display_name_for, format_value, label_for
from .views import CutReportTab, StandPages, TreatmentReportTab, VisualizeTab


class StandViewCoordinator(QWidget):
    _SLIDER_DEBOUNCE_MS = 150
    _BASIN_PALETTE_SIZE = 64

    def __init__(self, rx_units: list[RxUnit], state: StandViewState) -> None:
        super().__init__()
        self.rx_units = rx_units
        self.state = state
        self.unit_system = state.unit_system
        self._state_changed_callback: Callable[[str], None] | None = None
        self._landscape_invalidated_callback: Callable[[], None] | None = None
        self.treatment_engine = TreatmentEngine()
        self._mcs_prop_table = self._load_mcs_prop_table()
        self._treated_unit_ids: set[int] = set()
        self._pending_dbh_cutoff = self.state.dbh_cutoff
        self._dbh_update_timer = QTimer(self)
        self._dbh_update_timer.setSingleShot(True)
        self._dbh_update_timer.setInterval(self._SLIDER_DEBOUNCE_MS)
        self._dbh_update_timer.timeout.connect(self._flush_cutoff_change)
        self._raster_dirty = True
        self._treatment_report_dirty = True
        self._cut_report_dirty = True
        self._raster_colorbar = None
        self._raster_image = None
        self._hillshade_image = None
        self._raster_text = None
        self._raster_mode_key: tuple[int, bool] | None = None
        self._raster_scale_key: tuple[float | int | None, float | int | None] | None = None
        self._basin_colormap = self._build_basin_colormap()

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

    def set_state_changed_callback(self, callback: Callable[[str], None]) -> None:
        self._state_changed_callback = callback

    def set_landscape_invalidated_callback(self, callback: Callable[[], None]) -> None:
        self._landscape_invalidated_callback = callback

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
        self.sidebar.stand_controls.raster_mode.currentIndexChanged.connect(self._on_raster_mode_changed)
        self.sidebar.stand_controls.dbh_cutoff.valueChanged.connect(self._on_cutoff_changed)
        self.sidebar.stand_controls.dbh_cutoff.valueFinalized.connect(self._on_cutoff_finalized)
        self.sidebar.stand_controls.show_treatment_button.toggled.connect(self._on_show_treatment_toggled)
        self.pages.currentChanged.connect(self._on_page_changed)

    def _restore_state(self) -> None:
        self.sidebar.stand_controls.raster_mode.setCurrentIndex(self.state.raster_mode)
        self.sidebar.stand_controls.dbh_cutoff.setValue(int(round(dbh_to_display(self.state.dbh_cutoff, self.unit_system))))
        self.sidebar.stand_controls.show_treatment_button.setChecked(self.state.show_treatment)
        self.pages.setCurrentIndex(self.state.active_page)

        if self.rx_units:
            self.select_unit(self.state.selected_unit_index)
        self.refresh_all(trigger="restore_state")

    def _on_unit_changed(self, current: object, previous: object) -> None:
        del previous
        row = current.row() if getattr(current, "isValid", lambda: False)() else 0
        self.state.selected_unit_index = max(0, row)
        self._mark_all_pages_dirty()
        self.refresh_all(trigger="unit_changed")
        self._notify_state_changed("unit_changed")

    def _on_targets_changed(self) -> None:
        self.sidebar.model.refresh()
        self._mark_all_pages_dirty()
        self.refresh_all(trigger="targets_changed")
        self._notify_landscape_invalidated()
        self._notify_state_changed("targets_changed")

    def _on_raster_mode_changed(self, index: int) -> None:
        self.state.raster_mode = index
        self._raster_dirty = True
        self.refresh_current_page(trigger="raster_mode_changed")
        self._notify_state_changed("raster_mode_changed")

    def _on_cutoff_changed(self, value: int) -> None:
        self._pending_dbh_cutoff = dbh_from_display(float(value), self.unit_system)
        self.state.dbh_cutoff = self._pending_dbh_cutoff
        self._treatment_report_dirty = True
        self._dbh_update_timer.start()
        self._notify_state_changed("cutoff_changed")

    def _on_cutoff_finalized(self, value: int) -> None:
        self._pending_dbh_cutoff = dbh_from_display(float(value), self.unit_system)
        self.state.dbh_cutoff = self._pending_dbh_cutoff
        self._flush_cutoff_change(trigger="cutoff_finalized")

    def _flush_cutoff_change(self, trigger: str = "cutoff_debounced") -> None:
        self._dbh_update_timer.stop()
        self.state.dbh_cutoff = self._pending_dbh_cutoff
        self.refresh_current_page(trigger=trigger)

    def _on_show_treatment_toggled(self, checked: bool) -> None:
        if checked and self.rx_units:
            unit = self.current_unit()
            self.treatment_engine.do_treatment(unit, self.dbh_min(), self.dbh_max())
            self._treated_unit_ids.add(id(unit))
            self.sidebar.model.refresh()
            self._notify_landscape_invalidated()
        self.state.show_treatment = checked
        self._mark_all_pages_dirty()
        self.refresh_all(trigger="show_treatment_toggled")
        self._notify_state_changed("show_treatment_toggled")

    def _on_page_changed(self, index: int) -> None:
        self.state.active_page = index
        self.refresh_current_page(trigger="page_changed")
        self._notify_state_changed("page_changed")

    def refresh_all(self, trigger: str = "refresh_all") -> None:
        if not self.rx_units:
            self._mark_all_pages_dirty()
            self.refresh_current_page(trigger=trigger)
            return
        self.sidebar.structure_info.update_for_unit(self.current_unit())
        self.refresh_current_page(trigger=trigger)

    def refresh_current_page(self, trigger: str = "refresh_current_page") -> None:
        if self.state.active_page == 1:
            if not self._treatment_report_dirty:
                return
            self.update_treatment_report(trigger=trigger)
        elif self.state.active_page == 2:
            if not self._cut_report_dirty:
                return
            self.update_cut_report(trigger=trigger)
        else:
            if not self._raster_dirty:
                return
            self.update_raster_canvas(trigger=trigger)

    def _mark_all_pages_dirty(self) -> None:
        self._raster_dirty = True
        self._treatment_report_dirty = True
        self._cut_report_dirty = True

    def update_raster_canvas(self, trigger: str = "update_raster_canvas") -> None:
        self._update_raster_canvas_impl()
        self._raster_dirty = False

    def _update_raster_canvas_impl(self) -> None:
        figure = self.visualize_tab.raster_figure
        axes = self.visualize_tab.raster_axes
        self.visualize_tab.raster_axes = axes
        self.visualize_tab.sync_raster_colorbar_axes()
        axes.cla()
        self._raster_image = None
        self._hillshade_image = None
        self._raster_text = None

        if not self.rx_units:
            axes.text(0.5, 0.5, "No units available", ha="center", va="center")
            self._clear_raster_colorbar()
        else:
            unit = self.current_unit()
            data = self.current_raster_array()
            hillshade = unit.get_treat_hillshade() if self.state.show_treatment else unit.get_hillshade()
            colorbar_label = None
            colorbar_ticks = None
            colorbar_ticklabels = None

            if data.size == 0:
                self._raster_text = axes.text(0.5, 0.5, "No raster data", ha="center", va="center")
                self._clear_raster_colorbar()

            elif self.state.raster_mode == 0:
                display_data = array_to_display("canopy_height", data, self.unit_system)
                self._raster_image = axes.imshow(display_data, cmap="coolwarm", vmin=0)
                if hillshade.size:
                    self._hillshade_image = axes.imshow(hillshade, cmap="Greys", alpha=0.5)
                colorbar_label = display_name_for("canopy_height", self.unit_system)
                self._set_dynamic_raster_limits(display_data, lower_bound=0.0)
            elif self.state.raster_mode == 1:
                basin_display = self._categorical_basin_display(data)
                self._raster_image = axes.imshow(
                    basin_display,
                    cmap=self._basin_colormap,
                    interpolation="nearest",
                    vmin=0,
                    vmax=self._BASIN_PALETTE_SIZE - 1,
                )
                if hillshade.size:
                    self._hillshade_image = axes.imshow(hillshade, cmap="Greys", alpha=0.5)
                self._clear_raster_colorbar()
            elif self.state.raster_mode == 2:
                colors = ("white", "#7bc043", "#fdf498", "#f37736", "#ee4035")
                cmap = LinearSegmentedColormap.from_list("Clump Colors", colors, 5)
                boundaries = (-0.5, 0.5, 1.5, 4.5, 9.5, 99)
                norm = BoundaryNorm(boundaries, len(boundaries))
                self._raster_image = axes.imshow(data, cmap=cmap, norm=norm)
                if hillshade.size:
                    self._hillshade_image = axes.imshow(hillshade, cmap="Greys", alpha=0.5)
                colorbar_label = "Clump Map (Clump bins)"
                colorbar_ticks = [0, 1, 3, 7, 55]
                colorbar_ticklabels = ["0", "1", "2-4", "4-9", "10+"]
            else:
                self._raster_image = axes.imshow(data, cmap="Greys")
                colorbar_label = "Hillshade"
                self._set_dynamic_raster_limits(data)

            suffix = " (Treated)" if self.state.show_treatment else ""
            axes.set_title(f"{unit.name} {self.sidebar.stand_controls.raster_mode.currentText()}{suffix}")
            if self.state.raster_mode == 1:
                pass
            elif self._raster_image is not None and colorbar_label is not None:
                self._update_raster_colorbar(colorbar_label, colorbar_ticks, colorbar_ticklabels)

        axes.set_xticks([])
        axes.set_yticks([])
        self.visualize_tab.raster_canvas.draw_idle()

    def update_treatment_report(self, trigger: str = "update_treatment_report") -> None:
        del trigger
        self._update_treatment_report_impl()
        self._treatment_report_dirty = False

    def _update_treatment_report_impl(self) -> None:
        report = self.treatment_report_tab
        report.current_ba_axes.clear()
        report.current_mcs_axes.clear()
        report.displayed_ba_axes.clear()
        report.displayed_mcs_axes.clear()

        if not self.rx_units:
            report.current_label.setText("Current\n-")
            report.displayed_label.setText("Post-Treatment\n-")
            report.target_label.setText("Target\n-")
            report.displayed_mcs_prop.setText("")
            for canvas in (
                report.current_ba_canvas,
                report.current_mcs_canvas,
                report.displayed_ba_canvas,
                report.displayed_mcs_canvas,
            ):
                canvas.draw_idle()
            return

        unit = self.current_unit()
        has_treatment = self._has_treatment_results(unit)
        current_points = unit.get_taos()
        current_ba = self._tree_ba_distribution(current_points)
        current_csd = self._safe_array(unit.get_clump_sizes())

        if has_treatment:
            treated_points = unit.get_treat_taos()
            treated_ba = self._tree_ba_distribution(treated_points)
            treated_csd = self._safe_array(unit.get_treat_clump_sizes())
        else:
            treated_ba = current_ba
            treated_csd = current_csd

        report.current_label.setText(
            "Current\n"
            f"{self._format_structure_summary(unit.currentStructure, self.unit_system)}"
        )
        report.displayed_label.setText(
            "Post-Treatment\n"
            f"{self._format_structure_summary(unit.treatedStructure, self.unit_system)}"
        )
        report.target_label.setText(
            "Target\n"
            f"{self._format_structure_summary(unit.targetStructure, self.unit_system)}"
        )

        self._configure_density_axes(
            report.current_ba_axes,
            f"Basal Area ({self._density_ba_unit_label()})",
            "Kernel Density",
            "Pre-Treatment Tree Basal Area Distribution",
        )
        self._configure_density_axes(
            report.displayed_ba_axes,
            f"Basal Area ({self._density_ba_unit_label()})",
            "Kernel Density",
            "Post-Treatment Tree Basal Area Distribution",
        )
        self._configure_density_axes(
            report.current_mcs_axes,
            "Clump Size (n trees)",
            "Kernel Density",
            "Pre-Treatment Clump Size Distribution",
        )
        self._configure_density_axes(
            report.displayed_mcs_axes,
            "Clump Size (n trees)",
            "Kernel Density",
            "Post-Treatment Clump Size Distribution",
        )

        self._draw_density(report.current_ba_axes, current_ba)
        self._draw_density(report.current_mcs_axes, current_csd)
        report.displayed_mcs_prop.setText("")

        self._draw_density(report.displayed_ba_axes, treated_ba)
        self._draw_density(report.displayed_mcs_axes, treated_csd)
        self._sync_axis_limits(report.current_ba_axes, report.displayed_ba_axes)
        self._sync_axis_limits(report.current_mcs_axes, report.displayed_mcs_axes)
        report.displayed_mcs_prop.setText(self._mcs_prop_text(unit.treatedStructure))

        report.current_ba_canvas.draw_idle()
        report.current_mcs_canvas.draw_idle()
        report.displayed_ba_canvas.draw_idle()
        report.displayed_mcs_canvas.draw_idle()

    def update_cut_report(self, trigger: str = "update_cut_report") -> None:
        del trigger
        self._update_cut_report_impl()
        self._cut_report_dirty = False

    def _has_treatment_results(self, unit: RxUnit) -> bool:
        return id(unit) in self._treated_unit_ids

    def _update_cut_report_impl(self) -> None:
        report = self.cut_report_tab
        report.cut_axes.clear()

        if not self.rx_units:
            report.cut_summary.setText("Cut Trees:\nNo units available.")
            report.cut_canvas.draw_idle()
            return

        unit = self.current_unit()
        cut_points = unit.get_cut_taos()
        self._configure_density_axes(
            report.cut_axes,
            f"Basal Area ({self._density_ba_unit_label()})",
            "Kernel Density",
            "Cut Trees Basal Area Distribution",
        )

        if cut_points.ndim != 2 or cut_points.shape[0] == 0:
            report.cut_summary.setText("Cut Trees:\nBA:")
            report.cut_axes.text(0.5, 0.5, "No cut-tree data yet", ha="center", va="center")
        else:
            cut_ba = self._tree_ba_distribution(cut_points)
            self._draw_density(report.cut_axes, cut_ba)
            ba_per_area = self._cut_ba_per_area(cut_ba, unit.areaHa)
            report.cut_summary.setText(
                f"Cut Trees:\n\nBA:\t{ba_per_area:.3f} {label_for('ba', self.unit_system)}"
            )

        report.cut_canvas.draw_idle()

    def _notify_state_changed(self, reason: str) -> None:
        if self._state_changed_callback is not None:
            self._state_changed_callback(reason)

    def _notify_landscape_invalidated(self) -> None:
        if self._landscape_invalidated_callback is not None:
            self._landscape_invalidated_callback()

    def _clear_raster_colorbar(self) -> None:
        if self._raster_colorbar is not None:
            self.visualize_tab.sync_raster_colorbar_axes()
            self.visualize_tab.raster_colorbar_axes.cla()
            self.visualize_tab.raster_colorbar_axes.set_visible(False)
            self._raster_colorbar = None
        self._raster_mode_key = None
        self._raster_scale_key = None

    def _update_raster_colorbar(
        self,
        label: str,
        ticks: list[float] | None = None,
        ticklabels: list[str] | None = None,
    ) -> None:
        if self._raster_image is None:
            self._clear_raster_colorbar()
            return

        mode_key = (self.state.raster_mode, self.state.show_treatment)
        clim = self._raster_image.get_clim()
        scale_key = (clim[0], clim[1])
        needs_rebuild = (
            self._raster_colorbar is None
            or self._raster_mode_key != mode_key
            or self._raster_scale_key != scale_key
        )
        if needs_rebuild:
            self._clear_raster_colorbar()
            self.visualize_tab.sync_raster_colorbar_axes()
            self._raster_colorbar = self.visualize_tab.raster_figure.colorbar(
                self._raster_image,
                cax=self.visualize_tab.raster_colorbar_axes,
                orientation="vertical",
                label=label,
                ticks=ticks,
            )
            self.visualize_tab.raster_colorbar_axes.set_visible(True)
        else:
            self._raster_colorbar.update_normal(self._raster_image)
            self._raster_colorbar.set_label(label)
            if ticks is not None:
                self._raster_colorbar.set_ticks(ticks)

        self._raster_mode_key = mode_key
        self._raster_scale_key = scale_key
        if self._raster_colorbar is not None:
            if ticklabels is not None:
                self._raster_colorbar.ax.set_yticklabels(ticklabels)
            self._raster_colorbar.ax.tick_params(labelsize=8)

    def _set_dynamic_raster_limits(self, data: np.ndarray, lower_bound: float | None = None) -> None:
        if self._raster_image is None:
            return
        min_value, max_value = self._array_bounds(data)
        if lower_bound is not None:
            min_value = lower_bound
        if min_value is None or max_value is None:
            min_value, max_value = (0.0, 1.0)
        if max_value <= min_value:
            max_value = min_value + 1.0
        self._raster_image.set_clim(float(min_value), float(max_value))

    @classmethod
    def _build_basin_colormap(cls) -> ListedColormap:
        hues = (np.arange(cls._BASIN_PALETTE_SIZE, dtype=float) * 0.6180339887498949) % 1.0
        saturations = np.where(np.arange(cls._BASIN_PALETTE_SIZE) % 2 == 0, 0.70, 0.55)
        values = np.where(np.arange(cls._BASIN_PALETTE_SIZE) % 4 < 2, 0.92, 0.78)
        hsv = np.column_stack((hues, saturations, values))
        colors = hsv_to_rgb(hsv)
        colors[0] = (0.65, 0.65, 0.65)
        cmap = ListedColormap(colors, name="Basin Colors")
        cmap.set_bad((1.0, 1.0, 1.0, 1.0))
        return cmap

    def _categorical_basin_display(self, data: np.ndarray) -> np.ma.MaskedArray:
        array = np.asarray(data)
        if array.size == 0:
            return np.ma.masked_array(array, mask=np.ones_like(array, dtype=bool))

        if not np.issubdtype(array.dtype, np.integer):
            integer_array = array.astype(np.int64, copy=False)
        else:
            integer_array = array

        dtype_sentinel = np.iinfo(integer_array.dtype).min
        int32_sentinel = np.int64(np.iinfo(np.int32).min)
        mask = (integer_array == dtype_sentinel) | (integer_array == int32_sentinel)
        hashed = np.zeros(integer_array.shape, dtype=np.uint8)
        removed_mask = (~mask) & (integer_array == 1) if self.state.show_treatment else np.zeros_like(mask, dtype=bool)
        valid_tree_mask = (~mask) & (~removed_mask)
        if np.any(valid_tree_mask):
            ids = integer_array[valid_tree_mask].astype(np.uint64, copy=False)
            ids ^= ids >> np.uint64(33)
            ids *= np.uint64(0xff51afd7ed558ccd)
            ids ^= ids >> np.uint64(33)
            ids *= np.uint64(0xc4ceb9fe1a85ec53)
            ids ^= ids >> np.uint64(33)
            hashed[valid_tree_mask] = (np.mod(ids, self._BASIN_PALETTE_SIZE - 1) + 1).astype(np.uint8, copy=False)
        hashed[removed_mask] = 0

        return np.ma.masked_array(hashed, mask=mask)

    @staticmethod
    def _array_bounds(data: np.ndarray) -> tuple[float | None, float | None]:
        array = np.asarray(data)
        if array.size == 0:
            return (None, None)
        if np.issubdtype(array.dtype, np.integer):
            sentinel = np.iinfo(array.dtype).min
            array = array[array != sentinel]
        else:
            array = array[np.isfinite(array)]
        if array.size == 0:
            return (None, None)
        return (float(np.min(array)), float(np.max(array)))

    @staticmethod
    def _configure_density_axes(
        axes: object,
        xlabel: str,
        ylabel: str,
        title: str,
    ) -> None:
        axes.cla()
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        axes.set_title(title)

    @staticmethod
    def _safe_array(values: object) -> np.ndarray:
        array = np.asarray(values, dtype=float)
        return array[np.isfinite(array)]

    def _tree_ba_distribution(self, points: object) -> np.ndarray:
        point_array = np.asarray(points, dtype=float)
        if point_array.ndim != 2 or point_array.shape[0] == 0 or point_array.shape[1] < 5:
            return np.zeros(0, dtype=float)

        dbh_cm = point_array[:, 4]
        valid_dbh = dbh_cm[np.isfinite(dbh_cm) & (dbh_cm > 0.0)]
        if valid_dbh.size == 0:
            return np.zeros(0, dtype=float)

        dbh_m = valid_dbh / 100.0
        ba_m2 = np.pi * np.square(dbh_m / 2.0)
        if self.unit_system.value == "imperial":
            return ba_m2 * FEET_PER_METER * FEET_PER_METER
        return ba_m2

    def _cut_ba_per_area(self, cut_ba: np.ndarray, area_ha: float) -> float:
        if cut_ba.size == 0 or not np.isfinite(area_ha) or area_ha <= 0.0:
            return 0.0
        if self.unit_system.value == "imperial":
            area_acres = area_ha / (4046.8564224 / 10000.0)
            if area_acres <= 0.0:
                return 0.0
            return float(np.sum(cut_ba) / area_acres)
        return float(np.sum(cut_ba) / area_ha)

    @staticmethod
    def _draw_density(axes: object, values: np.ndarray) -> None:
        if values.size == 0:
            axes.text(0.5, 0.5, "No data available", ha="center", va="center", transform=axes.transAxes)
            return
        if values.size == 1 or np.allclose(values, values[0]):
            axes.axvline(values[0], color="#4a6fa5", linewidth=2)
            return

        density = gaussian_kde(values)
        density.covariance_factor = lambda: 0.25
        density._compute_covariance()
        x_values = np.linspace(0.0, float(np.max(values)), 200)
        axes.plot(x_values, density(x_values), color="#2f5d8a", linewidth=2)

    @staticmethod
    def _sync_axis_limits(primary_axes: object, secondary_axes: object) -> None:
        x_range = (
            min(primary_axes.get_xlim()[0], secondary_axes.get_xlim()[0]),
            max(primary_axes.get_xlim()[1], secondary_axes.get_xlim()[1]),
        )
        y_range = (
            min(primary_axes.get_ylim()[0], secondary_axes.get_ylim()[0]),
            max(primary_axes.get_ylim()[1], secondary_axes.get_ylim()[1]),
        )
        primary_axes.set_xlim(x_range)
        primary_axes.set_ylim(y_range)
        secondary_axes.set_xlim(x_range)
        secondary_axes.set_ylim(y_range)

    def _mcs_prop_text(self, structure: StructureSummary) -> str:
        if self._mcs_prop_table is None:
            return ""
        rows = self._mcs_prop_table
        mask = np.logical_and(rows["MCS_min"] < structure.mcs, rows["MCS_max"] >= structure.mcs)
        if not np.any(mask):
            return ""
        row = rows[mask][0]
        fields = ("1", "2to4", "5to9", "9to14", "15to30", "gt30")
        props = "\t".join(str(row[field]) for field in fields)
        return "Clumps:  1\t2-4\t5-9\t10-14\t15-29\t30+\nProps:    " + props

    @staticmethod
    def _load_mcs_prop_table() -> np.ndarray | None:
        candidates = (
            Path(__file__).resolve().parents[2] / "resources" / "mcs_prop.csv",
            Path(__file__).resolve().parents[2] / "mcs_prop.csv",
        )
        for candidate in candidates:
            if candidate.exists():
                return np.genfromtxt(candidate, delimiter=",", names=True, dtype=None, encoding="utf-8")
        return None

    def _density_ba_unit_label(self) -> str:
        return "ft^2/ac" if self.unit_system.value == "imperial" else "m^2/ha"

    @staticmethod
    def _format_structure_summary(
        structure: StructureSummary | None,
        unit_system: object,
        empty_text: str = "-",
    ) -> str:
        if structure is None:
            return empty_text
        return (
            f"{label_for('tph', unit_system)}: {format_value('tph', structure.tph, unit_system)}\n"
            f"BA: {format_value('ba', structure.ba, unit_system)}\n"
            f"MCS: {format_value('mcs', structure.mcs, unit_system)}\n"
            f"CC: {format_value('cc', structure.cc, unit_system)}"
        )
