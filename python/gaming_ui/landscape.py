from __future__ import annotations

import csv
from collections.abc import Callable
from pathlib import Path
from typing import Any, TYPE_CHECKING, cast

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Polygon as MplPolygon
import seaborn
from scipy.spatial import ConvexHull
from PySide6.QtCore import QTimer
from PySide6.QtGui import QResizeEvent, QShowEvent
from PySide6.QtWidgets import QVBoxLayout, QWidget

from rxgaming_core import ProjectSettings, RxUnit
from .state import StandViewState
from .units import array_to_display, display_name_for, MetricKind

from shapely.geometry import GeometryCollection, MultiPolygon, Polygon as ShapelyPolygon

if TYPE_CHECKING:
    from matplotlib.backend_bases import MouseEvent, Event
class LandscapeReferenceTab(QWidget):
    _REFERENCE_COLOR = "#f28263"
    _REFERENCE_FILL_ALPHA = 0.35
    _HULL_FILL_ALPHA = 0.32
    _INTERSECTION_COLOR = "#657fb5"
    _INTERSECTION_FILL_ALPHA = 0.45
    _MCS_LOG_THRESHOLD = 10.0

    def __init__(self, rx_units: list[RxUnit], project_settings: ProjectSettings, state: StandViewState) -> None:
        super().__init__()
        self._rx_units = rx_units
        self._state = state
        self._unit_system = state.unit_system
        self._on_unit_selected: Callable[[int], None] | None = None
        self._metric_names : list[MetricKind] = ["ba", "mcs", "cc"]
        self._metric_labels = ["Basal Area", "Mean Clump Size", "Canopy Cover"]
        self._reference = self._load_reference_points(project_settings.refDataPath)
        self._decision_spaces = self._load_decision_spaces()
        self._cached_geometry = self._build_cached_geometry()

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ba_ax = self.figure.add_subplot(1, 3, 1)
        self.mcs_ax = self.figure.add_subplot(1, 3, 2)
        self.cc_ax = self.figure.add_subplot(1, 3, 3)

        self._axes = [self.ba_ax, self.mcs_ax, self.cc_ax]
        self._unit_lines: list[Line2D] = []
        self._target_markers: list[Line2D] = []
        self._treated_markers: list[Line2D] = []
        self._annotations: list[Any] = []
        self._current_to_target_arrows: list[FancyArrowPatch] = []
        self._target_to_treated_arrows: list[FancyArrowPatch] = []
        self._unit_hull_patches: list[list[list[MplPolygon]]] = []
        self._intersection_patches: list[list[list[MplPolygon]]] = []
        self._empty_texts: list[Any] = []
        self._awaiting_post_layout_draw = True

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)

        self.canvas.mpl_connect("motion_notify_event", self._on_hover)
        self.canvas.mpl_connect("button_press_event", self._on_click)
        self._initialize_artists()
        self._render_reference_backgrounds()
        self._render_decision_space_patches()
        self._sync_legend()

    def set_unit_selected_callback(self, callback: Callable[[int], None]) -> None:
        self._on_unit_selected = callback

    def refresh(self, trigger: str = "refresh") -> None:
        if self._awaiting_post_layout_draw:
            if not self._is_layout_ready():
                return
            self._awaiting_post_layout_draw = False
        self._refresh_impl()

    def showEvent(self, event: QShowEvent) -> None:
        super().showEvent(event)
        if self._awaiting_post_layout_draw and self._is_layout_ready():
            QTimer.singleShot(0, self._perform_post_layout_draw)

    def resizeEvent(self, event: QResizeEvent) -> None:
        super().resizeEvent(event)
        if self._awaiting_post_layout_draw and event.size().width() > 0 and event.size().height() > 0 and self.isVisible():
            QTimer.singleShot(0, self._perform_post_layout_draw)

    def _refresh_impl(self) -> None:
        if not self._rx_units:
            for axis, metric_name, metric_label in zip(self._axes, self._metric_names, self._metric_labels):
                axis.set_title(metric_label)
                axis.set_xlabel(display_name_for("tph", self._unit_system))
                axis.set_ylabel(display_name_for(metric_name, self._unit_system))
            self._show_no_units_state()
            self.figure.tight_layout(rect=(0.0, 0.0, 0.92, 1.0))
            self.canvas.draw_idle()
            return

        self._hide_no_units_state()
        unit_tph = self.array_to_display(np.asarray([unit.currentStructure.tph for unit in self._rx_units], dtype=float), "tph")

        for metric_index, (axis, metric_name, metric_label) in enumerate(zip(self._axes, self._metric_names, self._metric_labels)):
            unit_values = self.array_to_display(
                np.asarray([getattr(unit.currentStructure, metric_name) for unit in self._rx_units], dtype=float),
                metric_name,
            )
            self._unit_lines[metric_index].set_data(unit_tph, unit_values)
            axis.set_title(metric_label)
            axis.set_xlabel(display_name_for("tph", self._unit_system))
            axis.set_ylabel(display_name_for(metric_name, self._unit_system))
            axis.relim()
            axis.autoscale_view()

        self._hide_hover_state()

        max_mcs = max([self._cached_geometry["reference_mcs_max"]] + [float(unit.currentStructure.mcs) for unit in self._rx_units])
        if max_mcs >= self._MCS_LOG_THRESHOLD:
            self.mcs_ax.set_yscale("log")
            self.mcs_ax.set_ylim(bottom=1)
        else:
            self.mcs_ax.set_yscale("linear")
        self.figure.tight_layout(rect=(0.0, 0.0, 0.92, 1.0))
        self.canvas.draw_idle()

    def _perform_post_layout_draw(self) -> None:
        if not self._awaiting_post_layout_draw:
            return
        if not self._is_layout_ready():
            return
        self.refresh(trigger="post_layout_draw")
        self.canvas.draw()

    def _is_layout_ready(self) -> bool:
        return self.isVisible() and self.width() > 0 and self.height() > 0 and self.canvas.width() > 0 and self.canvas.height() > 0

    def _initialize_artists(self) -> None:
        self._unit_lines = []
        self._target_markers = []
        self._treated_markers = []
        self._annotations = []
        self._current_to_target_arrows = []
        self._target_to_treated_arrows = []
        self._unit_hull_patches = []
        self._intersection_patches = []
        self._empty_texts = []

        for axis, metric_name, metric_label in zip(self._axes, self._metric_names, self._metric_labels):
            axis.set_title(metric_label)
            axis.set_xlabel(display_name_for("tph", self._unit_system))
            axis.set_ylabel(display_name_for(metric_name, self._unit_system))
            self._unit_lines.append(axis.plot([], [], "b^", label="Units")[0])
            self._target_markers.append(axis.plot([], [], "gs", label="Target")[0])
            self._treated_markers.append(axis.plot([], [], "mo", label="Treated")[0])
            annotation = axis.annotate(
                "",
                xy=(0.0, 0.0),
                xytext=(-20, 20),
                textcoords="offset points",
                bbox={"boxstyle": "round", "fc": "w"},
                arrowprops={"arrowstyle": "->"},
            )
            annotation.set_visible(False)
            self._annotations.append(annotation)
            current_to_target = FancyArrowPatch(
                (0.0, 0.0),
                (0.0, 0.0),
                arrowstyle="->",
                color="#2e7d32",
                mutation_scale=14,
                linewidth=2.0,
            )
            current_to_target.set_visible(False)
            current_to_target.set_zorder(6)
            axis.add_patch(current_to_target)
            self._current_to_target_arrows.append(current_to_target)
            target_to_treated = FancyArrowPatch(
                (0.0, 0.0),
                (0.0, 0.0),
                arrowstyle="->",
                color="#884ea0",
                mutation_scale=14,
                linewidth=2.0,
            )
            target_to_treated.set_visible(False)
            target_to_treated.set_zorder(6)
            axis.add_patch(target_to_treated)
            self._target_to_treated_arrows.append(target_to_treated)
            empty_text = axis.text(0.5, 0.5, "No units available", ha="center", va="center", transform=axis.transAxes)
            empty_text.set_visible(False)
            self._empty_texts.append(empty_text)
            self._unit_hull_patches.append([])
            self._intersection_patches.append([])

    def _sync_legend(self) -> None:
        self.figure.legends.clear()
        legend_handles = [
            Line2D([], [], color=self._REFERENCE_COLOR, linewidth=8, alpha=self._REFERENCE_FILL_ALPHA, label="Reference"),
            Line2D([], [], color="blue", marker="^", linestyle="None", label="Units"),
            Line2D([], [], color="green", marker="s", linestyle="None", label="Targets"),
            Line2D([], [], color="magenta", marker="o", linestyle="None", label="Treated"),
        ]
        self.figure.legend(legend_handles, ["Reference", "Units", "Targets", "Treated"], loc="right")

    def _render_reference_backgrounds(self) -> None:
        tph_values = self.array_to_display(np.asarray(self._reference.get("tph", []), dtype=float), "tph")
        for axis, metric_name in zip(self._axes, self._metric_names):
            metric_values = self.array_to_display(np.asarray(self._reference.get(metric_name, []), dtype=float), metric_name)
            valid_mask = np.isfinite(tph_values) & np.isfinite(metric_values)
            if np.count_nonzero(valid_mask) < 3:
                continue
            seaborn.kdeplot(
                ax=axis,
                x=tph_values[valid_mask],
                y=metric_values[valid_mask],
                cmap="Oranges",
                fill=True,
                bw_adjust=0.5,
                thresh=0.05,
                levels=8,
            )

    def _render_decision_space_patches(self) -> None:
        for metric_index in range(len(self._metric_names)):
            axis = self._axes[metric_index]
            for unit_geometry in self._cached_geometry["units"][metric_index]:
                hull_patches = self._create_polygon_patches(
                    axis,
                    unit_geometry["hull"],
                    facecolor=self._REFERENCE_COLOR,
                    edgecolor=self._REFERENCE_COLOR,
                    alpha=self._HULL_FILL_ALPHA,
                    zorder=3,
                )
                intersection_patches = self._create_polygon_patches(
                    axis,
                    unit_geometry["intersection"],
                    facecolor= self._INTERSECTION_COLOR,
                    edgecolor= self._INTERSECTION_COLOR,
                    alpha=self._INTERSECTION_FILL_ALPHA,
                    zorder=4,
                )
                self._unit_hull_patches[metric_index].append(hull_patches)
                self._intersection_patches[metric_index].append(intersection_patches)

    def _show_no_units_state(self) -> None:
        for empty_text in self._empty_texts:
            empty_text.set_visible(True)
        for unit_line, target, treated in zip(self._unit_lines, self._target_markers, self._treated_markers):
            unit_line.set_data([], [])
            target.set_data([], [])
            treated.set_data([], [])

    def _hide_no_units_state(self) -> None:
        for empty_text in self._empty_texts:
            empty_text.set_visible(False)

    def _selected_unit(self) -> RxUnit:
        if not self._rx_units:
            raise IndexError("No units available")
        index = max(0, min(self._state.selected_unit_index, len(self._rx_units) - 1))
        return self._rx_units[index]

    def _hide_hover_state(self) -> None:
        for annotation in self._annotations:
            annotation.set_visible(False)
        for target, treated, current_arrow, treated_arrow in zip(
            self._target_markers,
            self._treated_markers,
            self._current_to_target_arrows,
            self._target_to_treated_arrows,
        ):
            target.set_visible(False)
            treated.set_visible(False)
            current_arrow.set_visible(False)
            treated_arrow.set_visible(False)
        for metric_patch_groups in self._unit_hull_patches + self._intersection_patches:
            for patches in metric_patch_groups:
                for patch in patches:
                    patch.set_visible(False)
        for target_marker in self._target_markers:
            target_marker.set_data([], [])
        for treated_marker in self._treated_markers:
            treated_marker.set_data([], [])

    def _show_hover_state(self, unit_index: int, metric_index: int) -> None:
        self._hide_hover_state()
        unit = self._rx_units[unit_index]
        metric_name = self._metric_names[metric_index]
        x_value = self.array_to_display(np.asarray([unit.currentStructure.tph], dtype=float), "tph")[0]
        y_value = self.array_to_display(np.asarray([getattr(unit.currentStructure, metric_name)], dtype=float), metric_name)[0]
        annotation = self._annotations[metric_index]
        annotation.xy = (x_value, y_value)
        annotation.set_text(unit.name)
        annotation.set_visible(True)

        for patch in self._unit_hull_patches[metric_index][unit_index]:
            patch.set_visible(True)
        for patch in self._intersection_patches[metric_index][unit_index]:
            patch.set_visible(True)

        target = unit.targetStructure
        target_x = self.array_to_display(np.asarray([target.tph], dtype=float), "tph")[0]
        target_y = self.array_to_display(np.asarray([getattr(target, metric_name)], dtype=float), metric_name)[0]
        self._target_markers[metric_index].set_data([target_x], [target_y])
        self._target_markers[metric_index].set_visible(True)

        current_arrow = self._current_to_target_arrows[metric_index]
        current_arrow.set_positions((x_value, y_value), (target_x, target_y))
        current_arrow.set_visible(True)

        treated = unit.treatedStructure
        if treated is not None:
            treated_x = self.array_to_display(np.asarray([treated.tph], dtype=float), "tph")[0]
            treated_y = self.array_to_display(np.asarray([getattr(treated, metric_name)], dtype=float), metric_name)[0]
            self._treated_markers[metric_index].set_data([treated_x], [treated_y])
            self._treated_markers[metric_index].set_visible(True)
            treated_arrow = self._target_to_treated_arrows[metric_index]
            treated_arrow.set_positions((target_x, target_y), (treated_x, treated_y))
            treated_arrow.set_visible(True)

    def _on_hover(self, event: Event) -> None:
        event = cast("MouseEvent", event)
        inaxes = getattr(event, "inaxes", None)
        self._hide_hover_state()
        if inaxes not in self._axes:
            self.canvas.draw_idle()
            return

        axis_index = self._axes.index(inaxes)
        contains, info = self._unit_lines[axis_index].contains(event)
        if not contains:
            self.canvas.draw_idle()
            return

        point_index = info["ind"][0]
        self._show_hover_state(point_index, axis_index)
        self.canvas.draw_idle()

    def _on_click(self, event: Event) -> None:
        event = cast("MouseEvent", event)
        if getattr(event, "button", None) != 1:
            return
        for line in self._unit_lines:
            contains, info = line.contains(event)
            if contains:
                point_index = info["ind"][0]
                if self._on_unit_selected is not None:
                    self._on_unit_selected(point_index)
                self._hide_hover_state()
                self.canvas.draw_idle()
                break

    def _load_decision_spaces(self) -> list[list[Any]]:
        decision_spaces: list[list[Any]] = []
        for unit in self._rx_units:
            decision_spaces.append(list(unit.get_simulated_structures(self._state.dbh_cutoff)))
        return decision_spaces

    def _build_cached_geometry(self) -> dict[str, Any]:
        reference_tph = self.array_to_display(np.asarray(self._reference.get("tph", []), dtype=float), "tph")
        reference_metrics = {
            "ba": self.array_to_display(np.asarray(self._reference.get("ba", []), dtype=float), "ba"),
            "mcs": self.array_to_display(np.asarray(self._reference.get("mcs", []), dtype=float), "mcs"),
            "cc": self.array_to_display(np.asarray(self._reference.get("cc", []), dtype=float), "cc"),
        }
        reference_hulls = {
            metric_name: self._convex_hull_polygons(reference_tph, values)
            for metric_name, values in reference_metrics.items()
        }

        unit_geometry: list[list[dict[str, list[np.ndarray]]]] = []
        for metric_name in self._metric_names:
            metric_units: list[dict[str, list[np.ndarray]]] = []
            for structures in self._decision_spaces:
                tph_values = self.array_to_display(
                    np.asarray([structure.tph for structure in structures], dtype=float),
                    "tph",
                )
                metric_values = self.array_to_display(
                    np.asarray([getattr(structure, metric_name) for structure in structures], dtype=float),
                    metric_name,
                )
                hull = self._convex_hull_polygons(tph_values, metric_values)
                intersection = self._intersect_polygons(hull, reference_hulls[metric_name])
                metric_units.append({"hull": hull, "intersection": intersection})
            unit_geometry.append(metric_units)

        reference_mcs_max = 0.0
        if reference_metrics["mcs"].size:
            reference_mcs_max = float(np.nanmax(reference_metrics["mcs"]))
        return {
            "reference_hulls": reference_hulls,
            "reference_mcs_max": reference_mcs_max,
            "units": unit_geometry,
        }

    @staticmethod
    def _load_reference_points(path_text: str) -> dict[str, list[float]]:
        reference = {"tph": [], "ba": [], "mcs": [], "cc": []}
        path = Path(path_text)
        if not path.is_file():
            return reference

        with path.open("r", newline="", encoding="utf-8-sig") as fp:
            reader = csv.DictReader(fp)
            for row in reader:
                try:
                    reference["tph"].append(float(row["tph"]))
                    reference["ba"].append(float(row["ba"]))
                    reference["mcs"].append(float(row["mcs"]))
                    reference["cc"].append(float(row["cc"]))
                except (KeyError, TypeError, ValueError):
                    continue
        return reference

    @staticmethod
    def _finite_pairs(x_values: np.ndarray, y_values: np.ndarray) -> np.ndarray:
        pairs = np.column_stack((x_values, y_values))
        return pairs[np.isfinite(pairs).all(axis=1)]

    def _convex_hull_polygons(self, x_values: np.ndarray, y_values: np.ndarray) -> list[np.ndarray]:
        points = self._finite_pairs(x_values, y_values)
        if points.shape[0] < 3:
            return []
        unique_points = np.unique(points, axis=0)
        if unique_points.shape[0] < 3:
            return []
        try:
            hull = ConvexHull(unique_points)
        except Exception:
            return []
        return [unique_points[hull.vertices]]

    def _intersect_polygons(self, lhs: list[np.ndarray], rhs: list[np.ndarray]) -> list[np.ndarray]:
        if not lhs or not rhs or ShapelyPolygon is None:
            return []
        lhs_poly = self._to_shapely(lhs)
        rhs_poly = self._to_shapely(rhs)
        if lhs_poly is None or rhs_poly is None:
            return []
        intersection = lhs_poly.intersection(rhs_poly)
        return self._from_shapely(intersection)

    @staticmethod
    def _to_shapely(polygons: list[np.ndarray]) -> Any:
        if ShapelyPolygon is None:
            return None
        shapes = []
        for polygon in polygons:
            if polygon.shape[0] >= 3:
                shapes.append(ShapelyPolygon(polygon))
        if not shapes:
            return None
        geometry = shapes[0]
        for shape in shapes[1:]:
            geometry = geometry.union(shape)
        return geometry

    @staticmethod
    def _from_shapely(geometry: Any) -> list[np.ndarray]:
        if geometry is None or getattr(geometry, "is_empty", True):
            return []
        if ShapelyPolygon is not None and isinstance(geometry, ShapelyPolygon):
            return [np.asarray(geometry.exterior.coords[:-1], dtype=float)]
        if MultiPolygon is not None and isinstance(geometry, MultiPolygon):
            return [np.asarray(part.exterior.coords[:-1], dtype=float) for part in geometry.geoms if not part.is_empty]
        if GeometryCollection is not None and isinstance(geometry, GeometryCollection):
            polygons: list[np.ndarray] = []
            for part in geometry.geoms:
                polygons.extend(LandscapeReferenceTab._from_shapely(part))
            return polygons
        return []

    @staticmethod
    def _create_polygon_patches(
        axis: Any,
        polygons: list[np.ndarray],
        *,
        facecolor: str,
        edgecolor: str,
        alpha: float,
        zorder: int,
    ) -> list[MplPolygon]:
        patches: list[MplPolygon] = []
        for polygon in polygons:
            if polygon.shape[0] < 3:
                continue
            patch = MplPolygon(polygon, closed=True, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)
            patch.set_visible(False)
            patch.set_zorder(zorder)
            axis.add_patch(patch)
            patches.append(patch)
        return patches

    def array_to_display(self, values: np.ndarray, metric_kind: MetricKind) -> np.ndarray:
        return array_to_display(metric_kind, values, self._unit_system)
