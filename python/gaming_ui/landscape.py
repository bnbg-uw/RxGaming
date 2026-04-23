from __future__ import annotations

import csv
from collections.abc import Callable
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from PySide6.QtWidgets import QVBoxLayout, QWidget

from rxgaming_core import ProjectSettings, RxUnit
from .state import StandViewState
from .units import array_to_display, display_name_for


class LandscapeReferenceTab(QWidget):
    def __init__(self, rx_units: list[RxUnit], project_settings: ProjectSettings, state: StandViewState) -> None:
        super().__init__()
        self._rx_units = rx_units
        self._state = state
        self._unit_system = state.unit_system
        self._on_unit_selected: Callable[[int], None] | None = None
        self._reference = self._load_reference_points(project_settings.refDataPath)

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ba_ax = self.figure.add_subplot(1, 3, 1)
        self.mcs_ax = self.figure.add_subplot(1, 3, 2)
        self.cc_ax = self.figure.add_subplot(1, 3, 3)

        self._axes = [self.ba_ax, self.mcs_ax, self.cc_ax]
        self._metric_names = ["ba", "mcs", "cc"]
        self._metric_labels = ["Basal Area", "Mean Clump Size", "Canopy Cover"]
        self._unit_lines: list[Line2D] = []
        self._annotations = []

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)

        self.canvas.mpl_connect("motion_notify_event", self._on_hover)
        self.canvas.mpl_connect("button_press_event", self._on_click)
        self.refresh()

    def set_unit_selected_callback(self, callback: Callable[[int], None]) -> None:
        self._on_unit_selected = callback

    def refresh(self) -> None:
        self.figure.legends.clear()
        for axis in self._axes:
            axis.clear()

        ref_tph = self.array_to_display(np.asarray(self._reference.get("tph", []), dtype=float), "tph")
        ref_series = {
            "ba": self.array_to_display(np.asarray(self._reference.get("ba", []), dtype=float), "ba"),
            "mcs": np.asarray(self._reference.get("mcs", []), dtype=float),
            "cc": np.asarray(self._reference.get("cc", []), dtype=float),
        }
        unit_tph = self.array_to_display(
            np.asarray([unit.currentStructure.tph for unit in self._rx_units], dtype=float), "tph"
        )

        self._unit_lines = []
        self._annotations = []
        if not self._rx_units:
            for axis, metric_name, metric_label in zip(self._axes, self._metric_names, self._metric_labels):
                axis.set_title(metric_label)
                axis.set_xlabel(display_name_for("tph", self._unit_system))
                axis.set_ylabel(display_name_for(metric_name, self._unit_system))
                axis.text(0.5, 0.5, "No units available", ha="center", va="center")
            self.figure.tight_layout()
            self.canvas.draw_idle()
            return

        for axis, metric_name, metric_label in zip(
            self._axes, self._metric_names, self._metric_labels
        ):
            series = ref_series[metric_name]
            if ref_tph.size and series.size:
                axis.scatter(ref_tph, series, color="#f28263", alpha=0.15, s=12, label="Reference")

            unit_values = np.asarray(
                [getattr(unit.currentStructure, metric_name) for unit in self._rx_units],
                dtype=float,
            )
            unit_values = self.array_to_display(unit_values, metric_name)
            unit_line = axis.plot(unit_tph, unit_values, "b^", label="Units")[0]
            self._unit_lines.append(unit_line)

            selected = self._selected_unit()
            target = selected.targetStructure
            axis.plot(
                [self.array_to_display(np.asarray([target.tph], dtype=float), "tph")[0]],
                [self.array_to_display(np.asarray([getattr(target, metric_name)], dtype=float), metric_name)[0]],
                "gs",
                label="Target",
            )
            if selected.treatedStructure is not None:
                treated = selected.treatedStructure
                axis.plot(
                    [self.array_to_display(np.asarray([treated.tph], dtype=float), "tph")[0]],
                    [self.array_to_display(np.asarray([getattr(treated, metric_name)], dtype=float), metric_name)[0]],
                    "mo",
                    label="Treated",
                )
                axis.annotate(
                    "",
                    xy=(
                        self.array_to_display(np.asarray([treated.tph], dtype=float), "tph")[0],
                        self.array_to_display(np.asarray([getattr(treated, metric_name)], dtype=float), metric_name)[0],
                    ),
                    xytext=(
                        self.array_to_display(np.asarray([target.tph], dtype=float), "tph")[0],
                        self.array_to_display(np.asarray([getattr(target, metric_name)], dtype=float), metric_name)[0],
                    ),
                    arrowprops={"arrowstyle": "->", "color": "#884ea0"},
                )

            axis.annotate(
                "",
                xy=(
                    self.array_to_display(np.asarray([target.tph], dtype=float), "tph")[0],
                    self.array_to_display(np.asarray([getattr(target, metric_name)], dtype=float), metric_name)[0],
                ),
                xytext=(
                    self.array_to_display(np.asarray([selected.currentStructure.tph], dtype=float), "tph")[0],
                    self.array_to_display(
                        np.asarray([getattr(selected.currentStructure, metric_name)], dtype=float), metric_name
                    )[0],
                ),
                arrowprops={"arrowstyle": "->", "color": "#2e7d32"},
            )

            axis.set_title(metric_label)
            axis.set_xlabel(display_name_for("tph", self._unit_system))
            axis.set_ylabel(display_name_for(metric_name, self._unit_system))

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

        max_mcs = max(
            [float(np.max(ref_series["mcs"])) if ref_series["mcs"].size else 0.0]
            + [float(unit.currentStructure.mcs) for unit in self._rx_units]
        )
        if max_mcs >= 10:
            self.mcs_ax.set_yscale("log")
            self.mcs_ax.set_ylim(bottom=1)

        legend_handles = [
            Line2D([], [], color="#f28263", marker="o", linestyle="None", alpha=0.4, label="Reference"),
            Line2D([], [], color="blue", marker="^", linestyle="None", label="Units"),
            Line2D([], [], color="green", marker="s", linestyle="None", label="Targets"),
            Line2D([], [], color="magenta", marker="o", linestyle="None", label="Treated"),
        ]
        self.figure.legend(legend_handles, ["Reference", "Units", "Targets", "Treated"], loc="right")
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _selected_unit(self) -> RxUnit:
        if not self._rx_units:
            raise IndexError("No units available")
        index = max(0, min(self._state.selected_unit_index, len(self._rx_units) - 1))
        return self._rx_units[index]

    def _on_hover(self, event: object) -> None:
        inaxes = getattr(event, "inaxes", None)
        for annotation in self._annotations:
            annotation.set_visible(False)
        if inaxes not in self._axes:
            self.canvas.draw_idle()
            return

        axis_index = self._axes.index(inaxes)
        contains, info = self._unit_lines[axis_index].contains(event)
        if not contains:
            self.canvas.draw_idle()
            return

        point_index = info["ind"][0]
        unit = self._rx_units[point_index]
        x_data = self._unit_lines[axis_index].get_xdata()
        y_data = self._unit_lines[axis_index].get_ydata()
        annotation = self._annotations[axis_index]
        annotation.xy = (x_data[point_index], y_data[point_index])
        annotation.set_text(unit.name)
        annotation.set_visible(True)
        self.canvas.draw_idle()

    def _on_click(self, event: object) -> None:
        if getattr(event, "button", None) != 1:
            return
        for axis_index, line in enumerate(self._unit_lines):
            contains, info = line.contains(event)
            if contains:
                point_index = info["ind"][0]
                if self._on_unit_selected is not None:
                    self._on_unit_selected(point_index)
                self._annotations[axis_index].set_visible(False)
                break

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

    def array_to_display(self, values: np.ndarray, metric_kind: str) -> np.ndarray:
        return array_to_display(metric_kind, values, self._unit_system)
