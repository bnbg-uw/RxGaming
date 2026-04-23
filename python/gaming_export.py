from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from PySide6.QtWidgets import QFileDialog, QMessageBox, QWidget

from gaming_ui import GamingTabs


def export_current_view(tabs: GamingTabs, parent: QWidget) -> None:
    file_path = QFileDialog.getSaveFileName(parent, "Export current view", "", "PNG files (*.png)")[0]
    if not file_path:
        return

    output_path = Path(file_path)
    if output_path.suffix.lower() != ".png":
        output_path = output_path.with_suffix(".png")

    tabs.current_export_widget().grab().save(str(output_path), "PNG")


def export_raster(tabs: GamingTabs, parent: QWidget) -> None:
    file_path = QFileDialog.getSaveFileName(parent, "Export raster data", "", "NumPy files (*.npy)")[0]
    if not file_path:
        return
    output_path = Path(file_path)
    if output_path.suffix.lower() != ".npy":
        output_path = output_path.with_suffix(".npy")
    np.save(output_path, tabs.current_raster_array())


def export_features(tabs: GamingTabs, parent: QWidget) -> None:
    _export_points_csv(tabs, parent, "Export point data")


def export_treelist(tabs: GamingTabs, parent: QWidget) -> None:
    _export_points_csv(tabs, parent, "Export treelist")


def _export_points_csv(tabs: GamingTabs, parent: QWidget, title: str) -> None:
    file_path = QFileDialog.getSaveFileName(parent, title, "", "CSV files (*.csv)")[0]
    if not file_path:
        return

    output_path = Path(file_path)
    if output_path.suffix.lower() != ".csv":
        output_path = output_path.with_suffix(".csv")

    points = np.asarray(tabs.current_points())
    if points.ndim != 2 or points.shape[1] == 0:
        _show_warning(parent, "No point data is available to export.")
        return

    headers = ["x", "y", "area", "height", "crown", "dbh"]
    with output_path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(headers[: points.shape[1]])
        writer.writerows(points.tolist())


def _show_warning(parent: QWidget, text: str) -> None:
    msg = QMessageBox(parent)
    msg.setIcon(QMessageBox.Icon.Warning)
    msg.setText(text)
    msg.setWindowTitle("Export unavailable")
    msg.setStandardButtons(QMessageBox.StandardButton.Ok)
    msg.exec()
