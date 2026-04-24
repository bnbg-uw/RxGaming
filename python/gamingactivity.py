"""
Gaming activity UI shell.
"""

from __future__ import annotations

from typing import Any

from PySide6.QtGui import QAction
from PySide6.QtWidgets import QMessageBox, QVBoxLayout, QWidget

from activity import Activity
from gaming_export import export_current_view, export_features, export_raster, export_treelist
from gaming_ui import GamingTabs
from rxgaming_core import ProjectArea, ProjectSettings
from widgets import QMainWindowRx


class GamingActivity(Activity):
    """Main activity for viewing units and running treatments."""

    def on_start(self, saved_state: dict[str, Any], **kwargs: Any) -> None:
        del kwargs
        self.set_window(QMainWindowRx())
        self.central_widget = QWidget()
        self.window.setCentralWidget(self.central_widget)

        self.saved_state = saved_state
        print("before load ps")
        self.project_settings = self._load_project_settings(saved_state)
        print("before load pa")
        self.project_area = self._load_project_area(saved_state, self.project_settings)
        self.tab_widget = GamingTabs(self.project_settings, self.project_area, saved_state)

        layout = QVBoxLayout()
        layout.addWidget(self.tab_widget)
        self.central_widget.setLayout(layout)

        self.window.setWindowTitle(f"{self.project_settings.name} Gaming. Rxgaming tool version: {Activity.version}")
        self.window.setGeometry(200, 200, 1800, 1000)
        self.window.setMinimumSize(1200, 800)
        self.window.showMaximized()

        self._create_actions()
        self._create_menu()

    def save(self) -> dict[str, Any]:
        return {
            "ProjectSettings": self._serialize_project_settings(self.project_settings),
            "LastActivity": type(self),
        }

    def menu_exit(self) -> None:
        self.stop()

    def export_tif(self) -> None:
        export_current_view(self.tab_widget, self.window)

    def export_rasters(self) -> None:
        export_raster(self.tab_widget, self.window)

    def export_features(self) -> None:
        export_features(self.tab_widget, self.window)

    def export_treelists(self) -> None:
        export_treelist(self.tab_widget, self.window)

    @staticmethod
    def license() -> None:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Information)
        msg.setText(
            "Copyright (C) 2024  University of Washington\n\n"
            "This program is free software: you can redistribute it and/or modify it under the terms "
            "of the GNU General Public License as published by the Free Software Foundation, either "
            "version 3 of the License, or (at your option) any later version.\n\n"
            "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; "
            "without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. "
            "See the GNU General Public License for more details.\n\n"
            "You should have received a copy of the GNU General Public License along with this program. "
            "If not, see https://www.gnu.org/licenses/."
        )
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()

    def _create_actions(self) -> None:
        self.exit_action = QAction("&Exit", self.window)
        self.exit_action.setShortcut("Ctrl+Q")
        self.exit_action.triggered.connect(self.menu_exit)

        self.export_tif_action = QAction('&Export window image ("*.png")', self.window)
        self.export_tif_action.triggered.connect(self.export_tif)

        self.export_rasters_action = QAction('&Export raster data ("*.npy")', self.window)
        self.export_rasters_action.triggered.connect(self.export_rasters)

        self.export_features_action = QAction('&Export point data ("*.csv")', self.window)
        self.export_features_action.triggered.connect(self.export_features)

        self.export_treelists_action = QAction('&Export treelists ("*.csv")', self.window)
        self.export_treelists_action.triggered.connect(self.export_treelists)

        self.license_action = QAction("&License", self.window)
        self.license_action.triggered.connect(self.license)

    def _create_menu(self) -> None:
        main_menu = self.window.menuBar()
        file_menu = main_menu.addMenu("File")
        export_menu = file_menu.addMenu("Export")

        file_menu.addAction(self.license_action)
        file_menu.addAction(self.exit_action)

        export_menu.addAction(self.export_rasters_action)
        export_menu.addAction(self.export_features_action)
        export_menu.addAction(self.export_tif_action)
        export_menu.addAction(self.export_treelists_action)

    @staticmethod
    def _serialize_project_settings(project_settings: ProjectSettings) -> dict[str, Any]:
        return {
            "name": project_settings.name,
            "unitPolyPath": project_settings.unitPolyPath,
            "refDataPath": project_settings.refDataPath,
            "mcsPropPath": project_settings.mcsPropPath,
            "fiaPath": project_settings.fiaPath,
            "lidarPath": project_settings.lidarPath,
            "unitName": project_settings.unitName,
            "savePath": project_settings.savePath,
            "nThread": project_settings.nThread,
        }

    @classmethod
    def _load_project_settings(cls, saved_state: dict[str, Any]) -> ProjectSettings:
        project_settings = saved_state["ProjectSettings"]
        if isinstance(project_settings, ProjectSettings):
            return project_settings
        if isinstance(project_settings, dict):
            return ProjectSettings(**project_settings)
        raise TypeError("Saved ProjectSettings must be an rxgaming_core.ProjectSettings or a serialized settings dict.")

    @staticmethod
    def _load_project_area(saved_state: dict[str, Any], project_settings: ProjectSettings) -> ProjectArea:
        print("in load project area")
        project_area = saved_state.get("ProjectArea")
        if isinstance(project_area, ProjectArea):
            print("found project area")
            return project_area
        print("constructing project area")
        return ProjectArea(project_settings)
