"""
Gaming activity UI shell.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from PySide6.QtGui import QAction
from PySide6.QtWidgets import QMessageBox, QVBoxLayout, QWidget

from activity import Activity
from gaming_export import (
    export_basins_raster,
    export_chm_raster,
    export_clumpmap_raster,
    export_current_view,
    export_features,
    export_georeferenced_raster,
    export_treelist,
)
from gaming_ui import GamingTabs
from persistence import ProjectSnapshotSessionPersistence
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
        self.project_settings_form = saved_state.get("ProjectSettingsForm")
        self.project_snapshot_path = saved_state.get("ProjectSnapshotPath")
        self.project_settings = self._load_project_settings(saved_state)
        self.project_area = self._load_project_area(saved_state, self.project_settings)
        self._persistence = self._build_persistence()
        self.tab_widget = GamingTabs(
            self.project_settings,
            self.project_area,
            saved_state,
            persistence=self._persistence,
        )

        layout = QVBoxLayout()
        layout.addWidget(self.tab_widget)
        self.central_widget.setLayout(layout)

        self.window.setWindowTitle(f"{self.project_settings.name} Gaming. Rxgaming tool version: {Activity.version}")
        self.window.setGeometry(200, 200, 1800, 1000)
        self.window.setMinimumSize(1200, 800)
        self.window.showMaximized()

        self._create_actions()
        self._create_menu()

        if self._persistence is not None and self.project_snapshot_path is not None:
            self._persistence.initialize_snapshot(self.tab_widget.session_state)

    def save(self) -> dict[str, Any]:
        return {
            "ProjectSettings": self.project_settings,
            "ProjectArea": self.project_area,
            "ProjectSettingsForm": self.project_settings_form,
            "ProjectSnapshotPath": self.project_snapshot_path,
            "SessionState": self.tab_widget.session_state.to_dict(),
            "LastActivity": type(self),
        }

    def menu_exit(self) -> None:
        self.stop()

    def save_project(self) -> None:
        if self._persistence is None:
            self.save_project_as()
            return
        try:
            self._persistence.save_full_project(self.tab_widget.session_state)
            self._notify_save_success("Project saved successfully.")
        except Exception as exc:
            self._notify_save_failure(str(exc))

    def save_project_as(self) -> None:
        from PySide6.QtWidgets import QFileDialog

        selected_path = QFileDialog.getExistingDirectory(
            self.window,
            "Save project as...",
            str(self._default_project_root()),
        )
        if selected_path == "":
            return

        selected_root = Path(selected_path)
        if self._project_root_has_existing_save(selected_root) and not self._confirm_overwrite_project(selected_root):
            return
        self.project_snapshot_path = str(selected_root)
        self._persistence = self._build_persistence()
        if self._persistence is None:
            self._notify_save_failure("Could not create a project persistence target.")
            return

        try:
            self._persistence.save_full_project(self.tab_widget.session_state)
            self._notify_save_success("Project saved successfully.")
        except Exception as exc:
            self._notify_save_failure(str(exc))

    def export_tif(self) -> None:
        export_current_view(self.tab_widget, self.window)

    def export_chm_raster(self) -> None:
        export_chm_raster(self.tab_widget, self.window)

    def export_basins_raster(self) -> None:
        export_basins_raster(self.tab_widget, self.window)

    def export_clumpmap_raster(self) -> None:
        export_clumpmap_raster(self.tab_widget, self.window)

    def export_georeferenced_raster(self) -> None:
        export_georeferenced_raster(self.tab_widget, self.window)

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
        self.save_action = QAction("&Save Project", self.window)
        self.save_action.setShortcut("Ctrl+S")
        self.save_action.triggered.connect(self.save_project)

        self.save_as_action = QAction("Save Project &As", self.window)
        self.save_as_action.triggered.connect(self.save_project_as)

        self.exit_action = QAction("&Exit", self.window)
        self.exit_action.setShortcut("Ctrl+Q")
        self.exit_action.triggered.connect(self.menu_exit)

        self.export_tif_action = QAction('&Export window image ("*.png")', self.window)
        self.export_tif_action.triggered.connect(self.export_tif)

        self.export_chm_raster_action = QAction('&Export CHM ("*.tif")', self.window)
        self.export_chm_raster_action.triggered.connect(self.export_chm_raster)

        self.export_basins_raster_action = QAction('&Export basins ("*.tif")', self.window)
        self.export_basins_raster_action.triggered.connect(self.export_basins_raster)

        self.export_clumpmap_raster_action = QAction('&Export clumpmap ("*.tif")', self.window)
        self.export_clumpmap_raster_action.triggered.connect(self.export_clumpmap_raster)

        self.export_georeferenced_raster_action = QAction('&Export georeferenced raster image ("*.tif")', self.window)
        self.export_georeferenced_raster_action.triggered.connect(self.export_georeferenced_raster)

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
        file_menu.addAction(self.save_action)
        file_menu.addAction(self.save_as_action)
        file_menu.addAction(self.exit_action)

        export_menu.addAction(self.export_chm_raster_action)
        export_menu.addAction(self.export_basins_raster_action)
        export_menu.addAction(self.export_clumpmap_raster_action)
        export_menu.addAction(self.export_georeferenced_raster_action)
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
        project_area = saved_state.get("ProjectArea")
        if isinstance(project_area, ProjectArea):
            return project_area
        return ProjectArea(project_settings)

    def _build_persistence(self) -> ProjectSnapshotSessionPersistence | None:
        if not self.project_snapshot_path:
            return None
        return ProjectSnapshotSessionPersistence(
            self.project_snapshot_path,
            app_version=Activity.version,
            project_settings=self.project_settings,
            project_area=self.project_area,
            form_state=self.project_settings_form if isinstance(self.project_settings_form, dict) else None,
        )

    def _default_project_root(self) -> Path:
        if self.project_snapshot_path:
            return Path(self.project_snapshot_path)
        save_path = getattr(self.project_settings, "savePath", "")
        if save_path:
            return Path(save_path)
        return Path.cwd()

    @staticmethod
    def _project_root_has_existing_save(project_root: Path) -> bool:
        canonical_files = (
            project_root / "project.json",
            project_root / "settings.json",
            project_root / "session.json",
            project_root / "projectarea.h5",
        )
        return any(path.exists() for path in canonical_files)

    def _confirm_overwrite_project(self, project_root: Path) -> bool:
        result = QMessageBox.question(
            self.window,
            "Overwrite existing project?",
            (
                "This folder already contains a saved project.\n\n"
                f"Saving here will overwrite project files in:\n{project_root}"
            ),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        return result == QMessageBox.StandardButton.Yes

    @staticmethod
    def _notify_save_success(text: str) -> None:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Information)
        msg.setText(text)
        msg.setWindowTitle("Save successful")
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()

    @staticmethod
    def _notify_save_failure(text: str) -> None:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Warning)
        msg.setText("Saving the project failed.")
        msg.setInformativeText(text)
        msg.setWindowTitle("Save failed")
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()
