"""
    Copyright (C) 2024  University of Washington
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.

Bryce Bartl-Geller
University of Washington Forest Resilience Lab
12/6/2024

projectsettingsactivity.py
Handles user input and hands it off to ``ProjectSettings`` for processing.
"""

from __future__ import annotations

from pathlib import Path
import traceback
from typing import Any

from PySide6.QtCore import QTimer, Slot
from PySide6.QtWidgets import (
    QCheckBox,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QTextBrowser,
    QVBoxLayout,
)

from activity import Activity, SavedState, WindowMode
from gamingactivity import GamingActivity
from persistence import build_form_state, write_project_settings_file
from rxgaming_core import (
    ProjectArea,
    ProjectSettings,
    finish_project_area_build,
    poll_project_area_build,
    start_project_area_build,
)
from widgets import QFileSelectionLineEdit

PROGRESS_POLL_INTERVAL_MS = 100


class ProjectSettingsActivity(Activity):
    """Collect project inputs and build the native project area in the background."""

    def on_start(self, saved_state: SavedState, **kwargs: Any) -> None:
        self.save_file_location = ""
        self.prop_table_path = kwargs["prop_table_path"]
        self.fia_path = kwargs["fia_path"]
        self.worker = None
        self._closed = False
        self._last_progress_line = ""
        self._progress_timer = QTimer(self.window)
        self._progress_timer.setInterval(PROGRESS_POLL_INTERVAL_MS)
        self._progress_timer.timeout.connect(self._poll_worker_progress)

        assert isinstance(self.prop_table_path, str)
        assert isinstance(self.fia_path, str)

        self.prj_name_edit = QLineEdit()
        self.unit_poly_path_edit = QFileSelectionLineEdit(filter="ESRI Shapefile (*.shp)")
        self.reference_data_path_edit = QFileSelectionLineEdit(filter="Comma-separated values (*.csv)")
        self.lidar_data_path_edit = QFileSelectionLineEdit(file_type=QFileSelectionLineEdit.FileType.Directory)
        self.unit_name_edit = QLineEdit()
        self.threads_edit = QSpinBox()
        self.threads_edit.setMinimum(1)
        self.auto_save_checkbox = QCheckBox()
        self.auto_save_line_edit = QFileSelectionLineEdit(
            caption="Choose project folder...",
            file_type=QFileSelectionLineEdit.FileType.Directory,
        )
        self.auto_save_line_edit.setEnabled(False)

        form_layout = QFormLayout()
        form_layout.addRow("Project name", self.prj_name_edit)
        form_layout.addRow("Treatment unit polygon shapefile", self.unit_poly_path_edit)
        form_layout.addRow("Lidar data root folder", self.lidar_data_path_edit)
        form_layout.addRow("Reference data table (optional)", self.reference_data_path_edit)
        form_layout.addRow("Unit name field (optional)", self.unit_name_edit)
        form_layout.addRow("Number of threads to process on:", self.threads_edit)
        form_layout.addRow("Auto-Save?", self.auto_save_checkbox)
        form_layout.addRow("Project Folder", self.auto_save_line_edit)

        self.start_button = QPushButton("Start")
        self.save_button = QPushButton("Save settings")
        self.save_as_button = QPushButton("Save settings as")

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.start_button)
        button_layout.addWidget(self.save_button)
        button_layout.addWidget(self.save_as_button)

        outer_layout = QVBoxLayout()
        outer_layout.addLayout(form_layout)
        outer_layout.addLayout(button_layout)

        self.text_output = QTextBrowser()
        outer_layout.addWidget(self.text_output)

        self.window.setLayout(outer_layout)
        self.window.setWindowTitle("Project settings. Rxgaming tool version: {0}".format(Activity.version))

        self.auto_save_checkbox.stateChanged.connect(self.save_changed)
        self.start_button.clicked.connect(self.start_clicked)
        self.save_button.clicked.connect(self.save_clicked)
        self.save_as_button.clicked.connect(self.save_as_clicked)

        form_state = saved_state.get("ProjectSettingsForm")
        if isinstance(form_state, dict):
            self._apply_form_state(form_state)
        value = saved_state.get("save_file_location")
        if value is not None:
            self.save_file_location = value

    def save(self) -> SavedState:
        return {
            "ProjectSettingsForm": self._collect_form_state(),
            "save_file_location": self.save_file_location,
        }

    def start_clicked(self, checked: bool = False) -> None:
        del checked

        auto_save_text = self.auto_save_line_edit.text().strip()
        auto_save_path = Path(auto_save_text) if self.auto_save_checkbox.isChecked() else None
        unit_poly_path = Path(self.unit_poly_path_edit.text())
        lidar_data_path = Path(self.lidar_data_path_edit.text())
        reference_text = self.reference_data_path_edit.text()
        reference_data_path = Path(reference_text) if reference_text else None

        if auto_save_path is not None:
            if auto_save_text == "":
                self.notify_exception("Enter a project folder path before enabling Auto-Save.")
                return
            if not auto_save_path.parent.exists():
                self.notify_exception("The parent folder for the project folder does not exist. Enter a valid folder path before continuing.")
                return

        if not unit_poly_path.exists():
            self.notify_exception("The unit polygon file path does not exist. Enter a valid file path before continuing.")
            return
        if unit_poly_path.suffix.lower() != ".shp":
            self.notify_exception("The unit polygon expected file type is shapefile. Enter a valid file before continuing.")
            return

        if not lidar_data_path.exists():
            self.notify_exception("The lidar dataset path provided does not exist. Enter a valid file path before continuing.")
            return
        if not lidar_data_path.is_dir():
            self.notify_exception("The lidar dataset is expected to be a folder. Enter a valid folder before continuing.")
            return

        if reference_data_path is not None and reference_data_path.suffix.lower() != ".csv":
            self.notify_exception("The reference dataset expected file type is csv. Enter a valid csv before continuing.")
            return

        self.text_output.clear()
        self._append_progress("Processing started:")
        self.start_button.setEnabled(False)
        self._last_progress_line = ""

        self.worker = ProjectBuildWorker(
            self.prj_name_edit.text(),
            str(unit_poly_path),
            str(reference_data_path) if reference_data_path is not None else "",
            self.prop_table_path,
            self.fia_path,
            str(lidar_data_path),
            self.unit_name_edit.text(),
            self.auto_save_line_edit.text() if self.auto_save_checkbox.isChecked() else "",
            self.threads_edit.value(),
        )
        try:
            self.worker.start()
        except Exception as exc:
            self.start_button.setEnabled(True)
            message = "{0}\n\nDebugging Traceback:\n{1}".format(exc, traceback.format_exc())
            self.notify_exception(message)
            self.worker = None
            return
        self._progress_timer.start()

    def start_gamingactivity(self, result: Any) -> None:
        if not isinstance(result, dict):
            self.notify_exception(str(result))
            return

        project_settings = result.get("project_settings")
        project_area = result.get("project_area")
        if not isinstance(project_settings, ProjectSettings) or not isinstance(project_area, ProjectArea):
            self.notify_exception("Project build did not return a valid ProjectSettings and ProjectArea.")
            return

        try:
            autosave_path = self.auto_save_line_edit.text() if self.auto_save_checkbox.isChecked() else None
            Activity.start_activity(
                GamingActivity,
                None,
                {
                    "ProjectSettings": project_settings,
                    "ProjectArea": project_area,
                    "ProjectSettingsForm": self._collect_form_state(),
                    "ProjectSnapshotPath": autosave_path,
                    "SessionState": {},
                },
                WindowMode.Sibling,
            )
            self.stop()
        except Exception as exc:
            message = "{0}\n\nDebugging Traceback:\n{1}".format(exc, traceback.format_exc())
            self.notify_exception(message)

    def save_changed(self, state: int) -> None:
        self.auto_save_line_edit.setEnabled(bool(state))

    def save_clicked(self, checked: bool = False) -> None:
        del checked
        if self.save_file_location:
            save_path = Path(self.save_file_location)
            if save_path.parent.exists():
                write_project_settings_file(save_path, self._collect_form_state(), app_version=Activity.version)
                self.notify_save_success("Settings file saved successfully.")
                return
        self.save_as_clicked()

    def save_as_clicked(self, checked: bool = False) -> None:
        del checked
        from PySide6.QtWidgets import QFileDialog

        selected_path = QFileDialog.getSaveFileName(
            self.window,
            "Save settings file as...",
            str(self._default_settings_path()),
            "RxGaming settings files (*.json)",
        )[0]
        if selected_path == "":
            return

        save_path = Path(selected_path)
        write_project_settings_file(save_path, self._collect_form_state(), app_version=Activity.version)
        self.save_file_location = str(save_path)
        self.notify_save_success("Settings file saved successfully.")

    @Slot(str)
    def _append_progress(self, text: str) -> None:
        if not text or self._closed:
            return
        if self.text_output.toPlainText():
            self.text_output.insertPlainText("\n")
        self.text_output.insertPlainText(text)
        self.text_output.ensureCursorVisible()

    @Slot()
    def _poll_worker_progress(self) -> None:
        if self._closed or self.worker is None:
            return
        try:
            snapshot = self.worker.poll_snapshot()
        except Exception as exc:
            self._stop_progress_timer()
            self.start_button.setEnabled(True)
            self._dispose_worker()
            self.notify_exception("{0}\n\nDebugging Traceback:\n{1}".format(exc, traceback.format_exc()))
            return

        rendered_line = self._render_progress_snapshot(snapshot)
        if rendered_line and rendered_line != self._last_progress_line:
            self._append_progress(rendered_line)
            self._last_progress_line = rendered_line

        status = str(getattr(snapshot, "status", ""))
        if status not in {"succeeded", "failed"}:
            return

        self._stop_progress_timer()
        worker = self.worker
        assert worker is not None
        try:
            project_area = worker.finalize()
            project_settings = worker.project_settings
            if project_settings is None or not isinstance(project_area, ProjectArea):
                raise TypeError("Project build did not return a valid ProjectSettings and ProjectArea.")
        except Exception as exc:
            self.start_button.setEnabled(True)
            self._dispose_worker()
            self.notify_exception("{0}\n\nDebugging Traceback:\n{1}".format(exc, traceback.format_exc()))
            return

        self.start_button.setEnabled(True)
        self._append_progress("Project area ready.")
        result = {
            "project_settings": project_settings,
            "project_area": project_area,
        }
        self._dispose_worker()
        self.start_gamingactivity(result)

    def _render_progress_snapshot(self, snapshot: Any) -> str:
        message = str(getattr(snapshot, "message", "") or getattr(snapshot, "error", ""))
        completed = getattr(snapshot, "completed", -1)
        total = getattr(snapshot, "total", -1)

        if isinstance(completed, int) and isinstance(total, int) and total > 0 and completed >= 0:
            suffix = f"[{completed}/{total}]"
            if message:
                if suffix not in message:
                    return f"{message} {suffix}"
                return message
            return f"Progress {suffix}"

        return message

    def _stop_progress_timer(self) -> None:
        if self._progress_timer.isActive():
            self._progress_timer.stop()

    def _dispose_worker(self) -> None:
        worker = self.worker
        self.worker = None
        if worker is not None:
            worker.dispose()

    def _on_window_close(self, window: Any) -> None:
        self._closed = True
        self._stop_progress_timer()
        self._dispose_worker()
        super()._on_window_close(window)

    @staticmethod
    def notify_save_success(text: str) -> None:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Information)
        msg.setText(text)
        msg.setWindowTitle("Result")
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()

    @staticmethod
    def notify_exception(text: str) -> None:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Warning)
        msg.setText(
            "There were errors with some of the parameters provided.\n"
            "Often a description of the error is provided below, "
            "with further in depth traceback reporting after. "
            "If you would like to report this error or believe you have encountered a bug, "
            "Please email this text to bnbg@uw.edu. Thank you:"
        )
        msg.setInformativeText(text)
        msg.setWindowTitle("Parameter error")
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()

    def _collect_form_state(self) -> dict[str, Any]:
        return build_form_state(
            project_name=self.prj_name_edit.text(),
            unit_poly_path=self.unit_poly_path_edit.text(),
            reference_data_path=self.reference_data_path_edit.text(),
            lidar_data_path=self.lidar_data_path_edit.text(),
            unit_name=self.unit_name_edit.text(),
            threads=self.threads_edit.value(),
            auto_save_enabled=self.auto_save_checkbox.isChecked(),
            auto_save_path=self.auto_save_line_edit.text(),
        )

    def _apply_form_state(self, form_state: dict[str, Any]) -> None:
        self.prj_name_edit.setText(str(form_state.get("project_name", "")))
        self.unit_poly_path_edit.set_text(str(form_state.get("unit_poly_path", "")))
        self.reference_data_path_edit.set_text(str(form_state.get("reference_data_path", "")))
        self.lidar_data_path_edit.set_text(str(form_state.get("lidar_data_path", "")))
        self.unit_name_edit.setText(str(form_state.get("unit_name", "")))
        self.threads_edit.setValue(int(form_state.get("threads", 1)))
        self.auto_save_checkbox.setChecked(bool(form_state.get("auto_save_enabled", False)))
        self.auto_save_line_edit.set_text(str(form_state.get("auto_save_path", "")))

    def _default_settings_path(self) -> Path:
        if self.save_file_location:
            return Path(self.save_file_location)
        return Path.cwd() / "settings.json"


class ProjectBuildWorker:
    def __init__(
        self,
        prj_name: str,
        unit: str,
        ref: str,
        prop_table: str,
        fia: str,
        lidar: str,
        unit_name: str,
        save_path: str,
        threads: int,
    ):
        self.prj_name = prj_name
        self.unit = unit
        self.ref = ref
        self.prop_table = prop_table
        self.fia = fia
        self.lidar = lidar
        self.unit_name = unit_name
        self.save_path = save_path
        self.threads = threads
        self.handle = None
        self.project_settings: ProjectSettings | None = None

    def start(self) -> None:
        self.project_settings = ProjectSettings(
            self.prj_name,
            self.unit,
            self.ref,
            self.prop_table,
            self.fia,
            self.lidar,
            self.unit_name,
            self.save_path,
            self.threads,
        )
        self.handle = start_project_area_build(self.project_settings)

    def poll_snapshot(self) -> Any:
        if self.handle is None:
            raise RuntimeError("No project build is currently running.")
        return poll_project_area_build(self.handle)

    def finalize(self) -> ProjectArea:
        if self.handle is None:
            raise RuntimeError("No project build is currently running.")
        return finish_project_area_build(self.handle)

    def dispose(self) -> None:
        self.handle = None
