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
from typing import Any, Optional

from PySide6.QtCore import QObject, QThread, Signal, Slot
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
from rxgaming_core import ProjectArea, ProjectSettings
from widgets import QFileSelectionLineEdit

PROGRESS_UI_UPDATE_INTERVAL = 25

try:
    from rxgaming_core import build_project_area_with_progress
except ImportError:
    def build_project_area_with_progress(ps: ProjectSettings, callback: Any = None) -> ProjectArea:
        del callback
        return ProjectArea(ps)


class ProjectSettingsActivity(Activity):
    """Collect project inputs and build the native project area in the background."""

    def on_start(self, saved_state: SavedState, **kwargs: Any) -> None:
        self.save_file_location = ""
        self.prop_table_path = kwargs["prop_table_path"]
        self.fia_path = kwargs["fia_path"]
        self.worker_thread = None
        self.worker = None
        self._closed = False
        self._pending_worker_result = None
        self._worker_finished = False

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
        self._pending_worker_result = None
        self._worker_finished = False

        self.worker_thread = QThread()
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
        self.worker.moveToThread(self.worker_thread)
        self.worker_thread.started.connect(self.worker.run)
        self.worker.progress_message.connect(self._append_progress)
        self.worker.finished.connect(self._handle_worker_result)
        self.worker.finished.connect(self.worker_thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker_thread.finished.connect(self.worker_thread.deleteLater)
        self.worker_thread.finished.connect(self._handle_worker_thread_finished)
        self.worker_thread.start()

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

    @Slot(object)
    def _handle_worker_result(self, result: Any) -> None:
        self._pending_worker_result = result
        self._worker_finished = True
        if self.worker_thread is None or not self.worker_thread.isRunning():
            self._finalize_worker_result()

    @Slot()
    def _handle_worker_thread_finished(self) -> None:
        self._clear_worker_references()
        if self._worker_finished:
            self._finalize_worker_result()

    def _finalize_worker_result(self) -> None:
        self.start_button.setEnabled(True)
        if self._closed:
            return

        result = self._pending_worker_result
        self._pending_worker_result = None
        self._worker_finished = False

        if isinstance(result, dict):
            self._append_progress("Project area ready.")
            self.start_gamingactivity(result)
            return

        self.notify_exception(str(result))

    @Slot()
    def _clear_worker_references(self) -> None:
        self.worker_thread = None
        self.worker = None

    def _shutdown_worker(self) -> None:
        thread = self.worker_thread
        if thread is None:
            return
        if thread.isRunning():
            thread.quit()
            thread.wait(5000)
        self._clear_worker_references()

    def _on_window_close(self, window: Any) -> None:
        self._closed = True
        self._shutdown_worker()
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


class ProjectBuildWorker(QObject):
    progress_message = Signal(str)
    finished = Signal(object)

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
        *args: Any,
        **kwargs: Any
    ):
        QObject.__init__(self, *args, **kwargs)
        self.prj_name = prj_name
        self.unit = unit
        self.ref = ref
        self.prop_table = prop_table
        self.fia = fia
        self.lidar = lidar
        self.unit_name = unit_name
        self.save_path = save_path
        self.threads = threads

    @staticmethod
    def _format_progress_update(event: Any) -> str:
        stage = str(getattr(event, "stage", ""))
        message = str(getattr(event, "message", ""))
        completed = getattr(event, "completed", -1)
        total = getattr(event, "total", -1)

        if stage == "unit_start":
            return ""

        if stage in {"unit_complete", "unit_skipped"} and isinstance(completed, int) and isinstance(total, int):
            if completed <= 0 or total <= 0:
                return ""
            if completed != total and completed % PROGRESS_UI_UPDATE_INTERVAL != 0:
                return ""
            return f"{message} [{completed}/{total}]"

        return message

    @Slot()
    def run(self) -> None:
        try:
            project_settings = ProjectSettings(
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

            def on_progress(event: Any) -> None:
                message = self._format_progress_update(event)
                if message:
                    self.progress_message.emit(message)

            project_area = build_project_area_with_progress(project_settings, None)
            self.finished.emit(
                {
                    "project_settings": project_settings,
                    "project_area": project_area,
                }
            )
        except Exception as exc:
            message = "{0}\n\nDebugging Traceback:\n{1}".format(exc, traceback.format_exc())
            self.finished.emit(message)
