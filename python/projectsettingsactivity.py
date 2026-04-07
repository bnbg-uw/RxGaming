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
import pickle
from queue import Empty, Queue
import sys
import traceback

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

from activity import Activity, SavedState, SaveStateActivity, WindowMode
from gamingactivity import GamingActivity
from rxgaming_core import ProjectSettings
from widgets import QFileSelectionLineEdit


class ProjectSettingsActivity(Activity):
    """Collect project inputs and build a ``ProjectSettings`` instance."""

    def on_start(self, saved_state: SavedState, prop_table_path: str) -> None:
        self.save_file_location = ""
        self.prop_table_path = prop_table_path

        self.prj_name_edit = QLineEdit()
        self.unit_poly_path_edit = QFileSelectionLineEdit(filter="ESRI Shapefile (*.shp)")
        self.reference_data_path_edit = QFileSelectionLineEdit(filter="Comma-separated values (*.csv)")
        self.lidar_data_path_edit = QFileSelectionLineEdit(file_type=QFileSelectionLineEdit.FileType.Directory)
        self.unit_name_edit = QLineEdit()
        self.threads_edit = QSpinBox()
        self.threads_edit.setMinimum(1)
        self.auto_save_checkbox = QCheckBox()
        self.auto_save_line_edit = QFileSelectionLineEdit(filter="*.dat", new_file=True)
        self.auto_save_line_edit.setEnabled(False)
        self.fia_path = ""
        self.proj_db_path = ""

        form_layout = QFormLayout()
        form_layout.addRow("Project name", self.prj_name_edit)
        form_layout.addRow("Treatment unit polygon shapefile", self.unit_poly_path_edit)
        form_layout.addRow("Lidar data root folder", self.lidar_data_path_edit)
        form_layout.addRow("Reference data table (optional)", self.reference_data_path_edit)
        form_layout.addRow("Unit name field (optional)", self.unit_name_edit)
        form_layout.addRow("Number of threads to process on:", self.threads_edit)
        form_layout.addRow("Auto-Save?", self.auto_save_checkbox)
        form_layout.addRow("Save File Location", self.auto_save_line_edit)

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
        self.window.setWindowTitle(f"Project settings. Rxgaming tool version: {Activity.version}")

        self.auto_save_checkbox.stateChanged.connect(self.save_changed)
        self.start_button.clicked.connect(self.start_clicked)
        self.save_button.clicked.connect(self.save_clicked)
        self.save_as_button.clicked.connect(self.save_as_clicked)

        if (value := saved_state.get("ProjectSettingsActivity.prj_name")) is not None:
            self.prj_name_edit.setText(value)
        if (value := saved_state.get("ProjectSettingsActivity.unit_poly_path")) is not None:
            self.unit_poly_path_edit.set_text(value)
        if (value := saved_state.get("ProjectSettingsActivity.reference_data_path")) is not None:
            self.reference_data_path_edit.set_text(value)
        if (value := saved_state.get("ProjectSettingsActivity.lidar_data_path")) is not None:
            self.lidar_data_path_edit.set_text(value)
        if (value := saved_state.get("ProjectSettingsActivity.unit_name")) is not None:
            self.unit_name_edit.setText(value)
        if (value := saved_state.get("ProjectSettingsActivity.threads")) is not None:
            self.threads_edit.setValue(value)
        if (value := saved_state.get("ProjectSettingsActivity.auto_save")) is not None:
            self.auto_save_checkbox.setChecked(value)
        if (value := saved_state.get("ProjectSettingsActivity.auto_save_line")) is not None:
            self.auto_save_line_edit.set_text(value)
        if (value := saved_state.get("save_file_location")) is not None:
            self.save_file_location = value

    def save(self) -> SavedState:
        return {
            "ProjectSettingsActivity.prj_name": self.prj_name_edit.text(),
            "ProjectSettingsActivity.unit_poly_path": self.unit_poly_path_edit.text(),
            "ProjectSettingsActivity.reference_data_path": self.reference_data_path_edit.text(),
            "ProjectSettingsActivity.lidar_data_path": self.lidar_data_path_edit.text(),
            "ProjectSettingsActivity.unit_name": self.unit_name_edit.text(),
            "ProjectSettingsActivity.threads": self.threads_edit.value(),
            "ProjectSettingsActivity.auto_save": self.auto_save_checkbox.isChecked(),
            "ProjectSettingsActivity.auto_save_line": self.auto_save_line_edit.text(),
        }

    def start_clicked(self, checked: bool = False) -> None:
        auto_save_path = Path(self.auto_save_line_edit.text()) if self.auto_save_checkbox.isChecked() else None
        unit_poly_path = Path(self.unit_poly_path_edit.text())
        lidar_data_path = Path(self.lidar_data_path_edit.text())
        reference_text = self.reference_data_path_edit.text()
        reference_data_path = Path(reference_text) if reference_text else None

        if auto_save_path is not None and not auto_save_path.parent.exists():
            self.notify_exception("The Auto-Save file path does not exist. Enter a valid file path before continuing.")
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

        output_queue: Queue[str] = Queue()
        sys.stdout = WriteStream(output_queue)
        print("Processing started:")

        self.stream_thread = QThread()
        self.output_receiver = Receiver(output_queue)
        self.output_receiver.text_stream.connect(self.text_output.append)
        self.output_receiver.moveToThread(self.stream_thread)
        self.stream_thread.started.connect(self.output_receiver.run)
        self.stream_thread.start()

        self.worker_thread = QThread()
        self.worker = Worker(
            project_name=self.prj_name_edit.text(),
            unit_poly_path=str(unit_poly_path),
            reference_data_path=str(reference_data_path) if reference_data_path is not None else "",
            lidar_data_path=str(lidar_data_path),
            unit_name=self.unit_name_edit.text(),
            threads=self.threads_edit.value(),
            prop_table_path=self.prop_table_path,
            fia_path=self.fia_path,
            proj_db_path=self.proj_db_path,
            save_path=self.auto_save_line_edit.text() if self.auto_save_checkbox.isChecked() else "",
        )
        self.worker.moveToThread(self.worker_thread)
        self.worker_thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.worker_thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker_thread.finished.connect(self.worker_thread.deleteLater)
        self.worker.finished.connect(self._cleanup_stream_thread)
        self.worker.finished.connect(self.start_gaming_activity)
        self.worker.finished.connect(self._restore_stdout)

        self.start_button.setEnabled(False)
        self.worker.finished.connect(lambda _result: self.start_button.setEnabled(True))
        self.worker_thread.start()

    def start_gaming_activity(self, project_settings: ProjectSettings | str) -> None:
        if not isinstance(project_settings, ProjectSettings):
            self.notify_exception(project_settings)
            return

        try:
            autosave_path = self.auto_save_line_edit.text() if self.auto_save_checkbox.isChecked() else None
            Activity.start_activity(
                GamingActivity,
                None,
                {"ProjectSettings": project_settings, "Autosave_path": autosave_path},
                WindowMode.Sibling,
            )
            self.stop()
        except Exception as exc:
            message = f"{exc}\n\nDebugging Traceback:\n{traceback.format_exc()}"
            self.notify_exception(message)

    def save_changed(self, state: int) -> None:
        self.auto_save_line_edit.setEnabled(bool(state))

    def save_clicked(self, checked: bool = False) -> None:
        if self.save_file_location:
            save_path = Path(self.save_file_location)
            if save_path.is_file():
                with save_path.open("wb") as fp:
                    pickle.dump(self.save(), fp)
                self.notify_save_success()
                return
        self.save_as_clicked()

    def save_as_clicked(self, checked: bool = False) -> None:
        SaveStateActivity.prompt_and_save(SaveStateActivity, self.save())
        self.notify_save_success()

    def _cleanup_stream_thread(self, _result: ProjectSettings | str) -> None:
        if self.stream_thread.isRunning():
            self.output_receiver.stop()
            self.stream_thread.quit()
            self.output_receiver.deleteLater()
            self.stream_thread.deleteLater()

    @staticmethod
    def _restore_stdout(_result: ProjectSettings | str) -> None:
        sys.stdout = sys.__stdout__
        print("stdout restored")

    @staticmethod
    def notify_save_success() -> None:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Information)
        msg.setText("Save Successful!")
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


class WriteStream:
    def __init__(self, queue: Queue[str]) -> None:
        self.queue = queue

    def write(self, text: str) -> None:
        self.queue.put(text)

    def flush(self) -> None:
        pass


class Receiver(QObject):
    text_stream = Signal(str)

    def __init__(self, queue: Queue[str]) -> None:
        super().__init__()
        self.queue = queue
        self.running = False

    def stop(self) -> None:
        self.running = False

    @Slot()
    def run(self) -> None:
        self.running = True
        while self.running:
            try:
                text = self.queue.get_nowait()
            except Empty:
                QThread.msleep(50)
                continue

            if text:
                sys.__stdout__.write(text)
                self.text_stream.emit(text)


class Worker(QObject):
    finished = Signal(object)

    def __init__(
        self,
        project_name: str,
        unit_poly_path: str,
        reference_data_path: str,
        lidar_data_path: str,
        unit_name: str,
        threads: int,
        prop_table_path: str,
        fia_path: str,
        proj_db_path: str,
        save_path: str,
    ) -> None:
        super().__init__()
        self.project_name = project_name
        self.unit_poly_path = unit_poly_path
        self.reference_data_path = reference_data_path
        self.lidar_data_path = lidar_data_path
        self.unit_name = unit_name
        self.threads = threads
        self.prop_table_path = prop_table_path
        self.fia_path = fia_path
        self.proj_db_path = proj_db_path
        self.save_path = save_path

    @Slot()
    def run(self) -> None:
        try:
            project_settings = ProjectSettings(
                self.project_name,
                self.unit_poly_path,
                self.reference_data_path,
                self.prop_table_path,
                self.fia_path,
                self.proj_db_path,
                self.lidar_data_path,
                self.unit_name,
                self.save_path,
                self.threads,
            )
            self.finished.emit(project_settings)
        except Exception as exc:
            message = f"{exc}\n\nDebugging Traceback:\n{traceback.format_exc()}"
            self.finished.emit(message)

