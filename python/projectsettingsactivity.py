"""
    Copyright (C) 2024  University of Washington
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.

Bryce Bartl-Geller
University of Washington Forest Resilience Lab
12/6/2024

projectsettingsactivity.py
Handles user input and hands it of the ProjectSettings in projectsettings.py for processing.
"""

import os.path
import pickle
import sys
import traceback
from queue import Queue
from PyQt5.QtWidgets import (QFormLayout, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QMessageBox, QCheckBox,
                             QSpinBox, QTextBrowser)
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot, QThread
from QtWidgets import QFileSelectionLineEdit
from activity import Activity, WindowMode, SaveStateActivity
from projectsettings import ProjectSettings
from gamingactivity import GamingActivity


# ProjectSettingsActivity class
#  Handles user input for getting required files, then kicks off GamingActivity
class ProjectSettingsActivity(Activity):
    def on_start(self, saved_state, prop_table_path, dll_path):
        self.save_file_location = ""

        # MCS Prop table path:
        self.prop_table_path = prop_table_path
        self.dll_path = dll_path

        # Initialize UI
        self.prj_name_edit = QLineEdit()
        #self.prj_poly_path_edit = QFileSelectionLineEdit(filter="ESRI Shapefile (*.shp)")
        self.unit_poly_path_edit = QFileSelectionLineEdit(filter="ESRI Shapefile (*.shp)")
        self.reference_data_path_edit = QFileSelectionLineEdit(filter="Comma-separated values (*.csv)")
        self.lidar_data_path_edit = QFileSelectionLineEdit(file_type=QFileSelectionLineEdit.FileType.Directory)
        self.unit_name_edit = QLineEdit()
        self.threads_edit = QSpinBox()
        self.threads_edit.setMinimum(1)
        self.auto_save = QCheckBox()
        self.auto_save_line = QFileSelectionLineEdit(filter="*.dat", new_file=True)
        self.auto_save_line.setEnabled(False)

        self.form_layout = QFormLayout()
        self.form_layout.addRow("Project name", self.prj_name_edit)
        #self.form_layout.addRow("Project area polygon shapefile", self.prj_poly_path_edit)
        self.form_layout.addRow("Treatment unit polygon shapefile", self.unit_poly_path_edit)
        self.form_layout.addRow("Lidar data root folder", self.lidar_data_path_edit)
        self.form_layout.addRow("Reference data table (optional)", self.reference_data_path_edit)
        self.form_layout.addRow("Unit name field (optional)", self.unit_name_edit)
        self.form_layout.addRow("Number of threads to process on:", self.threads_edit)
        self.form_layout.addRow("Auto-Save?", self.auto_save)
        self.form_layout.addRow("Save File Location", self.auto_save_line)

        self.start_button = QPushButton("Start")
        self.save_button = QPushButton("Save settings")
        self.save_as_button = QPushButton("Save settings as")

        self.button_layout = QHBoxLayout()
        self.button_layout.addWidget(self.start_button)
        self.button_layout.addWidget(self.save_button)
        self.button_layout.addWidget(self.save_as_button)


        self.outer_layout = QVBoxLayout()
        self.outer_layout.addLayout(self.form_layout)
        self.outer_layout.addLayout(self.button_layout)

        self.text_output = QTextBrowser()
        self.outer_layout.addWidget(self.text_output)

        self.window.setLayout(self.outer_layout)
        self.window.setWindowTitle("Project settings. Rxgaming tool version: " + Activity.version)
        
        # Setup callbacks
        self.auto_save.stateChanged.connect(self.save_changed)
        self.start_button.clicked.connect(self.start_clicked)
        self.save_button.clicked.connect(self.save_clicked)
        self.save_as_button.clicked.connect(self.save_as_clicked)

        # Load saved state
        if "ProjectSettingsActivity.prj_name" in saved_state:
            self.prj_name_edit.setText(saved_state["ProjectSettingsActivity.prj_name"])
        #if "ProjectSettingsActivity.prj_poly_path" in saved_state:
        #    self.prj_poly_path_edit.setText(saved_state["ProjectSettingsActivity.prj_poly_path"])
        if "ProjectSettingsActivity.unit_poly_path" in saved_state:
            self.unit_poly_path_edit.setText(saved_state["ProjectSettingsActivity.unit_poly_path"])
        if "ProjectSettingsActivity.reference_data_path" in saved_state:
            self.reference_data_path_edit.setText(saved_state["ProjectSettingsActivity.reference_data_path"])
        if "ProjectSettingsActivity.lidar_data_path" in saved_state:
            self.lidar_data_path_edit.setText(saved_state["ProjectSettingsActivity.lidar_data_path"])
        if "ProjectSettingsActivity.unit_name" in saved_state:
            self.unit_name_edit.setText(saved_state['ProjectSettingsActivity.unit_name'])
        if "ProjectSettingsActivity.threads" in saved_state:
            self.threads_edit.setValue(saved_state['ProjectSettingsActivity.threads'])
        if "ProjectSettingsActivity.auto_save" in saved_state:
            self.auto_save.setChecked(saved_state['ProjectSettingsActivity.auto_save'])
        if "ProjectSettingsActivity.auto_save_line" in saved_state:
            self.auto_save_line.setText(saved_state['ProjectSettingsActivity.auto_save_line'])
        if "save_file_location" in saved_state:
            self.save_file_location = saved_state['save_file_location']

    def save(self):
        saved_state = dict()
        saved_state["ProjectSettingsActivity.prj_name"] = self.prj_name_edit.text()
        #saved_state["ProjectSettingsActivity.prj_poly_path"] = self.prj_poly_path_edit.text()
        saved_state["ProjectSettingsActivity.unit_poly_path"] = self.unit_poly_path_edit.text()
        saved_state["ProjectSettingsActivity.reference_data_path"] = self.reference_data_path_edit.text()
        saved_state["ProjectSettingsActivity.lidar_data_path"] = self.lidar_data_path_edit.text()
        saved_state["ProjectSettingsActivity.unit_name"] = self.unit_name_edit.text()
        saved_state["ProjectSettingsActivity.threads"] = self.threads_edit.value()
        saved_state["ProjectSettingsActivity.auto_save"] = self.auto_save.isChecked()
        saved_state['ProjectSettingsActivity.auto_save_line'] = self.auto_save_line.text()
        return saved_state

    # TODO split to worker thread
    def start_clicked(self, button):

        if self.auto_save.isChecked():
            if not os.path.exists(os.path.dirname(self.auto_save_line.text())):
                self.notify_exception("The Auto-Save file path does not exist.  Enter a valid file path before continuing")
                return

        if not os.path.exists(self.unit_poly_path_edit.text()):
            self.notify_exception("The unit polygon file path does not exist. Enter a valid file path before continuing")
            return
        if self.unit_poly_path_edit.text()[len(self.unit_poly_path_edit.text()) - 4:] != ".shp":
            self.notify_exception("The unit polygon expected file type is shapefile. Enter a valid file before continuing")
            return

        if not os.path.exists(self.lidar_data_path_edit.text()):
            self.notify_exception("The lidar dataset path provided does not exist. Enter a valid file path before continuing")
            return
        if not os.path.isdir(self.lidar_data_path_edit.text()):
            self.notify_exception("The lidar dataset is expected to be a folder. Enter a valid folder before continuing")
            return

        if self.reference_data_path_edit.text() != "":
            if self.reference_data_path_edit.text()[len(self.reference_data_path_edit.text()) - 4:] != ".csv":
                self.notify_exception("The reference dataset expected file type is csv. Enter a valid csv before continuing")
                return

        q = Queue()
        sys.stdout = WriteStream(q)
        print("Processing started:")
        def append_text(text):
            self.text_output.append(text)

        self.stream_thread = QThread()
        self.r = Receiver(q)
        self.r.text_stream.connect(append_text)
        self.r.moveToThread(self.stream_thread)
        self.stream_thread.started.connect(self.r.run)
        self.stream_thread.start()

        def restore_stdout():
            sys.stdout = sys.__stdout__
            print("stdout restored")

        self.ps_thread = QThread()
        project_settings = None
        self.w = Worker(
            self.prj_name_edit.text(),
            self.unit_poly_path_edit.text(),
            self.reference_data_path_edit.text(),
            self.lidar_data_path_edit.text(),
            self.unit_name_edit.text(),
            self.threads_edit.value(),
            self.prop_table_path,
            self.dll_path
        )
        self.w.moveToThread(self.ps_thread)
        self.ps_thread.started.connect(self.w.run)
        self.w.finished.connect(self.ps_thread.quit)
        self.w.finished.connect(self.w.deleteLater)
        self.ps_thread.finished.connect(self.ps_thread.deleteLater)
        def cleanup_stream_thread():
            print("cstreamthred1")
            if self.stream_thread.isRunning():
                print("cstreamthred2")
                self.r.stop()
                self.stream_thread.quit()
                self.r.deleteLater()
                self.stream_thread.deleteLater()
                print("cstreamthred3")

        self.w.finished.connect(cleanup_stream_thread)
        self.w.finished.connect(lambda: print(project_settings))

        self.ps_thread.start()

        self.start_button.setEnabled(False)
        self.w.finished.connect(
            lambda: self.start_button.setEnabled(True)
        )
        self.w.finished.connect(self.start_gamingactivity)
        self.w.finished.connect(
            lambda: print("projectsettings done")
        )
        self.w.finished.connect(restore_stdout)

    def start_gamingactivity(self, project_settings):
        if isinstance(project_settings, ProjectSettings):
            try:
                autosave_path = self.auto_save_line.text() if self.auto_save.isChecked() else None
                Activity.Start_Activity(
                    GamingActivity,
                    None,
                    {"ProjectSettings": project_settings, "Autosave_path": autosave_path},
                    WindowMode.Sibling,
                    dll_path=self.dll_path
                )
                self.stop()
            except Exception as e:
                txt = str(e) + "\n\nDebugging Traceback:\n" + traceback.format_exc()
                self.notify_exception(txt)
        else:
            self.notify_exception(project_settings)


    def save_changed(self, state):
        if state == 2:
            self.auto_save_line.setEnabled(True)
        else:
            self.auto_save_line.setEnabled(False)

    def save_clicked(self):
        if self.save_file_location != "":
            if os.path.isfile(self.save_file_location):
                with open(self.save_file_location, 'wb') as fp:
                    pickle.dump(self.save(), fp)
                self.notify_save_success()
                return
        self.save_as_clicked()

    def save_as_clicked(self):
        SaveStateActivity.prompt_and_save(SaveStateActivity, self.save())
        self.notify_save_success()

    @staticmethod
    def notify_save_success():
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText("Save Successful!")
        msg.setWindowTitle("Result")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    @staticmethod
    def notify_exception(text):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setText("There were errors with some of the parameters provided.\n"
                    "Often a description of the error is provided below, "
                    "with further in depth traceback reporting after. "
                    "If you would like to report this error or believe you have encountered a bug, "
                    "Please email this text to bnbg@uw.edu. Thank you:")
        msg.setInformativeText(text)
        msg.setWindowTitle("Parameter error")
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()


class WriteStream:
    def __init__(self, queue):
        self.queue = queue

    def write(self, text):
        self.queue.put(text)

    def flush(self):
        pass


class Receiver(QObject):
    text_stream = pyqtSignal(str)

    def __init__(self, queue, *args, **kwargs):
        QObject.__init__(self, *args, **kwargs)
        self.queue = queue
        self.running = False

    def stop(self):
        self.running = False
    @pyqtSlot()
    def run(self):
        self.running = True
        while self.running:
            if not self.queue.empty():
                text = self.queue.get()
                if len(text):
                    sys.__stdout__.write(text)
                    self.text_stream.emit(text)
            else:
                QThread.msleep(50)


class Worker(QObject):
    finished = pyqtSignal(object)

    def __init__(self, prj_name, unit, ref, lidar, unit_name, threads, prop_table, dll, *args, **kwargs):
        QObject.__init__(self, *args, **kwargs)

        self.prj_name = prj_name
        self.unit = unit
        self.ref = ref
        self.lidar = lidar
        self.unit_name = unit_name
        self.threads = threads
        self.prop_table = prop_table
        self.dll = dll

    def run(self):
        try:
            ps = ProjectSettings(
                self.prj_name,
                self.unit,
                self.ref,
                self.lidar,
                self.unit_name,
                self.threads,
                self.prop_table,
                self.dll
            )
            self.finished.emit(ps)
        except Exception as e:
            txt = str(e) + "\n\nDebugging Traceback:\n" + traceback.format_exc()
            self.finished.emit(txt)

