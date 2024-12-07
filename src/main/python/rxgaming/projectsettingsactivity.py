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
from PyQt5.QtWidgets import (QFormLayout, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QMessageBox, QCheckBox,
                             QSpinBox)
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
                msg = QMessageBox()
                msg.setText("The Auto-Save file path does not exist.  Enter a valid file path before continuing")
                msg.setWindowTitle("Auto-Save error")
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec_()
                return

        try:
            project_settings = ProjectSettings(self.prj_name_edit.text(), self.unit_poly_path_edit.text(),
                                               self.reference_data_path_edit.text(), self.lidar_data_path_edit.text(),
                                               self.unit_name_edit.text(), self.threads_edit.value(),
                                               self.prop_table_path, self.dll_path)
            print("projectsettings done")
        except Exception as e:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("There were errors with some of the parameters provided:")
            msg.setInformativeText(traceback.format_exc())
            msg.setWindowTitle("Parameter error")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
            return

        autosave_path = self.auto_save_line.text() if self.auto_save.isChecked() else None
        Activity.Start_Activity(GamingActivity, None, {"ProjectSettings": project_settings, "Autosave_path": autosave_path},
                                WindowMode.Sibling, dll_path=self.dll_path)
        self.stop()

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

