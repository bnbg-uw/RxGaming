import os.path
import pickle
import sys
import traceback
from PyQt5.QtWidgets import QFormLayout, QVBoxLayout, QLineEdit, QPushButton, QMessageBox, QCheckBox
from QtWidgets import QFileSelectionLineEdit
from activity import Activity, WindowMode, SaveStateActivity
from projectsettings import ProjectSettings
from gamingactivity import GamingActivity


# ProjectSettingsActivity class
#  Handles user input for getting required files, then kicks off GamingActivity
class ProjectSettingsActivity(Activity):
    def on_start(self, saved_state, prop_table_path, dll_path):
        # MCS Prop table path:
        self.prop_table_path = prop_table_path
        self.dll_path = dll_path

        # Initialize UI
        self.prj_name_edit = QLineEdit()
        self.prj_poly_path_edit = QFileSelectionLineEdit(filter="ESRI Shapefile (*.shp)")
        self.unit_poly_path_edit = QFileSelectionLineEdit(filter="ESRI Shapefile (*.shp)")
        self.reference_data_path_edit = QFileSelectionLineEdit(filter="Comma-separated values (*.csv)")
        self.lidar_data_path_edit = QFileSelectionLineEdit(file_type=QFileSelectionLineEdit.FileType.Directory)
        self.allometry_coefficients_edit = QLineEdit()
        #self.allometry_coefficients_edit.setText("-5.1, 2.5")
        self.unit_name_edit = QLineEdit()
        self.auto_save = QCheckBox()
        self.auto_save_line = QFileSelectionLineEdit(filter="*.dat", new_file=True)
        self.auto_save_line.setEnabled(False)

        self.form_layout = QFormLayout()
        self.form_layout.addRow("Project name", self.prj_name_edit)
        self.form_layout.addRow("Project area polygon shapefile", self.prj_poly_path_edit)
        self.form_layout.addRow("Treatment unit polygon shapefile", self.unit_poly_path_edit)
        self.form_layout.addRow("Reference data table", self.reference_data_path_edit)
        self.form_layout.addRow("Lidar data root folder", self.lidar_data_path_edit)
        #self.form_layout.addRow("Allometric Coefficients", self.allometry_coefficients_edit)
        self.form_layout.addRow("Unit name field:", self.unit_name_edit)
        self.form_layout.addRow("Auto-Save?", self.auto_save)
        self.form_layout.addRow("Save File Location", self.auto_save_line)
        self.start_button = QPushButton("Start")

        self.outer_layout = QVBoxLayout()
        self.outer_layout.addLayout(self.form_layout)
        self.outer_layout.addWidget(self.start_button)
        self.window.setLayout(self.outer_layout)
        self.window.setWindowTitle("Project settings. Rxgaming tool version: " + Activity.version)
        
        # Setup callbacks
        self.auto_save.stateChanged.connect(self.save_changed)
        self.start_button.clicked.connect(self.start_clicked)
        
        # Load saved state
        if "ProjectSettingsActivity.prj_name" in saved_state:
            self.prj_name_edit.setText(saved_state["ProjectSettingsActivity.prj_name"])
        if "ProjectSettingsActivity.prj_poly_path" in saved_state:
            self.prj_poly_path_edit.setText(saved_state["ProjectSettingsActivity.prj_poly_path"])
        if "ProjectSettingsActivity.unit_poly_path" in saved_state:
            self.unit_poly_path_edit.setText(saved_state["ProjectSettingsActivity.unit_poly_path"])
        if "ProjectSettingsActivity.reference_data_path" in saved_state:
            self.reference_data_path_edit.setText(saved_state["ProjectSettingsActivity.reference_data_path"])
        if "ProjectSettingsActivity.lidar_data_path" in saved_state:
            self.lidar_data_path_edit.setText(saved_state["ProjectSettingsActivity.lidar_data_path"])
        #if "ProjectSettingsActivity.allometry_coefficients" in saved_state:
        #    self.allometry_coefficients_edit.setText(saved_state['ProjectSettingsActivity.allometry_coefficients'])
        if "ProjectSettingsActivity.unit_name" in saved_state:
            self.unit_name_edit.setText(saved_state['ProjectSettingsActivity.unit_name'])
        if "ProjectSettingsActivity.auto_save" in saved_state:
            self.auto_save.setChecked(saved_state['ProjectSettingsActivity.auto_save'])
        if "ProjectSettingsActivity.auto_save_line" in saved_state:
            self.auto_save_line.setText(saved_state['ProjectSettingsActivity.auto_save_line'])

    def save(self):
        saved_state = dict()
        saved_state["ProjectSettingsActivity.prj_name"] = self.prj_name_edit.text()
        saved_state["ProjectSettingsActivity.prj_poly_path"] = self.prj_poly_path_edit.text()
        saved_state["ProjectSettingsActivity.unit_poly_path"] = self.unit_poly_path_edit.text()
        saved_state["ProjectSettingsActivity.reference_data_path"] = self.reference_data_path_edit.text()
        saved_state["ProjectSettingsActivity.lidar_data_path"] = self.lidar_data_path_edit.text()
        #saved_state["ProjectSettingsActivity.allometry_coefficients"] = self.allometry_coefficients_edit.text()
        saved_state["ProjectSettingsActivity.unit_name"] = self.unit_name_edit.text()
        saved_state["ProjectSettingsActivity.auto_save"] = self.auto_save.isChecked()
        saved_state['ProjectSettingsActivity.auto_save_line'] = self.auto_save_line.text()
        return saved_state

    # TODO split to worker thread
    def start_clicked(self, button):
        # TODO reference dbs
        if self.auto_save.isChecked():
            if not os.path.exists(os.path.dirname(self.auto_save_line.text())):
                msg = QMessageBox()
                msg.setText("The Auto-Save file path does not exist.  Enter a valid file path before continuing")
                msg.setWindowTitle("Auto-Save error")
                msg.setStandardButtons(QMessageBox.Ok)
                msg.exec_()
                return

        if len(self.allometry_coefficients_edit.text()) == 0:
            allom_list = []
        else:
            allom_list = [float(i) for i in ''.join(self.allometry_coefficients_edit.text().split()).split(',')]

        try:
            project_settings = ProjectSettings(self.prj_name_edit.text(), self.prj_poly_path_edit.text(),
                                               self.unit_poly_path_edit.text(), self.reference_data_path_edit.text(),
                                               self.lidar_data_path_edit.text(), self.unit_name_edit.text(), allom_list,
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
