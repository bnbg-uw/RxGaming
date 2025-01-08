"""
    Copyright (C) 2024  University of Washington
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.

Bryce Bartl-Geller
University of Washington Forest Resilience Lab
12/6/2024

gamingactivity.py
This file is the UI for the "gaming" side of the tool, after project settings processing has been completed.
"""

# General
import numpy as np
import os

# QT
from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton, QTabWidget, QWidget, QLabel, QListView,
                             QDataWidgetMapper, QLineEdit, QComboBox, QFileDialog, QSizePolicy, QAction, QApplication,
                             QMessageBox)
from PyQt5.QtCore import QAbstractTableModel, Qt, QVariant, QModelIndex
from PyQt5 import QtPrintSupport, QtGui

# RxGaming
from activity import Activity, SaveStateActivity
from QtWidgets import SliderWithValue, QMainWindowRx  # custom widgets

# IO
from shapely.geometry import mapping, Point, Polygon as ShpPolygon  # for displaying on ref plots & for writing .shp
import raster  # our custom raster wrapper-not full featured but holds data in the way we expect and easy to send to c++
import fiona  # Read/write .shp and crs
import fiona.crs
import rasterio
import csv  # for writing treelists
import PIL.Image  # To read the temp file for the geotifs
from PIL import ImageDraw, ImageFont  # To add text to exported geotifs
import tempfile  # For writing geotifs
import pickle  # For saving the whole program

# Graphics
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas  # Allows us to redner plots in widgets.
import descartes  # For polygons on ref plots
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Polygon as MplPolygon
import matplotlib.ticker as ticker  # custom axis ticks
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm  # coloring the displayed raster
from scipy.spatial import ConvexHull  # polygons on reference plots
from scipy.stats import gaussian_kde  # report pages
import seaborn # reference 2d KDE plots


# GamingActivity class
# Handles user input and interactions with the projectsettings and surfaces data and treatment results
class GamingActivity(Activity):
    def on_start(self, saved_state, **kwargs):
        self.dll_path = kwargs['dll_path']  #find the link to the dll that does the treatment operations
        # Set up the mainwindow
        self.set_window(QMainWindowRx())
        self.central_widget = QWidget()
        self.window.setCentralWidget(self.central_widget)

        # set up the window settings, saved state, and tabs.
        self.saved_state = saved_state
        project_settings = saved_state["ProjectSettings"]
        self.tab_widget = Tabs(saved_state, self.dll_path)
        self.layout  = QVBoxLayout()
        self.layout.addWidget(self.tab_widget)
        self.central_widget.setLayout(self.layout)
        self.window.setWindowTitle(project_settings.get_name() + " Gaming. Rxgaming tool version: " + Activity.version)
        self.window.setGeometry(200, 200, 1800, 1000)
        self.window.setMinimumSize(1800, 1000)
        self.window.showMaximized()

        # Set up the menu bar
        self.save_action = QAction("&Save")
        self.save_action.setShortcut("Ctrl+S")
        self.save_action.triggered.connect(self.menu_save)

        self.save_as_action = QAction("&Save as")
        self.save_as_action.setShortcut("Ctrl+Shift+S")
        self.save_as_action.triggered.connect(self.menu_save_as)

        self.exit_action = QAction("&Exit")
        self.exit_action.setShortcut("Ctrl+Q")
        self.exit_action.triggered.connect(self.menu_exit)

        self.export_tif_action = QAction("&Export window- geotiff (\".tif\")")
        self.export_tif_action.triggered.connect(self.export_tif)

        self.export_rasters_action = QAction("&Export raster data (\".tif\")")
        self.export_rasters_action.triggered.connect(self.export_rasters)

        self.export_features_action = QAction("&Export shapefile data (\".shp\")")
        self.export_features_action.triggered.connect(self.export_features)

        self.export_treelists_action = QAction("&Export treelists (\".csv\")")
        self.export_treelists_action.triggered.connect(self.export_treelists)

        self.print_action = QAction("&Print")
        self.print_action.triggered.connect(self.print)

        self.license_action = QAction("&License")
        self.license_action.triggered.connect(self.license)

        self.main_menu = self.window.menuBar()
        self.file_menu = self.main_menu.addMenu("File")
        self.export_menu = self.file_menu.addMenu("Export")
        #self.file_menu.addAction(self.print_action)

        self.file_menu.addAction(self.save_action)
        self.file_menu.addAction(self.save_as_action)
        self.file_menu.addAction(self.license_action)
        self.file_menu.addAction(self.exit_action)

        self.export_menu.addAction(self.export_rasters_action)
        self.export_menu.addAction(self.export_features_action)
        self.export_menu.addAction(self.export_tif_action)
        self.export_menu.addAction(self.export_treelists_action)

        # If the user specified they wanted to autosave after the initial processing, do that here.
        if "Autosave_path" in saved_state:
            if saved_state['Autosave_path'] is not None:
                SaveStateActivity.write_file(file_path=self.saved_state['Autosave_path'], saved_state=self.save())
                self.tab_widget.project_settings.prj_area.dePickle(self.dll_path)

    # Implement the method called when the activity closes.
    def save(self):
        saved_state = dict()
        self.tab_widget.project_settings.prj_area.prepToPickle()

        saved_state['ProjectSettings'] = self.tab_widget.project_settings
        saved_state['GamingActivity.tabs.rx_units'] = self.tab_widget.rx_units
        saved_state['GamingActivity.tabs.decision_spaces'] = self.tab_widget.decision_spaces
        saved_state['GamingActivity.tabs.hills'] = self.tab_widget.hills
        saved_state['GamingActivity.tabs.model_index_row'] = self.tab_widget.stand_tab_list_view.currentIndex().row()
        saved_state['GamingActivity.tabs.viewmode'] = self.tab_widget.raster_viewmode.currentIndex()
        saved_state['GamingActivity.tabs.treatmethod'] = self.tab_widget.treatment_method.currentIndex()
        saved_state['GamingActivity.tabs.cut_range'] = self.tab_widget.cut_range.value()
        saved_state['LastActivity'] = type(self)
        if 'save_file_location' in self.saved_state:
            saved_state['save_file_location'] = self.saved_state['save_file_location']
        return saved_state

    # action for a user instantiated save.
    def menu_save(self):
        if 'save_file_location' in self.saved_state:
            print(self.saved_state['save_file_location'])
            if type(self.saved_state['LastActivity']) is not type(self):
                msg = QMessageBox()
                msg.setText("The file you are saving to is a project settings instance, would you like to overwrite?")
                msg.setWindowTitle("Overwrite save state?")
                msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
                msg.buttonClicked.connect(self.menu_save_connect)
                msg.exec_()
            else:
                if os.path.isfile(self.saved_state['save_file_location']):
                    self.menu_save_success()
                else:
                    self.menu_save_as()
        else:
            self.menu_save_as()

    def menu_save_success(self):
        with open(self.saved_state['save_file_location'], 'wb') as fp:
            pickle.dump(self.save(), fp)
        self.tab_widget.project_settings.prj_area.dePickle(self.dll_path)

    def menu_save_connect(self, i):
        print(i.text())
        if i.text() == "&Yes":
            self.menu_save_success()
        else:
            self.menu_save_as()

    def menu_save_as(self):
        SaveStateActivity.prompt_and_save(SaveStateActivity, self.save())
        self.tab_widget.project_settings.prj_area.dePickle(self.dll_path)

    def menu_exit(self):
        self.stop()

    def export_tif(self):
        self.tab_widget.export_tif()

    def export_rasters(self):
        self.tab_widget.export_raster()

    def export_features(self):
        self.tab_widget.export_features()

    def export_treelists(self):
        self.tab_widget.export_treelist()

    def print(self):
        self.tab_widget.print()

    @staticmethod
    def license():
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText("""
            Copyright (C) 2024  University of Washington\n
            This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.\n
            This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
            You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.\n
            """)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()


# The majority of the UI sits inside these tabs.
class Tabs(QTabWidget):
    def __init__(self, saved_state, dll_path):
        super(Tabs, self).__init__(None)
        self.layout = QVBoxLayout()

        # Rx decision space polygons will be held here eventually
        self.decision_spaces = []

        # Initialize tabs
        self.tabs = QTabWidget()
        self.stand_tab = QWidget()
        self.land_tab = QWidget()

        # Add tabs
        self.tabs.addTab(self.stand_tab, 'Stand View')
        self.tabs.addTab(self.land_tab, 'Landscape View')

        # Create stand view tab
        self.stand_tab.layout = QHBoxLayout()

        self.stand_tab_left_layout = QVBoxLayout()
        self.stand_tab_middle_layout = QVBoxLayout()
        self.stand_tab_right_layout = QVBoxLayout()

        self.stand_tab_unit_label = QLabel("UNITS")
        self.stand_tab_list_view = QListView()
        self.stand_tab_info_widget = StructureInfo()

        self.stand_tab_unit_label.setMaximumWidth(200)
        self.stand_tab_list_view.setMaximumWidth(200)
        self.stand_tab_info_widget.setMaximumWidth(200)

        self.stand_tab_left_layout.addWidget(self.stand_tab_unit_label)
        self.stand_tab_left_layout.addWidget(self.stand_tab_list_view)
        self.stand_tab_left_layout.addWidget(self.stand_tab_info_widget)

        # Set up treatment report widgets figures.
        #    these are various labels needed.
        self.report_label = QLabel("Treatment Report")
        self.current_label = QLabel("Current")
        self.displayed_label = QLabel("Post-Treatment")
        self.target_label = QLabel("Target")
        self.report_label.setStyleSheet("font: 24pt;")
        self.current_label.setStyleSheet('font: 16pt;')
        self.displayed_label.setStyleSheet('font: 16pt;')
        self.target_label.setStyleSheet('font: 16pt;')

        #  These are the kernel density estimates that show info about the selected stand.
        self.current_fig_ba = Figure()
        self.current_canvas_ba = FigureCanvas(self.current_fig_ba)
        self.report_current_ax_ba = self.current_fig_ba.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.current_fig_ba.subplots_adjust(left=0.1, right=0.97, top=.85, bottom=0.25)
        self.current_fig_mcs = Figure()
        self.current_canvas_mcs = FigureCanvas(self.current_fig_mcs)
        self.report_current_ax_mcs = self.current_fig_mcs.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.current_fig_mcs.subplots_adjust(left=0.1, right=0.97, top=.85, bottom=0.25)
        self.displayed_fig_ba = Figure()
        self.displayed_canvas_ba = FigureCanvas(self.displayed_fig_ba)
        self.displayed_ax_ba = self.displayed_fig_ba.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.displayed_fig_ba.subplots_adjust(left=0.1, right=0.97, top=.85, bottom=0.25)
        self.displayed_fig_mcs = Figure()
        self.displayed_canvas_mcs = FigureCanvas(self.displayed_fig_mcs)
        self.displayed_ax_mcs = self.displayed_fig_mcs.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.displayed_fig_mcs.subplots_adjust(left=0.1, right=0.97, top=.85, bottom=0.25)
        self.displayed_mcs_prop = QLabel("")
        self.displayed_mcs_prop.setStyleSheet('font: 16pt;')

        # Set up treatment window layouts and widgets.
        self.treat_display = QWidget()
        self.treat_display.setVisible(False)
        self.treat_display_layout = QGridLayout()
        self.treat_display.setLayout(self.treat_display_layout)

        # Put the previously created widgets in to one layout to be displayed together.
        self.treat_display_layout.addWidget(self.report_label, 0, 1)
        self.treat_display_layout.addWidget(self.current_label, 1, 0)
        self.treat_display_layout.addWidget(self.displayed_label, 2, 0)
        self.treat_display_layout.addWidget(self.target_label, 3, 0)
        self.treat_display_layout.addWidget(self.current_canvas_ba, 1, 1)
        self.treat_display_layout.addWidget(self.displayed_canvas_ba, 2, 1)
        self.treat_display_layout.addWidget(self.displayed_mcs_prop, 3, 2)
        self.treat_display_layout.addWidget(self.current_canvas_mcs, 1, 2)
        self.treat_display_layout.addWidget(self.displayed_canvas_mcs, 2, 2)
        self.treat_display_layout.setRowStretch(0, 0)
        self.treat_display_layout.setRowStretch(1, 2)
        self.treat_display_layout.setRowStretch(2, 2)
        self.treat_display_layout.setRowStretch(3, 1)
        self.displayed_mcs_prop.setMinimumHeight(10)
        self.displayed_canvas_mcs.setMinimumHeight(10)

        # set up cut trees report window
        self.cut_trees_page_label = QLabel("Cut trees info")
        self.cut_trees_page_label.setStyleSheet('font: 16pt;')
        self.cut_trees_label = QLabel("Cut Trees:")
        self.cut_trees_label.setStyleSheet('font: 16pt;')
        self.cut_trees_fig = Figure()
        self.cut_trees_canvas = FigureCanvas(self.cut_trees_fig)
        self.cut_trees_ax = self.cut_trees_fig.add_subplot(111)

        self.cut_trees_display = QWidget()
        self.cut_trees_display.setVisible(False)
        self.cut_trees_display_layout = QGridLayout()
        self.cut_trees_display.setLayout(self.cut_trees_display_layout)

        self.cut_trees_display_layout.addWidget(self.cut_trees_page_label, 0, 1)
        self.cut_trees_display_layout.addWidget(self.cut_trees_label, 1, 0)
        self.cut_trees_display_layout.addWidget(self.cut_trees_canvas, 1, 1)

        self.cut_trees_display_layout.setColumnStretch(1, 1)
        self.cut_trees_display_layout.setRowStretch(1, 1)

        # Set up raster report window.
        self.treatment = Treatment(None, None, None, None)  # For storing treat info upon multi threading.
        self.raster_widget = QWidget()
        self.raster_layout = QVBoxLayout()
        self.raster_widget.setLayout(self.raster_layout)
        self.raster_figure = Figure()
        self.raster_canvas = FigureCanvas(self.raster_figure)
        self.raster_ax = self.raster_figure.add_subplot(111)
        self.raster_cb = None  # For saving the color bar
        self.raster_layout.addWidget(self.raster_canvas)

        self.raster_button = QPushButton("Show Treatment")
        self.raster_button.setCheckable(True)
        self.raster_button.clicked[bool].connect(self.raster_button_clicked)
        self.raster_button.setStyleSheet("QPushButton:enabled {color: black;}\n\
                                         QPushButton:checked { background-color: rgb(80, 80, 80); \
                                         border: none;\
                                         color: white}")

        self.cut_range = SliderWithValue(Qt.Horizontal)
        self.cut_range.setMinimum(0)
        self.cut_range.setMaximum(120)
        self.cut_range.setValue(30)

        self.raster_viewmode = QComboBox()
        self.raster_viewmode.addItems(["Canopy Model", "Basins", "Clumps"])
        self.raster_viewmode.currentIndexChanged.connect(self.update_raster_canvas)

        self.treatment_method = QComboBox()
        self.treatment_method.addItems(["Add Clumps", 'DBH Thin', "Sum Squares"])
        self.treatment_method.setVisible(False)

        # Page increment and decrement.
        self.page_counter = 0
        self.next_page = QPushButton("----->")
        self.next_page.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.next_page.clicked.connect(self.increment_page)
        self.prev_page = QPushButton("<-----")
        self.prev_page.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.prev_page.clicked.connect(self.decrement_page)

        # Set up the ui at the bottom of the stand view page.
        self.stand_tab_settings = QWidget()
        self.stand_tab_settings.setMaximumHeight(130)
        self.stand_tab_settings_layout = QGridLayout()
        self.stand_tab_settings.setLayout(self.stand_tab_settings_layout)
        self.stand_tab_settings_wrapper = QVBoxLayout()

        self.cut_label = QLabel("DBH Cutoff (Inches):")
        self.cut_label.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.viewmode_label = QLabel("View Mode:")
        self.viewmode_label.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.treatment_label = QLabel("Thinning Method:")
        self.treatment_label.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.treatment_label.setVisible(False)

        # This is the UI that stays constant between pages.
        self.stand_tab_settings_layout.addWidget(self.prev_page, 0, 2)
        self.stand_tab_settings_layout.addWidget(self.next_page, 0, 3)
        self.stand_tab_settings_layout.addWidget(self.cut_label, 1, 0)
        self.stand_tab_settings_layout.addWidget(self.cut_range, 1, 1, 1, 4)
        self.stand_tab_settings_layout.addWidget(self.viewmode_label, 2, 0)
        self.stand_tab_settings_layout.addWidget(self.raster_viewmode, 2, 1, 1, 4)
        self.stand_tab_settings_layout.addWidget(self.treatment_label, 3, 0)
        self.stand_tab_settings_layout.addWidget(self.treatment_method, 3, 1, 1, 4)
        self.stand_tab_settings_layout.addWidget(self.raster_button, 4, 1, 1, 3)

        self.stand_tab_settings_wrapper.addWidget(self.stand_tab_settings)

        self.stand_tab_middle_layout.addWidget(self.raster_widget)
        self.stand_tab_middle_layout.addWidget(self.treat_display)
        self.stand_tab_middle_layout.addWidget(self.cut_trees_display)
        self.stand_tab_middle_layout.addLayout(self.stand_tab_settings_wrapper)

        self.stand_tab.layout.addLayout(self.stand_tab_left_layout)
        self.stand_tab.layout.addLayout(self.stand_tab_middle_layout)
        self.stand_tab.layout.addLayout(self.stand_tab_right_layout)
        self.stand_tab.setLayout(self.stand_tab.layout)

        # Create land view tab
        self.land_tab.layout = QVBoxLayout()
        self.reference_figure = Figure()
        self.ba_ax = self.reference_figure.add_subplot(1, 3, 1)
        self.mcs_ax = self.reference_figure.add_subplot(1, 3, 2)
        self.cc_ax = self.reference_figure.add_subplot(1, 3, 3)
        self.reference_canvas = FigureCanvas(self.reference_figure)

        # QT window layout and stuff
        self.land_tab.layout.addWidget(self.reference_canvas)
        self.land_tab.setLayout(self.land_tab.layout)

        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        # Trigger loading and prepping data from save state, and prepare for threading.
        self.load_and_prep_data(saved_state, dll_path)

        ba_label = r"$feet^2\ ac^{-1}$"
        tpa_label = r"$ trees\ ac^{-1}$"

        self.ba_ax.set_title("Basal Area")
        self.ba_ax.set_ylabel("Basal Area (%s)" % ba_label)
        self.ba_ax.set_xlabel("Density (%s)" % tpa_label)

        self.mcs_ax.set_title("Mean Clump Size")
        self.mcs_ax.set_ylabel("Mean Clump Size (n trees, log)")
        self.mcs_ax.set_xlabel("Density (%s)" % tpa_label)

        self.cc_ax.set_title("Canopy Cover")
        self.cc_ax.set_ylabel("Canopy Cover (%)")
        self.cc_ax.set_xlabel("Density (Trees %s)" % tpa_label)

        #set model for stand tab list view
        self.stand_tab_list_view.setModel(self.model)
        self.stand_tab_info_widget.set_model(self.model)
        self.stand_tab_list_view.selectionModel().currentChanged.connect(self.stand_tab_info_widget.set_selection)
        self.stand_tab_list_view.selectionModel().currentChanged.connect(self.reset_page)
        self.stand_tab_list_view.selectionModel().currentChanged.connect(self.reset_button)
        self.stand_tab_list_view.selectionModel().currentChanged.connect(self.update_raster_canvas)

        # set up plots for land tab
        self.unit_ba_points, = self.ba_ax.plot(self.current_tpa, self.current_ba, 'b^')
        self.unit_mcs_points, = self.mcs_ax.plot(self.current_tpa, self.current_mcs, 'b^')
        self.unit_cc_points, = self.cc_ax.plot(self.current_tpa, self.current_cc, 'b^')

        # if ref data is provided, 2d KDE plots of ref data.
        tmprefmcs = 0
        if self.ref_tpa is not None:
            seaborn.kdeplot(ax=self.ba_ax, x=self.ref_tpa[self.ref_ind['ba']], y=self.ref_ba[self.ref_ind['ba']],
                            cmap="Oranges", fill=True, bw_adjust=0.5)
            seaborn.kdeplot(ax=self.mcs_ax, x=self.ref_tpa[self.ref_ind['mcs']], y=self.ref_mcs[self.ref_ind['mcs']],
                            cmap="Oranges", fill=True, bw_adjust=0.5)
            seaborn.kdeplot(ax=self.cc_ax, x=self.ref_tpa[self.ref_ind['cc']], y=self.ref_cc[self.ref_ind['cc']],
                            cmap="Oranges", fill=True, bw_adjust=0.5)
            tmprefmcs = np.max(self.ref_mcs[self.ref_ind['mcs']])

        # Log mcs axes if there are large values.
        if max(np.max(self.current_mcs), tmprefmcs) >= 10:
            self.mcs_ax.set_yscale("log")
            self.mcs_ax.set_ylim(bottom=1)


        #These are points that will be modified based on rx targets and treatment results values
        self.target_points = [self.ba_ax.plot(0, 0, "gs")[0], self.mcs_ax.plot(0, 0, "gs")[0], self.cc_ax.plot(0, 0, "gs")[0]]
        self.treated_points = [self.ba_ax.plot(0, 0, "mo")[0], self.mcs_ax.plot(0, 0, "mo")[0], self.cc_ax.plot(0, 0, "mo")[0]]

        # And these are arrows that span between the points.
        self.cur_targ_arrow = [self.ba_ax.arrow(0, 0, 0, 0), self.mcs_ax.arrow(0, 0, 0, 0), self.cc_ax.arrow(0, 0, 0, 0)]
        self.targ_treated_arrow = [self.ba_ax.arrow(0, 0, 0, 0), self.mcs_ax.arrow(0, 0, 0, 0), self.cc_ax.arrow(0, 0, 0, 0)]

        # But they start invisible since we only display them when hovering.
        for i in range(3):
            self.cur_targ_arrow[i].set_visible(False)
            self.targ_treated_arrow[i].set_visible(False)

        # this renders the changes we made.
        self.reference_canvas.draw_idle()

        #ANNOTATIONS reference:
        # Set up annotations and hide them so they show up when hovering over them.
        # Connect an action that monitors the mouse position so they show up only when the mouse is hovering over.
        self.ba_annot = self.ba_ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                                            bbox=dict(boxstyle="round", fc='w'),
                                            arrowprops=dict(arrowstyle='->'))
        self.ba_annot.set_visible(True)
        self.reference_canvas.mpl_connect("motion_notify_event", self.ba_hover)

        self.mcs_annot = self.mcs_ax.annotate("", xy=(0, 0), xytext=(-20, 20),
                                              textcoords="offset points",
                                              bbox=dict(boxstyle="round", fc='w'),
                                              arrowprops=dict(arrowstyle='->'))
        self.mcs_annot.set_visible(False)
        self.reference_canvas.mpl_connect("motion_notify_event", self.mcs_hover)

        self.cc_annot = self.cc_ax.annotate("", xy=(0, 0), xytext=(-20, 20),
                                            textcoords="offset points",
                                            bbox=dict(boxstyle="round", fc='w'),
                                            arrowprops=dict(arrowstyle='->'))
        self.cc_annot.set_visible(False)
        self.reference_canvas.mpl_connect("motion_notify_event", self.cc_hover)

        self.reference_canvas.mpl_connect("button_press_event", self.rx_unit_pick)

        # Creating the legend for the points
        handles = [Patch(facecolor=plt.cm.Oranges(100)),
                   self.unit_cc_points, self.target_points[0], self.treated_points[0]]
        self.reference_figure.legend(handles, ("Reference", "Units", "Targets", "Treated"), loc='right')

        # I forget why but these have to be set invisible after the handles and legend are created...
        [x.set_visible(False) for x in self.target_points]
        [x.set_visible(False) for x in self.treated_points]


        # Load any info that is stored in the save state.
        if 'GamingActivity.tabs.model_index_row' in saved_state:
            row = saved_state['GamingActivity.tabs.model_index_row']
            index = self.model.index(row, 0)
            self.stand_tab_list_view.setCurrentIndex(index)
        else:
            index = self.model.index(0, 0)
            self.stand_tab_list_view.setCurrentIndex(index)

        if 'GamingActivity.tabs.viewmode' in saved_state:
            self.raster_viewmode.setCurrentIndex(saved_state['GamingActivity.tabs.viewmode'])
        if 'GamingActivity.tabs.treatmethod' in saved_state:
            self.treatment_method.setCurrentIndex(saved_state['GamingActivity.tabs.treatmethod'])
        if 'GamingActivity.tabs.cut_range' in saved_state:
            self.cut_range.setValue(saved_state['GamingActivity.tabs.cut_range'])
        self.reference_canvas.draw_idle()
        self.update_raster_canvas()

    # Take the data from the save state and load.  Or load for the first time if not in the save state.
    def load_and_prep_data(self, saved_state, dll_path):
        self.project_settings = saved_state["ProjectSettings"]

        # When dll type is float, it means that we save the convfactor in its place
        # (this allows the program to shut down the link to the dll correctly because dll's can't be pickled).
        if(type(self.project_settings.prj_area._tao_data.dll) is float):
            print("depickling")
            self.project_settings.prj_area.dePickle(dll_path)
            for rx in self.project_settings.prj_area.get_units():
                self.project_settings.prj_area.get_tao_data().exportRxToDll(rx)

        ref_db = self.project_settings.get_reference_data()

        if "GamingActivity.tabs.rx_units" in saved_state:
            self.rx_units = saved_state['GamingActivity.tabs.rx_units']
        else:
            self.rx_units = self.project_settings.prj_area.get_units()

        if "GamingActivity.tabs.hills" in saved_state:
            self.hills = saved_state["GamingActivity.tabs.hills"]
        else:
            self.hills = [None for _ in range(len(self.rx_units))]

        if "GamingActivity.tabs.model" in saved_state:
            self.model = saved_state["GamingActivity.tabs.model"]
        else:
            self.model = RxUnitTableModel(self.rx_units)

        # Get reference db info
        if ref_db is not None:
            self.ref_tpa = ref_db['tph'].astype(float) / 2.47105
            self.ref_ba = ref_db['ba'].astype(float) * 4.356
            self.ref_mcs = ref_db["mcs"].astype(float)
            self.ref_cc = ref_db['cc'].astype(float) * 100

            # Indices for each metric where tpa and the metric are both not NaN.
            self.ref_ind = {'tpa': np.where(~np.isnan(self.ref_tpa)), 'ba': np.where(~np.isnan(self.ref_ba)),
                            'mcs': np.where(~np.isnan(self.ref_mcs)), 'cc': np.where(~np.isnan(self.ref_cc))}
            for key in self.ref_ind:
                self.ref_ind[key] = np.intersect1d(self.ref_ind[key], self.ref_ind['tpa'])

            # Ref data chull
            # BA
            ref_ba_data = np.array([self.ref_tpa[self.ref_ind['ba']], self.ref_ba[self.ref_ind['ba']]]).transpose()
            ref_ba_chull = ConvexHull(ref_ba_data)
            ref_ba_shp_poly = ShpPolygon(list(zip(ref_ba_data[ref_ba_chull.vertices, 0],
                                                  ref_ba_data[ref_ba_chull.vertices, 1])))

            # MCS
            ref_mcs_data = np.array([self.ref_tpa[self.ref_ind['mcs']], self.ref_mcs[self.ref_ind['mcs']]]).transpose()
            ref_mcs_chull = ConvexHull(ref_mcs_data)
            ref_mcs_shp_poly = ShpPolygon(list(zip(ref_mcs_data[ref_mcs_chull.vertices, 0],
                                                   ref_mcs_data[ref_mcs_chull.vertices, 1])))

            # CC
            ref_cc_data = np.array([self.ref_tpa[self.ref_ind['cc']], self.ref_cc[self.ref_ind['cc']]]).transpose()
            ref_cc_chull = ConvexHull(ref_cc_data)
            ref_cc_shp_poly = ShpPolygon(list(zip(ref_cc_data[ref_cc_chull.vertices, 0],
                                                   ref_cc_data[ref_cc_chull.vertices, 1])))
        else:
            self.ref_tpa = None
            ref_ba_shp_poly = ShpPolygon(((0, 0), (0, -1), (-1, -1), (-1, 0)))
            ref_mcs_shp_poly = ShpPolygon(((0, 0), (0, -1), (-1, -1), (-1, 0)))
            ref_cc_shp_poly = ShpPolygon(((0, 0), (0, -1), (-1, -1), (-1, 0)))

        # Unit info setup
        self.unit_names = [rx.get_name() for rx in self.rx_units]

        self.current_ba = [rx_unit.get_current_structure().ba for rx_unit in self.project_settings.prj_area.get_units()]
        self.current_tpa = [rx_unit.get_current_structure().tpa for rx_unit in self.project_settings.prj_area.get_units()]
        self.current_mcs = [rx_unit.get_current_structure().mcs for rx_unit in self.project_settings.prj_area.get_units()]
        self.current_cc = [rx_unit.get_current_structure().cc for rx_unit in self.project_settings.prj_area.get_units()]

        # Get the decisions spaces for each unit.
        if 'GamingActivity.tabs.decision_spaces' in saved_state:
            self.decision_spaces = saved_state['GamingActivity.tabs.decision_spaces']
        else:
            self.decision_spaces = []
            for rx_unit in self.project_settings.prj_area.get_units():
                print(rx_unit._name)
                x = rx_unit.get_simulated_structures_dll()
                self.decision_spaces.append(x)

        # BA polygons
        self.unit_ba_polygons = []
        self.unit_ba_intersections = []
        for d_space in self.decision_spaces:
            tpa = [ss.tpa for ss in d_space]
            ba = [ss.ba for ss in d_space]
            unit_ba_data = np.array([tpa, ba]).transpose()
            unit_ba_data = unit_ba_data[~np.isnan(unit_ba_data).any(axis=1), :]
            unit_ba_chull = None
            try:
                unit_ba_chull = ConvexHull(unit_ba_data)
            except Exception as e:
                print(unit_ba_data)
                raise e

            poly = ShpPolygon(list(zip(unit_ba_data[unit_ba_chull.vertices, 0],
                                       unit_ba_data[unit_ba_chull.vertices, 1])))
            inter = poly.intersection(ref_ba_shp_poly)


            if poly.wkt == 'GEOMETRYCOLLECTION EMPTY' or poly.wkt == "POLYGON EMPTY" or inter.wkt == "POINT (0 0)":
                poly = MplPolygon([[0, 0]])
            else:
                poly = descartes.PolygonPatch(poly)
            if inter.wkt == 'GEOMETRYCOLLECTION EMPTY' or inter.wkt == "POLYGON EMPTY" or inter.wkt == "POINT (0 0)":
                inter = MplPolygon([[0, 0]])
            else:
                print(inter.wkt)
                print(inter)
                inter = descartes.PolygonPatch(inter)
            self.ba_ax.add_patch(poly)
            self.ba_ax.add_patch(inter)
            poly.set_visible(False)
            inter.set_visible(False)
            poly.set_alpha(0.6)
            poly.set_color("#f28263")
            inter.set_alpha(0.5)
            poly.set_zorder(3)
            inter.set_zorder(3)
            self.unit_ba_polygons.append(poly)
            self.unit_ba_intersections.append(inter)

        # MCS
        self.unit_mcs_polygons = []
        self.unit_mcs_intersections = []
        for d_space in self.decision_spaces:
            tpa, mcs = zip(*((ss.tpa, ss.mcs) for ss in d_space))
            unit_mcs_data = np.array([tpa, mcs]).transpose()
            unit_mcs_data = unit_mcs_data[~np.isnan(unit_mcs_data).any(axis=1), :]
            unit_mcs_chull = ConvexHull(unit_mcs_data)
            poly = ShpPolygon(list(zip(unit_mcs_data[unit_mcs_chull.vertices, 0],
                                       unit_mcs_data[unit_mcs_chull.vertices, 1])))
            inter = poly.intersection(ref_mcs_shp_poly)
            if poly.wkt == 'GEOMETRYCOLLECTION EMPTY' or poly.wkt == "POLYGON EMPTY" or inter.wkt == "POINT (0 0)":
                poly = MplPolygon([[0, 0]])
            else:
                poly = descartes.PolygonPatch(poly)
            if inter.wkt == 'GEOMETRYCOLLECTION EMPTY' or inter.wkt == "POLYGON EMPTY" or inter.wkt == "POINT (0 0)":
                inter = MplPolygon([[0, 0]])
            else:
                inter = descartes.PolygonPatch(inter)
            self.mcs_ax.add_patch(poly)
            self.mcs_ax.add_patch(inter)
            poly.set_visible(False)
            inter.set_visible(False)
            poly.set_alpha(0.6)
            poly.set_color("#f28263")
            inter.set_alpha(0.5)
            poly.set_zorder(3)
            inter.set_zorder(3)
            self.unit_mcs_polygons.append(poly)
            self.unit_mcs_intersections.append(inter)

        # CC
        self.unit_cc_polygons = []
        self.unit_cc_intersections = []
        for d_space in self.decision_spaces:
            tpa, cc = zip(*((ss.tpa, ss.cc) for ss in d_space))
            unit_cc_data = np.array([tpa, cc]).transpose()
            unit_cc_data = unit_cc_data[~np.isnan(unit_cc_data).any(axis=1), :]
            unit_cc_chull = ConvexHull(unit_cc_data)
            poly = ShpPolygon(list(zip(unit_cc_data[unit_cc_chull.vertices, 0],
                                       unit_cc_data[unit_cc_chull.vertices, 1])))
            inter = poly.intersection(ref_cc_shp_poly)
            if poly.wkt == 'GEOMETRYCOLLECTION EMPTY' or poly.wkt == "POLYGON EMPTY" or inter.wkt == "POINT (0 0)":
                poly = MplPolygon([[0, 0]])
            else:
                poly = descartes.PolygonPatch(poly)
            if inter.wkt == 'GEOMETRYCOLLECTION EMPTY' or inter.wkt == "POLYGON EMPTY" or inter.wkt == "POINT (0 0)":
                inter = MplPolygon([[0, 0]])
            else:
                inter = descartes.PolygonPatch(inter)
            self.cc_ax.add_patch(poly)
            self.cc_ax.add_patch(inter)
            poly.set_visible(False)
            inter.set_visible(False)
            poly.set_alpha(0.6)
            poly.set_color("#f28263")
            inter.set_alpha(0.5)
            poly.set_zorder(3)
            inter.set_zorder(3)
            self.unit_cc_polygons.append(poly)
            self.unit_cc_intersections.append(inter)

    # Move to the next page view on the treatment view tab.
    def increment_page(self):
        if self.page_counter == 2:
            self.page_counter = 0
        else:
            self.page_counter += 1
        self.update_raster_canvas()

    # Move to the previous view on the treatment view tab.
    def decrement_page(self):
        if self.page_counter == 0:
            self.page_counter = 2
        else:
            self.page_counter -= 1
        self.update_raster_canvas()

    # Set the treatment view tab.
    def set_page(self, i):
        if 0 <= i <= 2:
            self.page_counter = i
        else:
            raise ValueError("i is not a valid page index.")
        self.update_raster_canvas()

    def raster_button_clicked(self):
        self.set_page(0)
        self.update_raster_canvas()

    # Parent function for figuring out what to draw on the treat view.
    def update_raster_canvas(self):
        row = self.stand_tab_list_view.currentIndex().row()
        unit = self.rx_units[row]
        if self.page_counter == 0:
            self._draw_raster(unit)
        elif self.page_counter == 1:
            self._draw_report_page(unit)
        else:
            self._draw_cut_report_page(unit)

    # This draws the view that shows the treatment visualization rasters.
    def _draw_raster(self, unit):
        self.raster_widget.setVisible(True)
        self.treat_display.setVisible(False)
        self.cut_trees_display.setVisible(False)

        clicked = self.raster_button.isChecked()
        viewmode = self.raster_viewmode.currentIndex()

        if clicked: # i.e. display treated data
            method = self.treatment_method.currentText()
            kwargs = {'dbh_cutoff': self.cut_range.value()}
            if method == "Add Clumps":
                method = "add_clumps"
                kwargs['prop_table'] = self.project_settings.get_mcs_prop()
            if method == "Sum Squares":
                method = "sum_squares"
            if method == "DBH Thin":
                method = "dbh_thin"

            self.treatment.chm, self.treatment.hill, self.treatment.tao_pts, self.treatment.basin, _ = unit.get_treatment(method, **kwargs)

            clumps = unit.get_treat_clump_map()
            #in some edge cases taos in the tao list might be outside the treated chm...
            # We remove them here to avoid a crash later.
            rem_pts = self.treatment.tao_pts
            rem_pts = rem_pts[~np.logical_or(rem_pts[:, 0] < self.treatment.chm.xmin, rem_pts[:, 0] > self.treatment.chm.xmax), :]
            rem_pts = rem_pts[~np.logical_or(rem_pts[:, 1] < self.treatment.chm.ymin, rem_pts[:, 1] > self.treatment.chm.ymax), :]

            try:
                x = self.treatment.chm.colFromX(rem_pts[:, 0], False)
                y = self.treatment.chm.rowFromY(rem_pts[:, 1], False)
            except BaseException as e:
                print(rem_pts[np.logical_or(rem_pts[:, 0] < self.treatment.chm.xmin, rem_pts[:, 0] > self.treatment.chm.xmax), 0])
                print(rem_pts[np.logical_or(rem_pts[:, 1] < self.treatment.chm.ymin, rem_pts[:, 1] > self.treatment.chm.ymax), 1])
                raise e

            x = x[(0 <= x) & (x < self.treatment.chm.ncol) & (0 <= y) & (y < self.treatment.chm.nrow)]
            y = y[(0 <= x) & (x < self.treatment.chm.ncol) & (0 <= y) & (y < self.treatment.chm.nrow)]
        else:
            self.treatment.chm = unit.get_chm()
            self.treatment.hill = unit.get_hill()
            self.treatment.basin = unit.get_basin_map()
            self.treatment.tao_pts = unit.get_tao_points()
            clumps = unit.get_clump_map()
            x = []
            y = []

        self.raster_figure.clf()
        self.raster_ax = self.raster_figure.add_subplot(111)
        if viewmode == 0:
            img = self.raster_ax.imshow(self.treatment.chm.values, cmap='coolwarm', vmin=0)
            self.raster_ax.imshow(self.treatment.hill, cmap='Greys', alpha=0.5)
            cb_label = "Canopy Height Model (Meters)"
        elif viewmode == 1:
            img = self.raster_ax.imshow(self.treatment.basin.values, cmap='nipy_spectral')
            self.raster_ax.imshow(self.treatment.hill, cmap='Greys', alpha=0.5)
            cb_label = "Basin ID (unique values)"
        else:
            colors = ("white", "#7bc043", "#fdf498", "#f37736", "#ee4035")
            cm_name = "Clump Colors"
            n_bin = 5
            cm = LinearSegmentedColormap.from_list(cm_name, colors, n_bin)
            n_bin_ranges = (-0.5, 0.5, 1.5, 4.5, 9.5, 99)
            norm = BoundaryNorm(n_bin_ranges, len(n_bin_ranges))
            img = self.raster_ax.imshow(clumps.values, cmap=cm, norm=norm)
            self.raster_ax.imshow(self.treatment.hill, cmap='Greys', alpha=0.5)
            cb_label = "Clump Map (Clump bins)"

        # special legend for clump bins, otherwise normal legend.
        if viewmode == 2:
            self.raster_cb = self.raster_figure.colorbar(img, orientation="vertical", label=cb_label,
                                                         ticks=[0, 1, 3, 7, 55])
            self.raster_cb.ax.set_yticklabels(['0', '1', '2-4', '4-9', '10+'])
        else:
            self.raster_cb = self.raster_figure.colorbar(img, orientation="vertical", label=cb_label)

        # Finishing touches on the figure.
        self.raster_ax.scatter(x, y, color="r", s=0.005)
        xrange = self.treatment.chm.ncol
        yrange = self.treatment.chm.nrow
        self.raster_ax.set_xlim(-0.5 - xrange * 0.05, xrange + xrange * 0.05)
        self.raster_ax.set_ylim(yrange + yrange * 0.05, -0.5 - yrange * 0.05)
        if 'f' in unit.get_units():
            if xrange * self.treatment.chm.xres > 10000 or yrange * self.treatment.chm.yres > 10000:
                factor = 1 / 5280
                label = 'miles'
            else:
                factor = 1
                label = 'feet'
        else:
            if xrange * self.treatment.chm.xres > 3048 or yrange * self.treatment.chm.yres > 3048:
                factor = 1 / 1609.34
                label = 'miles'
            else:
                factor = 3.28084
                label = 'feet'

        self.raster_ax.set_xlabel("Distance, %s" % label)
        self.raster_ax.set_ylabel("Distance, %s" % label)

        # Set title and format axes to the units rather than pixels.
        result = unit.treatment_result
        if clicked:
            if result == 1:
                result = "; Failed due to diameter limit."
            else:
                result = ""
        else:
            result = ""
        self.raster_ax.set_title(
            self.project_settings.get_name() + ', ' + str(unit.get_name()) + ' ' + str(self.raster_viewmode.currentText()) + result)
        self.raster_canvas.draw_idle()

        self.raster_ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: round(x * self.treatment.chm.xres * factor)))
        self.raster_ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: round(x * self.treatment.chm.xres * factor)))

        for tick in self.raster_ax.get_xticklabels():
            tick.set_rotation(45)
        self.raster_canvas.draw_idle()

    # Draw the page that gives info on the pre treat and post treat conditions.
    def _draw_report_page(self, unit):
        self.raster_widget.setVisible(False)
        self.treat_display.setVisible(True)
        self.cut_trees_display.setVisible(False)

        self.current_label.setText('Current\n' + str(unit.get_current_structure()))
        if unit.get_treat_structure() is not None:
            self.displayed_label.setText('Post-Treatment\n' + str(unit.get_treat_structure()))
        else:
            self.displayed_label.setText("-")
        self.target_label.setText('Target\n' + str(unit.get_target_structure()))

        ba_label = r"$feet^2\ $"
        tpa_label = r"$ trees\ ac^{-1}$"

        self.displayed_ax_ba.cla()
        self.displayed_ax_ba.set_xlabel("Basal Area (%s)" % ba_label)
        self.displayed_ax_ba.set_ylabel('Kernal Density')
        self.displayed_ax_ba.set_title('Post-Treatment Tree Basal Area Distribution')

        self.displayed_ax_mcs.cla()
        self.displayed_ax_mcs.set_xlabel("Clump Size (n trees)")
        self.displayed_ax_mcs.set_ylabel('Kernal Density')
        self.displayed_ax_mcs.set_title('Post-Treatment Clump Size Distribution')

        self.report_current_ax_ba.cla()
        self.report_current_ax_ba.set_xlabel("Basal Area (%s)" % ba_label)
        self.report_current_ax_ba.set_ylabel('Kernal Density')
        self.report_current_ax_ba.set_title('Pre-Treatment Tree Basal Area Distribution')

        self.report_current_ax_mcs.cla()
        self.report_current_ax_mcs.set_xlabel("Clump Size (n trees)")
        self.report_current_ax_mcs.set_ylabel('Kernal Density')
        self.report_current_ax_mcs.set_title('Pre-Treatment Clump Size Distribution')

        cur_ba = unit.get_ba_dist(unit.get_tao_points())
        cur_ba_dens = gaussian_kde(cur_ba)
        cur_ba_x = np.linspace(0, max(cur_ba), 200)
        cur_ba_dens.covariance_factor = lambda: 0.25
        cur_ba_dens._compute_covariance()
        self.report_current_ax_ba.plot(cur_ba_x, cur_ba_dens(cur_ba_x))

        cur_csd = unit.get_csd()['clump_size']
        cur_mcs_dens = gaussian_kde(cur_csd)
        cur_mcs_x = np.linspace(0, max(cur_csd), 200)
        cur_mcs_dens.covariance_factor = lambda: 0.25
        cur_mcs_dens._compute_covariance()
        self.report_current_ax_mcs.plot(cur_mcs_x, cur_mcs_dens(cur_mcs_x))

        self.displayed_mcs_prop.setText('')

        # If there is a stored treatment.
        if unit.get_treat_points() is not None:
            disp_ba = unit.get_ba_dist(unit.get_treat_points())
            disp_ba_dens = gaussian_kde(disp_ba)
            disp_ba_x = np.linspace(0, max(disp_ba), 200)
            disp_ba_dens.covariance_factor = lambda: 0.25
            disp_ba_dens._compute_covariance()
            self.displayed_ax_ba.plot(disp_ba_x, disp_ba_dens(disp_ba_x))

            disp_csd = unit._get_csd(unit.get_treat_points())['clump_size']
            disp_mcs_dens = gaussian_kde(disp_csd)
            disp_mcs_x = np.linspace(0, max(disp_csd), 200)
            disp_mcs_dens.covarariance_factor = lambda: 0.25
            disp_mcs_dens._compute_covariance()
            self.displayed_ax_mcs.plot(disp_mcs_x, disp_mcs_dens(disp_mcs_x))

            ba_xrange = (min(self.report_current_ax_ba.get_xlim()[0], self.displayed_ax_ba.get_xlim()[0]),
                         max(self.report_current_ax_ba.get_xlim()[1], self.displayed_ax_ba.get_xlim()[1]))
            ba_yrange = (min(self.report_current_ax_ba.get_ylim()[0], self.displayed_ax_ba.get_ylim()[0]),
                         max(self.report_current_ax_ba.get_ylim()[1], self.displayed_ax_ba.get_ylim()[1]))

            mcs_xrange = (min(self.report_current_ax_mcs.get_xlim()[0], self.displayed_ax_mcs.get_xlim()[0]),
                          max(self.report_current_ax_mcs.get_xlim()[1], self.displayed_ax_mcs.get_xlim()[1]))
            mcs_yrange = (min(self.report_current_ax_mcs.get_ylim()[0], self.displayed_ax_mcs.get_ylim()[0]),
                          max(self.report_current_ax_mcs.get_ylim()[1], self.displayed_ax_mcs.get_ylim()[1]))

            self.displayed_ax_ba.set_xlim(ba_xrange)
            self.displayed_ax_ba.set_ylim(ba_yrange)
            self.displayed_ax_mcs.set_xlim(mcs_xrange)
            self.displayed_ax_mcs.set_ylim(mcs_yrange)
            self.report_current_ax_mcs.set_xlim(mcs_xrange)
            self.report_current_ax_mcs.set_ylim(mcs_yrange)
            self.report_current_ax_ba.set_xlim(ba_xrange)
            self.report_current_ax_ba.set_ylim(ba_yrange)

            mcs_table = self.project_settings.get_mcs_prop()
            mcs_prop = mcs_table[np.logical_and(mcs_table['MCS_min'] < unit.get_treat_structure().mcs,
                                                mcs_table['MCS_max'] >= unit.get_treat_structure().mcs)][0]
            self.displayed_mcs_prop.setText("Clumps:  1\t2-4\t5-9\t10-14\t15-29\t30+\n"
                                            "Props:    " + str(mcs_prop[2]) + '\t' + str(mcs_prop[3]) + '\t' + str(mcs_prop[4]) +
                                            '\t' + str(mcs_prop[5]) + '\t' + str(mcs_prop[6]) + '\t' + str(mcs_prop[7]))

    # Draw the page that shows the info on the cut trees.
    def _draw_cut_report_page(self, unit):
        self.raster_widget.setVisible(False)
        self.treat_display.setVisible(False)
        self.cut_trees_display.setVisible(True)

        self.cut_trees_label.setText('Cut Trees:\nBA:')
        ba_label = r"$feet^2\ $"

        self.cut_trees_ax.cla()
        self.cut_trees_ax.set_xlabel("Basal Area (%s)" % ba_label)
        self.cut_trees_ax.set_ylabel('Density')
        self.cut_trees_ax.set_title('Cut Trees Basal Area Distribution')

        pts = unit.cut_taos
        if pts is not None:
            dbh = unit._tao_data.get_dbh_from_height_dll(pts[:, 3])
            ba = 0.005454 * (dbh / 2.54)**2

            area = unit.get_unit_polygon().area
            if 'f' in unit.get_units():
                area /= 43560
            else:
                area /= 10000
                area *= 2.47105
            ba_ac = ba / area

            cut_ba_dens = gaussian_kde(ba)
            cut_ba_x = np.linspace(0, max(ba), 200)
            cut_ba_dens.covariance_factor = lambda: 0.25
            cut_ba_dens._compute_covariance()
            self.cut_trees_ax.plot(cut_ba_x, cut_ba_dens(cut_ba_x))
            self.cut_trees_label.setText('Cut Trees:\n'
                                         + '\nBA :\t' + str(np.round(ba_ac.sum(), 3)))

    # Function to be connected to a signal to reset the treatment raster button.
    def reset_button(self):
        self.raster_button.setChecked(False)

    # Function to connect to a signal to reset the page counter.
    def reset_page(self):
        self.set_page(0)

    # Update the hover annotation on the overall view
    @staticmethod
    def update_annot(ind, names, child, annot):
        x, y = child.get_data()
        annot.xy = (x[ind['ind'][0]], y[ind['ind'][0]])
        #annot.xy = (50, 50)
        text = '{}'.format("\n".join([str(names[n]) for n in ind['ind']]))
        annot.set_text(text)

    # This updates the hover targ and treated points and arrows.
    def hover_arrows(self, unit_idx, ax_idx, metric_idx):
        axes = [self.ba_ax, self.mcs_ax, self.cc_ax]

        current = self.rx_units[unit_idx].get_current_structure()
        target = self.rx_units[unit_idx].get_target_structure()
        self.target_points[ax_idx].set_data([target.tpa], [target[metric_idx]])
        self.target_points[ax_idx].set_visible(True)
        self.cur_targ_arrow[ax_idx].set_data(x=current.tpa, y=current[metric_idx], dx=target.tpa-current.tpa,
                                       dy=target[metric_idx] - current[metric_idx])
        self.cur_targ_arrow[ax_idx].set_visible(True)

        treated = self.rx_units[unit_idx].get_treated_structure()
        if treated is not None:
            self.treated_points[ax_idx].set_data([treated.tpa], [treated[metric_idx]])
            self.treated_points[ax_idx].set_visible(True)
            self.targ_treated_arrow[ax_idx].set_data(x=target.tpa, y=target[metric_idx], dx=treated.tpa - target.tpa,
                                           dy=treated[metric_idx] - target[metric_idx])
            self.targ_treated_arrow[ax_idx].set_visible(True)

    # This resets the arrows and points.
    def clean_arrows(self, ax_idx):
        self.target_points[ax_idx].set_visible(False)
        self.treated_points[ax_idx].set_visible(False)
        self.cur_targ_arrow[ax_idx].set_visible(False)
        self.targ_treated_arrow[ax_idx].set_visible(False)

    # hover prep for current ba.
    def ba_hover(self, event):
        #clean slate
        self.ba_annot.set_visible(False)
        for i in range(len(self.unit_ba_polygons)):
            self.unit_ba_polygons[i].set_visible(False)
            self.unit_ba_intersections[i].set_visible(False)
        self.clean_arrows(0)

        #check what to show
        if event.inaxes == self.ba_ax:
            cont_2, ind_2 = self.unit_ba_points.contains(event)
            if cont_2:
                self.update_annot(ind_2, self.unit_names, self.unit_ba_points, self.ba_annot)
                x = ind_2['ind'][0]
                self.unit_ba_polygons[x].set_visible(True)
                self.unit_ba_intersections[x].set_visible(True)
                self.ba_annot.set_visible(True)
                self.hover_arrows(x, 0, 1)

        self.reference_canvas.draw_idle()

    def mcs_hover(self, event):
        # clean slate
        self.mcs_annot.set_visible(False)
        for i in range(len(self.unit_ba_polygons)):
            self.unit_mcs_polygons[i].set_visible(False)
            self.unit_mcs_intersections[i].set_visible(False)
        self.clean_arrows(1)

        # check what to show
        if event.inaxes == self.mcs_ax:
            cont_2, ind_2 = self.unit_mcs_points.contains(event)
            if cont_2:
                self.update_annot(ind_2, self.unit_names, self.unit_mcs_points, self.mcs_annot)
                x = ind_2['ind'][0]
                self.unit_mcs_polygons[x].set_visible(True)
                self.unit_mcs_intersections[x].set_visible(True)
                self.mcs_annot.set_visible(True)
                self.hover_arrows(x, 1, 2)
        self.reference_canvas.draw_idle()

    def cc_hover(self, event):
        # clean slate
        self.cc_annot.set_visible(False)
        for i in range(len(self.unit_cc_polygons)):
            self.unit_cc_polygons[i].set_visible(False)
            self.unit_cc_intersections[i].set_visible(False)
        self.clean_arrows(2)

        # check what to show
        if event.inaxes == self.cc_ax:
            cont_2, ind_2 = self.unit_cc_points.contains(event)
            if cont_2:
                self.update_annot(ind_2, self.unit_names, self.unit_cc_points, self.cc_annot)
                x = ind_2['ind'][0]
                self.unit_cc_polygons[x].set_visible(True)
                self.unit_cc_intersections[x].set_visible(True)
                self.cc_annot.set_visible(True)
                self.hover_arrows(x, 2, 4)
        self.reference_canvas.draw_idle()

    # Function to connect to a signal from clicking on a unit.
    def rx_unit_pick(self, event):
        if event.button == 1:
            cont_ba, ind_ba = self.unit_ba_points.contains(event)
            cont_mcs, ind_mcs = self.unit_mcs_points.contains(event)
            cont_cc, ind_cc = self.unit_cc_points.contains(event)
            if cont_ba or cont_mcs or cont_cc:
                self.tabs.setCurrentIndex(0)
                if cont_ba:
                    ind = QModelIndex(self.model.index(ind_ba['ind'][0], 0))
                elif cont_mcs:
                    ind = QModelIndex(self.model.index(ind_mcs['ind'][0], 0))
                else:
                    ind = QModelIndex(self.model.index(ind_cc['ind'][0], 0))
                self.stand_tab_list_view.setCurrentIndex(ind)

    # This exports the current "raster canvas" whether it's the treatment view, or one of the report pages.
    def export_tif(self):
        if self.tabs.currentIndex() == 0:
            if self.page_counter == 0:
                row = self.stand_tab_list_view.currentIndex().row()
                unit = self.rx_units[row]
                filename = QFileDialog.getSaveFileName(self, "Enter a file name to export data to.", "", '*.tif')[0]
                with tempfile.TemporaryDirectory() as tmpdir:
                    self.raster_figure.savefig(tmpdir + "/tmp.tif", dpi=300)
                    img = PIL.Image.open(tmpdir + "/tmp.tif")
                    i = PIL.ImageDraw.Draw(img)
                    f = PIL.ImageFont.truetype("arial.ttf", 70)
                    i.text((100, 100), 'Current\n' + str(unit.get_current_structure()), font=f, fill=(0, 0, 0))
                    if unit.get_treat_structure() is not None:
                        i.text((100, 800), 'Post-Treatment\n' + str(unit.get_treat_structure()), font=f, fill=(0, 0, 0))
                    else:
                        i.text((100, 800), "-", font=f, fill=(0, 0, 0))
                    i.text((100, 1600), 'Target\n' + str(unit.get_target_structure()), font=f, fill=(0, 0, 0))
                    img.save(filename)
                    img.close()

                fig_dpi = self.raster_figure.dpi
                # Corners array:
                #   [[left, bottom]
                #    [right, top]]
                corners = self.raster_ax.bbox.get_points() * 300 / fig_dpi

                r = self.get_current_unit().get_chm()
                proj = r.projection

                e = r.extent
                xrange = e.xmax - e.xmin
                yrange = e.ymax - e.ymin

                e.xmax = e.xmax + xrange * 0.05
                e.xmin = e.xmin - xrange * 0.05
                e.ymin = e.ymin - yrange * 0.05
                e.ymax = e.ymax + yrange * 0.05

                xres = (e.xmax - e.xmin) / (corners[1, 0] - corners[0, 0])  # Positive  value, left side = 0
                yres = (e.ymax - e.ymin) / (corners[0, 1] - corners[1, 1])*-1  # Negative value, top = 0

                width, height = self.raster_canvas.get_width_height()
                height *= 300 / fig_dpi
                xmin = e.xmin - corners[0, 0] * xres
                ymax = e.ymax + (height - corners[1, 1]) * yres

                ds = rasterio.open(filename)
                data = ds.read()
                ht = ds.height
                wd = ds.width
                ext = raster.Extent(xmin, xmin+wd*xres, ymax-ht*yres, ymax)
                r = raster.Raster(ext, (xres, yres), wd, ht, proj, data, ds.dtypes[0])
                ds.close()
                r.writeRaster(filename)

            elif self.page_counter == 1:
                filename = QFileDialog.getSaveFileName(self, "Enter a file name to export data to.", "", '*.tif')[0]
                x = QApplication.primaryScreen()
                p = x.grabWindow(self.treat_display.winId())
                p.save(filename, 'tif')
            else:
                filename = QFileDialog.getSaveFileName(self, "Enter a file name to export data to.", "", '*.tif')[0]
                x = QApplication.primaryScreen()
                p = x.grabWindow(self.cut_trees_display.winId())
                p.save(filename, 'tif')
        elif self.tabs.currentIndex() == 1:
            filename = QFileDialog.getSaveFileName(self, "Enter a file name to export data to.", "", '*.tif')[0]
            x = QApplication.primaryScreen()
            p = x.grabWindow(self.land_tab.winId())
            p.save(filename, 'tif')



    # export the current raster depending on whether the treatment view is check or not.
    def export_raster(self):
        filename = QFileDialog.getSaveFileName(self, "Enter a file name to export raster data to.", "", "'tif")[0]
        if filename[len(filename) - 4:] != ".tif":
            filename = filename + ".tif"
        unit = self.get_current_unit()
        viewmode = self.raster_viewmode.currentIndex()
        if viewmode == 0:
            if self.raster_button.isChecked():
                r = unit.get_treat_chm()
            else:
                r = unit.get_chm()
        elif viewmode == 1:
            if self.raster_button.isChecked():
                r = unit._treat_basin
            else:
                r = unit._basin_map
        else:
            if self.raster_button.isChecked():
                r = unit.get_treat_clump_map()
            else:
                r = unit.get_clump_map()
        r.writeRaster(filename)

    # Export the tree points
    def export_features(self):
        filename = QFileDialog.getSaveFileName(self, "Enter a file name to export point data to.", "", "'.shp")[0]
        if filename[len(filename) - 4:] != ".shp":
            filename = filename + ".shp"
        unit = self.get_current_unit()
        if self.raster_button.isChecked():
            data = unit.get_treat_points()
        else:
            data = unit.get_tao_points()
        pts = [Point(row) for row in data[:, (0, 1, 3)]]

        schema = {
            'geometry': 'Point',
            'properties': {'x': 'float',
                           'y': 'float',
                           'area': 'float',
                           'height': 'float',
                           'crown': 'float',
                           }
        }

        wkt = unit.get_wkt()
        crs = fiona.crs.from_string(wkt)

        with fiona.open(filename, 'w', crs=crs, driver='ESRI Shapefile', schema=schema) as c:
            for i, pt in enumerate(pts):
                c.write({
                    'geometry': mapping(pt),
                    'properties': {'x': data[i, 0],
                                   'y': data[i, 1],
                                   'area': data[i, 2],
                                   'height': data[i, 3],
                                   'crown': data[i, 4]}
                })

    # export csv treelist
    def export_treelist(self):
        filename = QFileDialog.getSaveFileName(self, "Enter a file name to export tree list data to.", "", "'.csv")[0]
        if filename[len(filename) - 4:] != ".csv":
            filename = filename + ".csv"
        headers = ['x', 'y', 'area', 'height', 'crown']
        unit = self.get_current_unit()
        if self.raster_button.isChecked():
            data = unit.get_treat_points()
        else:
            data = unit.get_tao_points()
        np.savetxt(filename, data, delimiter=',')
        with open(filename, 'r') as f:
            r = list(csv.reader(f))
            r.insert(0, headers)
        with open(filename, 'w', newline='') as out:
            w = csv.writer(out)
            for line in r:
                w.writerow(line)

    # Not really working yet, hence it is not shown in the file menu yet.
    def print(self):
        printer = QtPrintSupport.QPrinter()
        printer.setPageSize(QtGui.QPagedPaintDevice.Letter)
        printer.setColorMode(QtPrintSupport.QPrinter.ColorMode.Color)
        printer.setOutputFormat(QtPrintSupport.QPrinter.PdfFormat)
        printer.setOutputFileName("E:/test.pdf")
        painter = QtGui.QPainter(printer)
        #painter.begin(printer)

        screen = self.raster_canvas.grab().scaled(
            printer.pageRect(QtPrintSupport.QPrinter.DevicePixel).size().toSize(),
            Qt.KeepAspectRatio)
        painter.drawPixmap(0, 0, screen)

        screen = self.reference_canvas.grab().scaled(
            printer.pageRect(QtPrintSupport.QPrinter.DevicePixel).size().toSize(),
            Qt.KeepAspectRatio)
        painter.drawPixmap(0, 500, screen)

        painter.end()

    # Helper function to return the unit currently displayed.
    def get_current_unit(self):
        return self.rx_units[self.stand_tab_list_view.currentIndex().row()]

# Table model to relate a list view gui with the target seting widget and the raster views.
class RxUnitTableModel(QAbstractTableModel):
    def __init__(self, rxunits=[], decision_space = [], parent=None):
        QAbstractTableModel.__init__(self, parent)

        # Get all the relevent data out of the units.
        rxunits_names = [rx.get_name() for rx in rxunits]
        rxunits_current_tpa = [rx.get_current_structure().tpa for rx in rxunits]
        rxunits_current_ba = [rx.get_current_structure().ba for rx in rxunits]
        rxunits_current_mcs = [rx.get_current_structure().mcs for rx in rxunits]
        rxunits_current_cc = [rx.get_current_structure().cc for rx in rxunits]
        rxunits_target_tpa = [rx.get_target_structure().tpa for rx in rxunits]
        rxunits_target_ba = [rx.get_target_structure().ba for rx in rxunits]
        rxunits_target_mcs = [rx.get_target_structure().mcs for rx in rxunits]
        rxunits_target_cc = [rx.get_target_structure().cc for rx in rxunits]

        self._rxunits = rxunits
        # self._decision_space = decision_space
        # self._locks = [[False, False, False, False] for _ in self._rxunits]

        # zip names and data twice, once for currents, and once as a temp for target.
        self._data = [[a, b, c, d, e, f, g, h, i] for a, b, c, d, e, f, g, h, i in zip(rxunits_names,
                       rxunits_current_tpa, rxunits_current_ba, rxunits_current_mcs, rxunits_current_cc,
                       rxunits_target_tpa, rxunits_target_ba, rxunits_target_mcs, rxunits_target_cc)]

    def rowCount(self, parent=QModelIndex()):
        return len(self._data)

    def columnCount(self, parent=QModelIndex()):
        return len(self._data[0])

    # This interfaces the datamodel and the querier
    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            row = index.row()
            col = index.column()
            dat = self._data[row][col]
            if role == Qt.DisplayRole:
                if type(dat) is str:
                    return dat
                else:
                    return str(round(dat, 3))
            if role == Qt.EditRole:
                if type(dat) is str:
                    return dat
                else:
                    return str(round(dat, 3))
            if role == Qt.ToolTipRole:
                return str(self._rxunits[row].get_current_structure())
        return QVariant()

    def flags(self, index):
        flags = super(RxUnitTableModel, self).flags(index)
        return flags

    # This updates the data in the model and the overall rx list
    def setData(self, index, value, role=Qt.EditRole):
        try:
            value = float(value)
            if role == Qt.EditRole and index.isValid():
                row = index.row()
                col = index.column()
                self._data[row][col] = value
                c = col - 5
                if c == 3: #skip OSI in the structure summary class.
                    c = 4
                self._rxunits[row].set_target_structure_by_idx(c, value)
                self.dataChanged.emit(index, index)
                return True
        except Exception as e:
            print(e)
            pass
        self.dataChanged.emit(index, index)
        return False


# This class is the widget that surfaces the target and current structure from the datamodel.
class StructureInfo(QWidget):
    def __init__(self, parent=None):
        super(StructureInfo, self).__init__(parent)
        self.mapper = None

        self.name = QLabel()
        self.current = QLabel("Current")
        self.target = QLabel("Targets")
        self.tpa = QLabel("TPA")
        self.ba = QLabel("BA")
        self.mcs = QLabel("MCS")
        self.cc = QLabel("CC")

        self.current_tpa = QLabel()
        self.current_ba = QLabel()
        self.current_mcs = QLabel()
        self.current_cc = QLabel()

        self.target_tpa = QLineEdit()
        self.target_ba = QLineEdit()
        self.target_mcs = QLineEdit()
        self.target_cc = QLineEdit()

        layout = QGridLayout()
        layout.addWidget(self.name, 0, 2, 1, 2)
        layout.addWidget(self.current, 1, 2)
        layout.addWidget(self.target, 1, 3)
        layout.addWidget(self.tpa, 2, 1)
        layout.addWidget(self.ba, 3, 1)
        layout.addWidget(self.mcs, 4, 1)
        layout.addWidget(self.cc, 5, 1)
        layout.addWidget(self.current_tpa, 2, 2)
        layout.addWidget(self.current_ba, 3, 2)
        layout.addWidget(self.current_mcs, 4, 2)
        layout.addWidget(self.current_cc, 5, 2)
        layout.addWidget(self.target_tpa, 2, 3)
        layout.addWidget(self.target_ba, 3, 3)
        layout.addWidget(self.target_mcs, 4, 3)
        layout.addWidget(self.target_cc, 5, 3)

        self.setLayout(layout)

    # Set mapping for the widget that reports the treatment conditions.
    def set_model(self, model):
        self.mapper = QDataWidgetMapper(self)

        self.mapper.setModel(model)
        self.mapper.addMapping(self.name, 0, b"text")
        self.mapper.addMapping(self.current_tpa, 1, b"text")
        self.mapper.addMapping(self.current_ba, 2, b"text")
        self.mapper.addMapping(self.current_mcs, 3, b"text")
        self.mapper.addMapping(self.current_cc, 4, b"text")
        self.mapper.addMapping(self.target_tpa, 5)
        self.mapper.addMapping(self.target_ba, 6)
        self.mapper.addMapping(self.target_mcs, 7)
        self.mapper.addMapping(self.target_cc, 8)
        self.mapper.toFirst()

    def set_selection(self, current, old):
        self.mapper.setCurrentModelIndex(current)


# Dummy class to store the treatment info (for multithreading in theory).
class Treatment:
    def __init__(self, chm, hill, basin, tao_pts):
        self.chm = chm
        self.hill = hill
        self.basin = basin
        self.tao_pts = tao_pts
