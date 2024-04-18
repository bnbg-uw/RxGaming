# General
import matplotlib.pyplot as plt
import numpy as np
import os

# QT
from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton, QTabWidget, QWidget, QLabel, QListView,
                             QDataWidgetMapper, QLineEdit, QComboBox, QFileDialog, QSizePolicy, QAction, QApplication)
from PyQt5.QtCore import QAbstractTableModel, Qt, QVariant, QModelIndex

# RxGaming
from activity import Activity, SaveStateActivity
from projectsettings import DllStorage
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
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas  # Allows us to redner plots in widgets.
import descartes # For polygons on ref plots
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Polygon as MplPolygon
import matplotlib.ticker as ticker  # custom axis ticks
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm  # coloring the displayed raster
from scipy.spatial import ConvexHull  # polygons on reference plots
from scipy.stats import gaussian_kde  # report pages
import seaborn


# GamingActivity class
# Handles user input and interactions with the projectsettings and surfaces data and treatment results
class GamingActivity(Activity):
    def on_start(self, saved_state, **kwargs):
        self.dll_path = kwargs['dll_path']
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

        # Set ups the menu bar
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

        self.main_menu = self.window.menuBar()
        self.file_menu = self.main_menu.addMenu("File")
        self.export_menu = self.file_menu.addMenu("Export")

        self.file_menu.addAction(self.save_action)
        self.file_menu.addAction(self.save_as_action)
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
        saved_state['LastActivity'] = type(self)
        if 'save_file_location' in self.saved_state:
            saved_state['save_file_location'] = self.saved_state['save_file_location']
        return saved_state

    # action for a user instantiated save.
    def menu_save(self):
        if 'save_file_location' in self.saved_state:
            print(self.saved_state['save_file_location'])
            if os.path.isfile(self.saved_state['save_file_location']):
                with open(self.saved_state['save_file_location'], 'wb') as fp:
                    pickle.dump(self.save(), fp)
                self.tab_widget.project_settings.prj_area.dePickle(self.dll_path)
            else:
                self.menu_save_as()
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

class Tabs(QTabWidget):
    def __init__(self, saved_state, dll_path):
        super(Tabs, self).__init__(None)
        self.layout = QVBoxLayout()

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
        self.report_label = QLabel("Treatment Report")
        self.current_label = QLabel("Current")
        self.displayed_label = QLabel("Post-Treatment")
        self.target_label = QLabel("Target")
        self.report_label.setStyleSheet("font: 24pt;")
        self.current_label.setStyleSheet('font: 16pt;')
        self.displayed_label.setStyleSheet('font: 16pt;')
        self.target_label.setStyleSheet('font: 16pt;')
        self.current_fig_ba = Figure()
        self.current_canvas_ba = FigureCanvas(self.current_fig_ba)
        self.report_current_ax_ba = self.current_fig_ba.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.current_fig_mcs = Figure()
        self.current_canvas_mcs = FigureCanvas(self.current_fig_mcs)
        self.report_current_ax_mcs = self.current_fig_mcs.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.displayed_fig_ba = Figure()
        self.displayed_canvas_ba = FigureCanvas(self.displayed_fig_ba)
        self.displayed_ax_ba = self.displayed_fig_ba.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.displayed_fig_mcs = Figure()
        self.displayed_canvas_mcs = FigureCanvas(self.displayed_fig_mcs)
        self.displayed_ax_mcs = self.displayed_fig_mcs.add_subplot(111, position=[0.15, 0.15, 0.75, 0.75])
        self.displayed_mcs_prop = QLabel("")
        self.displayed_mcs_prop.setStyleSheet('font: 16pt;')

        # Set up treatment window layouts and widgets.
        self.treat_display = QWidget()
        self.treat_display.setVisible(False)
        self.treat_display_layout = QGridLayout()
        self.treat_display.setLayout(self.treat_display_layout)

        self.treat_display_layout.addWidget(self.report_label, 0, 1)
        self.treat_display_layout.addWidget(self.current_label, 1, 0)
        self.treat_display_layout.addWidget(self.displayed_label, 2, 0)
        self.treat_display_layout.addWidget(self.target_label, 3, 0)
        self.treat_display_layout.addWidget(self.current_canvas_ba, 1, 1)
        self.treat_display_layout.addWidget(self.displayed_canvas_ba, 2, 1)
        self.treat_display_layout.addWidget(self.displayed_mcs_prop, 3, 2)
        self.treat_display_layout.addWidget(self.current_canvas_mcs, 1, 2)
        self.treat_display_layout.addWidget(self.displayed_canvas_mcs, 2, 2)

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
        self.raster_button.clicked[bool].connect(self.update_raster_canvas)
        self.raster_button.setStyleSheet("QPushButton:enabled {color: black;}\n\
                                         QPushButton:checked { background-color: rgb(80, 80, 80); \
                                         border: none;\
                                         color: white}")

        self.cut_range = SliderWithValue(Qt.Horizontal)
        self.cut_range.setMinimum(0)
        self.cut_range.setMaximum(120)
        self.cut_range.setValue(21)

        self.raster_viewmode = QComboBox()
        self.raster_viewmode.addItems(["Canopy Model", "Basins", "Clumps"])
        self.raster_viewmode.currentIndexChanged.connect(self.update_raster_canvas)

        self.treatment_method = QComboBox()
        self.treatment_method.addItems(["Add Clumps", 'DBH Thin', "Sum Squares"])

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
        self.upper_figure = Figure()
        self.current_ba_ax = self.upper_figure.add_subplot(1, 3, 1)
        self.current_mcs_ax = self.upper_figure.add_subplot(1, 3, 2)
        self.current_cc_ax = self.upper_figure.add_subplot(1, 3, 3)
        self.upper_canvas = FigureCanvas(self.upper_figure)

        self.lower_figure = Figure()
        self.lower_ba_ax = self.lower_figure.add_subplot(131)
        self.lower_mcs_ax = self.lower_figure.add_subplot(132)
        self.lower_cc_ax = self.lower_figure.add_subplot(133)
        self.lower_canvas = FigureCanvas(self.lower_figure)

        # QT window layout and stuff
        self.land_tab.layout.addWidget(self.upper_canvas)
        self.land_tab.layout.addWidget(self.lower_canvas)
        self.land_tab.setLayout(self.land_tab.layout)

        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

        # Trigger loading and prepping data from save state, and prepare for threading.
        self.load_and_prep_data(saved_state, dll_path)

        ba_label = r"$feet^2\ ac^{-1}$"
        tpa_label = r"$ trees\ ac^{-1}$"

        self.current_ba_ax.set_title("Current Basal Area")
        self.current_ba_ax.set_ylabel("Basal Area (%s)" % ba_label)
        self.current_ba_ax.set_xlabel("Density (%s)" % tpa_label)

        self.current_mcs_ax.set_title("Current Mean Clump Size")
        self.current_mcs_ax.set_ylabel("Mean Clump Size (n trees)")
        self.current_mcs_ax.set_xlabel("Density (%s)" % tpa_label)

        self.current_cc_ax.set_title("Current Canopy Cover")
        self.current_cc_ax.set_ylabel("Canopy Cover (%)")
        self.current_cc_ax.set_xlabel("Density (Trees %s)" % tpa_label)

        #set model for stand tab list view
        self.stand_tab_list_view.setModel(self.model)
        self.stand_tab_info_widget.set_model(self.model)
        self.stand_tab_list_view.selectionModel().currentChanged.connect(self.stand_tab_info_widget.set_selection)
        self.stand_tab_list_view.selectionModel().currentChanged.connect(self.reset_page)
        self.stand_tab_list_view.selectionModel().currentChanged.connect(self.reset_button)
        self.stand_tab_list_view.selectionModel().currentChanged.connect(self.update_raster_canvas)
        self.model.dataChanged.connect(self.update_lower_canvas)

        # set up plots for land tab
        seaborn.kdeplot(ax=self.current_ba_ax, x=self.ref_tpa[self.ref_ind['ba']], y=self.ref_ba[self.ref_ind['ba']],
                    cmap="Oranges", fill=True, bw_adjust=0.5)
        self.unit_ba_points, = self.current_ba_ax.plot(self.current_tpa, self.current_ba, 'b^')

        seaborn.kdeplot(ax=self.current_mcs_ax, x=self.ref_tpa[self.ref_ind['mcs']], y=self.ref_mcs[self.ref_ind['mcs']],
                        cmap="Oranges", fill=True, bw_adjust=0.5)
        self.unit_mcs_points, = self.current_mcs_ax.plot(self.current_tpa, self.current_mcs, 'b^')

        seaborn.kdeplot(ax=self.current_cc_ax, x=self.ref_tpa[self.ref_ind['cc']], y=self.ref_cc[self.ref_ind['cc']],
                        cmap="Oranges", fill=True, bw_adjust=0.5)
        self.unit_cc_points, = self.current_cc_ax.plot(self.current_tpa, self.current_cc, 'b^')

        self.upper_canvas.draw_idle()

        #ANNOTATIONS upper/current:
        self.current_ba_annot = self.current_ba_ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                                                            bbox=dict(boxstyle="round", fc='w'),
                                                            arrowprops=dict(arrowstyle='->'))
        self.current_ba_annot.set_visible(True)
        self.upper_canvas.mpl_connect("motion_notify_event", self.current_ba_hover)

        self.current_mcs_annot = self.current_mcs_ax.annotate("", xy=(0, 0), xytext=(-20, 20),
                                                              textcoords="offset points",
                                                              bbox=dict(boxstyle="round", fc='w'),
                                                              arrowprops=dict(arrowstyle='->'))
        self.current_mcs_annot.set_visible(False)
        self.upper_canvas.mpl_connect("motion_notify_event", self.current_mcs_hover)

        self.current_cc_annot = self.current_cc_ax.annotate("", xy=(0, 0), xytext=(-20, 20),
                                                              textcoords="offset points",
                                                              bbox=dict(boxstyle="round", fc='w'),
                                                              arrowprops=dict(arrowstyle='->'))
        self.current_cc_annot.set_visible(False)
        self.upper_canvas.mpl_connect("motion_notify_event", self.current_cc_hover)

        self.upper_canvas.mpl_connect("button_press_event", self.rx_unit_pick)

        #ANNOTATIONS lower/targets
        self.lower_cc_annot = self.lower_cc_ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                                                          bbox=dict(boxstyle="round", fc='w'),
                                                          arrowprops=dict(arrowstyle='->'))
        self.lower_cc_annot.set_visible(False)
        self.lower_mcs_annot = self.lower_mcs_ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                                                          bbox=dict(boxstyle="round", fc='w'),
                                                          arrowprops=dict(arrowstyle='->'))
        self.lower_mcs_annot.set_visible(False)
        self.lower_ba_annot = self.lower_ba_ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                                                        bbox=dict(boxstyle="round", fc='w'),
                                                        arrowprops=dict(arrowstyle='->'))
        self.lower_ba_annot.set_visible(True)

        self.lower_canvas.mpl_connect("motion_notify_event", self.lower_ba_hover)
        self.lower_canvas.mpl_connect("motion_notify_event", self.lower_mcs_hover)
        self.lower_canvas.mpl_connect("motion_notify_event", self.lower_cc_hover)

        handles = [Patch(facecolor=plt.cm.Oranges(100)),
                   self.unit_cc_points]
        self.upper_figure.legend(handles, ("Units", "Reference"), loc='right')

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
        self.upper_canvas.draw_idle()
        self.update_lower_canvas()
        self.update_raster_canvas()

    # TODO shift off of project_settings.prj_area.get_units() to self.rx_units.
    def load_and_prep_data(self, saved_state, dll_path):
        # Get ref data ready
        # This could be a record array instead of individual arrays??
        self.project_settings = saved_state["ProjectSettings"]
        if(type(self.project_settings.prj_area._tao_data.dll) is DllStorage):
            print("depickling")
            self.project_settings.prj_area.dePickle(dll_path)

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
        self.ref_tpa = ref_db['tph'].astype(float) / 2.47105
        self.ref_ba = ref_db['ba'].astype(float) * 4.356
        self.ref_mcs = ref_db["mcs"].astype(float)
        self.ref_cc = ref_db['cc'].astype(float)
        self.ref_names = ["_".join([str(tp), str(source), str(plot_id)]) for tp, source, plot_id in zip(ref_db['name'],
                                                                                         ref_db['id'],
                                                                                         ref_db['type'])]
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
        self.ref_ba_polygon = descartes.PolygonPatch(ref_ba_shp_poly)

        # MCS
        ref_mcs_data = np.array([self.ref_tpa[self.ref_ind['mcs']], self.ref_mcs[self.ref_ind['mcs']]]).transpose()
        ref_mcs_chull = ConvexHull(ref_mcs_data)
        ref_mcs_shp_poly = ShpPolygon(list(zip(ref_mcs_data[ref_mcs_chull.vertices, 0],
                                               ref_mcs_data[ref_mcs_chull.vertices, 1])))
        self.ref_mcs_polygon = descartes.PolygonPatch(ref_mcs_shp_poly)

        # CC
        ref_cc_data = np.array([self.ref_tpa[self.ref_ind['cc']], self.ref_cc[self.ref_ind['cc']]]).transpose()
        ref_cc_chull = ConvexHull(ref_cc_data)
        ref_cc_shp_poly = ShpPolygon(list(zip(ref_cc_data[ref_cc_chull.vertices, 0],
                                               ref_cc_data[ref_cc_chull.vertices, 1])))
        self.ref_cc_polygon = descartes.PolygonPatch(ref_cc_shp_poly)

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
                rx_unit._tao_data.exportRxToDll(rx_unit)
                x = rx_unit.get_simulated_structures_dll()
                self.decision_spaces.append(x)

        # Get the actionable space for each RxUnit (the intersection of ref conditions and possible treatments)
        #       This may be really inefficient, patch collections??? but I don't know how to only plot one polygon from
        #       a patch collection.
        # BA
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


            if poly.wkt == 'GEOMETRYCOLLECTION EMPTY' or poly.wkt == "POLYGON EMPTY":
                poly = MplPolygon([[0, 0]])
            else:
                poly = descartes.PolygonPatch(poly)
            if inter.wkt == 'GEOMETRYCOLLECTION EMPTY' or inter.wkt == "POLYGON EMPTY":
                inter = MplPolygon([[0, 0]])
            else:
                inter = descartes.PolygonPatch(inter)
            self.current_ba_ax.add_patch(poly)
            self.current_ba_ax.add_patch(inter)
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
            if poly.wkt == 'GEOMETRYCOLLECTION EMPTY' or poly.wkt == "POLYGON EMPTY":
                poly = MplPolygon([[0, 0]])
            else:
                poly = descartes.PolygonPatch(poly)
            if inter.wkt == 'GEOMETRYCOLLECTION EMPTY' or inter.wkt == "POLYGON EMPTY":
                inter = MplPolygon([[0, 0]])
            else:
                inter = descartes.PolygonPatch(inter)
            self.current_mcs_ax.add_patch(poly)
            self.current_mcs_ax.add_patch(inter)
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
            if poly.wkt == 'GEOMETRYCOLLECTION EMPTY' or poly.wkt == "POLYGON EMPTY":
                poly = MplPolygon([[0, 0]])
            else:
                poly = descartes.PolygonPatch(poly)
            if inter.wkt == 'GEOMETRYCOLLECTION EMPTY' or inter.wkt == "POLYGON EMPTY":
                inter = MplPolygon([[0, 0]])
            else:
                inter = descartes.PolygonPatch(inter)
            self.current_cc_ax.add_patch(poly)
            self.current_cc_ax.add_patch(inter)
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

    # update lower canvas on the overall view.
    def update_lower_canvas(self):
        target_tpa = [unit.get_target_structure().tpa for unit in self.rx_units]
        target_ba = [unit.get_target_structure().ba for unit in self.rx_units]
        target_mcs = [unit.get_target_structure().mcs for unit in self.rx_units]
        target_cc = [unit.get_target_structure().cc for unit in self.rx_units]

        self.lower_ba_ax.cla()
        self.lower_mcs_ax.cla()
        self.lower_cc_ax.cla()

        ba_label = r"$feet^2\ ac^{-1}$"
        tpa_label = r"$ trees\ ac^{-1}$"

        self.lower_ba_ax.set_title("Target Basal Area")
        self.lower_ba_ax.set_ylabel("Basal Area (%s)" % ba_label)
        self.lower_ba_ax.set_xlabel("Density (%s)" % tpa_label)

        self.lower_mcs_ax.set_title("Target Mean Clump Size")
        self.lower_mcs_ax.set_ylabel("Mean Clump Size (n trees)")
        self.lower_mcs_ax.set_xlabel("Density (%s)" % tpa_label)

        self.lower_cc_ax.set_title("Target Canopy Cover")
        self.lower_cc_ax.set_ylabel("Canopy Cover (%)")
        self.lower_cc_ax.set_xlabel("Density (%s)" % tpa_label)

        seaborn.kdeplot(ax=self.lower_ba_ax, x=self.ref_tpa[self.ref_ind['ba']], y=self.ref_ba[self.ref_ind['ba']],
                        cmap="Oranges", fill=True, bw_adjust=0.5)
        self.target_ba_points, = self.lower_ba_ax.plot(target_tpa, target_ba, 'b^')

        seaborn.kdeplot(ax=self.lower_mcs_ax, x=self.ref_tpa[self.ref_ind['mcs']], y=self.ref_mcs[self.ref_ind['mcs']],
                        cmap="Oranges", fill=True, bw_adjust=0.5)
        self.target_mcs_points, = self.lower_mcs_ax.plot(target_tpa, target_mcs, 'b^')

        seaborn.kdeplot(ax=self.lower_cc_ax, x=self.ref_tpa[self.ref_ind['cc']],y=self.ref_cc[self.ref_ind['cc']],
                        cmap="Oranges", fill=True, bw_adjust=0.5)
        self.target_cc_points, = self.lower_cc_ax.plot(target_tpa, target_cc, 'b^')

        self.lower_canvas.draw()

    # Parent function for figuring out what to draw on the treat view.
    def update_raster_canvas(self):
        row = self.stand_tab_list_view.currentIndex().row()
        unit = self.rx_units[row]
        unit._tao_data.exportRxToDll(unit)
        if self.page_counter == 0:
            self.raster_button.setEnabled(True)
            self._draw_raster(unit)
        elif self.page_counter == 1:
            self.raster_button.setEnabled(False)
            self._draw_report_page(unit)
        else:
            self.raster_button.setEnabled(False)
            self._draw_cut_report_page(unit)

    # This draws the view that shows the treatment visualization rasters.
    def _draw_raster(self, unit):
        self.raster_widget.setVisible(True)
        self.treat_display.setVisible(False)
        self.cut_trees_display.setVisible(False)

        clicked = self.raster_button.isChecked()
        viewmode = self.raster_viewmode.currentIndex()

        if clicked:
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
            colors = ("black", "#7bc043", "#fdf498", "#f37736", "#ee4035")
            cm_name = "Clump Colors"
            n_bin = 5
            cm = LinearSegmentedColormap.from_list(cm_name, colors, n_bin)
            n_bin_ranges = (-0.5, 0.5, 1.5, 4.5, 9.5, 99)
            norm = BoundaryNorm(n_bin_ranges, len(n_bin_ranges))
            img = self.raster_ax.imshow(clumps.values, cmap=cm, norm=norm)
            self.raster_ax.imshow(self.treatment.hill, cmap='Greys', alpha=0.5)
            cb_label = "Clump Map (Clump bins)"

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
        self.raster_ax.set_title(
            self.project_settings.get_name() + ', ' + str(unit.get_name()) + ' ' + str(self.raster_viewmode.currentText()))
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

        pts = unit.get_treat_points()
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

    # hover prep for current ba.
    def current_ba_hover(self, event):
        #clean slate
        self.current_ba_annot.set_visible(False)
        for i in range(len(self.unit_ba_polygons)):
            self.unit_ba_polygons[i].set_visible(False)
            self.unit_ba_intersections[i].set_visible(False)

        #check what to show
        if event.inaxes == self.current_ba_ax:
            cont_2, ind_2 = self.unit_ba_points.contains(event)
            if cont_2:
                self.update_annot(ind_2, self.unit_names, self.unit_ba_points, self.current_ba_annot)
                self.unit_ba_polygons[ind_2['ind'][0]].set_visible(True)
                self.unit_ba_intersections[ind_2['ind'][0]].set_visible(True)
                self.current_ba_annot.set_visible(True)
        self.upper_canvas.draw_idle()

    def current_mcs_hover(self, event):
        # clean slate
        self.current_mcs_annot.set_visible(False)
        for i in range(len(self.unit_ba_polygons)):
            self.unit_mcs_polygons[i].set_visible(False)
            self.unit_mcs_intersections[i].set_visible(False)

        # check what to show
        if event.inaxes == self.current_mcs_ax:
            cont_2, ind_2 = self.unit_mcs_points.contains(event)
            if cont_2:
                self.update_annot(ind_2, self.unit_names, self.unit_mcs_points, self.current_mcs_annot)
                self.unit_mcs_polygons[ind_2['ind'][0]].set_visible(True)
                self.unit_mcs_intersections[ind_2['ind'][0]].set_visible(True)
                self.current_mcs_annot.set_visible(True)
        self.upper_canvas.draw_idle()

    def current_cc_hover(self, event):
        # clean slate
        self.current_cc_annot.set_visible(False)
        for i in range(len(self.unit_cc_polygons)):
            self.unit_cc_polygons[i].set_visible(False)
            self.unit_cc_intersections[i].set_visible(False)

        # check what to show
        if event.inaxes == self.current_cc_ax:
            cont_2, ind_2 = self.unit_cc_points.contains(event)
            if cont_2:
                self.update_annot(ind_2, self.unit_names, self.unit_cc_points, self.current_cc_annot)
                self.unit_cc_polygons[ind_2['ind'][0]].set_visible(True)
                self.unit_cc_intersections[ind_2['ind'][0]].set_visible(True)
                self.current_cc_annot.set_visible(True)
        self.upper_canvas.draw_idle()

    def lower_ba_hover(self, event):
        # clean slate
        self.lower_ba_annot.set_visible(False)

        # check what to show
        if event.inaxes == self.lower_ba_ax:
            cont_1, ind_1 = self.targ_ref_ba_points.contains(event)
            cont_2, ind_2 = self.target_ba_points.contains(event)
            if cont_2:
                self.update_annot(ind_2, self.unit_names, self.target_ba_points, self.lower_ba_annot)
                self.lower_ba_annot.set_visible(True)
            elif cont_1:
                self.update_annot(ind_1, self.ref_names, self.targ_ref_ba_points, self.lower_ba_annot)
                self.lower_ba_annot.set_visible(True)
        self.lower_canvas.draw_idle()

    def lower_mcs_hover(self, event):
        vis = self.lower_mcs_annot.get_visible()
        if event.inaxes == self.lower_mcs_ax:
            cont_1, ind_1 = self.target_mcs_points.contains(event)
            cont_2, ind_2 = self.targ_ref_mcs_points.contains(event)
            if cont_1:
                self.update_annot(ind_1, self.unit_names, self.target_mcs_points, self.lower_mcs_annot)
                self.lower_mcs_annot.set_visible(True)
                self.unit_mcs_polygons[ind_1['ind'][0]].set_visible(True)
                self.lower_canvas.draw_idle()
            elif cont_2:
                self.update_annot(ind_2, self.ref_names, self.targ_ref_mcs_points, self.lower_mcs_annot)
                self.lower_mcs_annot.set_visible(True)
                self.lower_canvas.draw_idle()
            else:
                if vis:
                    self.lower_mcs_annot.set_visible(False)
                    for i in range(len(self.unit_mcs_polygons)):
                        self.unit_mcs_polygons[i].set_visible(False)
                    self.lower_canvas.draw_idle()
        else:
            self.lower_mcs_annot.set_visible(False)
            for i in range(len(self.unit_mcs_polygons)):
                self.unit_mcs_polygons[i].set_visible(False)
            self.lower_canvas.draw_idle()

    def lower_cc_hover(self, event):
        vis = self.lower_cc_annot.get_visible()
        if event.inaxes == self.lower_cc_ax:
            cont_1, ind_1 = self.target_cc_points.contains(event)
            cont_2, ind_2 = self.targ_ref_cc_points.contains(event)
            if cont_1:
                self.update_annot(ind_1, self.unit_names, self.target_cc_points, self.lower_cc_annot)
                self.lower_cc_annot.set_visible(True)
                self.unit_cc_polygons[ind_1['ind'][0]].set_visible(True)
                self.lower_canvas.draw_idle()
            elif cont_2:
                self.update_annot(ind_2, self.ref_names, self.targ_ref_cc_points, self.lower_cc_annot)
                self.lower_cc_annot.set_visible(True)
                self.lower_canvas.draw_idle()
            else:
                if vis:
                    self.lower_cc_annot.set_visible(False)
                    for i in range(len(self.unit_cc_polygons)):
                        self.unit_cc_polygons[i].set_visible(False)
                    self.lower_canvas.draw_idle()
        else:
            self.lower_cc_annot.set_visible(False)
            for i in range(len(self.unit_cc_polygons)):
                self.unit_cc_polygons[i].set_visible(False)
            self.lower_canvas.draw_idle()

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

    # This exports the current "raster canvas" whether its the treatment view, or one of the report pages.
    # TODO: Potentially split this off into it's own activity.
    def export_tif(self):
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
    # TODO: Add flags for when qdatawidgetmapper queries so only return strings then.
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
                if self.value_is_valid(value, row, col):
                    self._data[row][col] = value
                    c = col - 5
                    if c == 3:
                        c = 4
                    self._rxunits[row].set_target_structure_by_idx(c, value)
                    # self._locks[row][col - 5] = True
                    # self.update_and_guess_values(row)
                    self.dataChanged.emit(index, index)
                    return True
        except Exception as e:
            print(e)
            pass
        self.dataChanged.emit(index, index)
        return False

    # TODO Implement this.
    def value_is_valid(self, value, row, col):
        ''' ds = self._decision_space[row]
        vals = [ss[col] for ss in ds]
        if min(vals) <= value <= max(vals):
            return True
        else:
            return False '''
        return True

    # TODO. When impemented this would autofill unset targets based on set targets.
    def update_and_guess_values(self, row, col):
        # ts = self._rxunits[row].get_target_structure
        pass


# This class is the widget that srufaces the target and current structure from the datamodel.
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

    # Unused
    def get_checked(self):
        return list(self.tpa_checkbox.isChecked(), self.ba_checkbox.isChecked(), self.mcs_checkbox.isChecked(),
                    self.cc_checkbox.isChecked())


# Dummy class to store the treatment info (for multithreading in theory).
class Treatment:
    def __init__(self, chm, hill, basin, tao_pts):
        self.chm = chm
        self.hill = hill
        self.basin = basin
        self.tao_pts = tao_pts
