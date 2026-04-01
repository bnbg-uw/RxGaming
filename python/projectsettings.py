"""
    Copyright (C) 2024  University of Washington
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.

Bryce Bartl-Geller
University of Washington Forest Resilience Lab
12/6/2024

projectsettings.py
This file is the guts of the data and processing as well as the link to the dll part of the tool.
There's alot of important classes in here: ProjectSettings, LidarDataset, ProjectArea, RxUnit, Allometry, StructureSummary
"""

import fiona # shapefiles i/o
import shapely.geometry  # shapefile geometry.
import shapely.ops
import shapely.wkt
import collections
import numpy as np
import os
import copy

import raster  # our custom raster class
import ctypes # for our dll link
from datetime import datetime
import traceback


# This class more or less holds all the data.
class ProjectSettings:
    def __init__(self, name, unit_poly_path, ref_db_path, lidar_data_path, unit_name, nthreads,
                 prop_table_path, dll_path):
        self._name = name

        if ref_db_path != "":
            try:
                self._ref_data = np.genfromtxt(ref_db_path, delimiter=',', dtype=None, names=True, encoding=None)
            except Exception as exc:
                raise TypeError("Could not read reference data.  Expected a csv") from exc
        else:
            self._ref_data = None

        try:
            self._mcs_prop = np.genfromtxt(prop_table_path, delimiter=',', dtype=None, names=True, encoding=None)
        except Exception as exc:
            raise TypeError("Could not read MCS proportion data.") from exc

        print("RxGaming Dll Path Found: " + dll_path)
        self.prj_area = ProjectArea(unit_poly_path, lidar_data_path, unit_name, nthreads, dll_path)

    def get_name(self):
        return self._name

    def set_name(self, name):
        self._name = name

    def get_prj_area(self):
        return self.prj_area

    def set_prj_area(self, new_prj_area):
        if not isinstance(new_prj_area, ProjectArea):
            raise TypeError("new_prj_area should be an object of type ProjectArea")
        self.prj_area = new_prj_area

    def get_reference_data(self):
        return self._ref_data

    def set_reference_data(self, new_ref_data):
        if isinstance(new_ref_data, str):
            try:
                self._ref_data = np.genfromtxt(new_ref_data, delimiter=',', dtype=None, names=True, encoding=None)
            except Exception as exc:
                raise TypeError("Could not read reference data.  Expected a csv") from exc
        elif isinstance(new_ref_data, np.ndarray):
            self._ref_data = new_ref_data
        else:
            raise TypeError("new_ref_data should be a string to a new dataset or a dataset as numpy ndarray.")

    def get_mcs_prop(self):
        return self._mcs_prop

# These are nesting class, the project area holds the lidar data and the rxunits.
class ProjectArea:
    def __init__(self, unit_poly_path, lidar_data_path, unit_name, nthreads, dll_path):
        print("Creating project area.")
        print("Reading unit polygon.")
        with fiona.open(unit_poly_path, 'r') as shp:
            self._unit_wkt = shp.crs_wkt
            self._unit_poly = list(shp)

        print("Reading lidar dataset.")
        self._tao_data = LidarDataset(lidar_data_path, dll_path)

        # this cleans up bad geometry (self intersections, duplicated points, etc.).
        print("Cleaning unit polygon geometry.")
        shp = [shapely.geometry.shape(feature['geometry']).buffer(0) for feature in self._unit_poly
               if feature['geometry'] is not None]
        shp = [self._tao_data.reprojectPolygon(s, self._unit_wkt).buffer(0) for s in shp]
        self._unit_wkt = self._tao_data.get_wkt()

        print("Calculating Height to DBH allometry.")
        merged = shapely.unary_union(shp)
        self._tao_data.set_allometry(merged.convex_hull)
        shp = self._tao_data.doPreProcessing(shp, nthreads)

        if not len(shp):
            raise RuntimeError("No treatment units overlap lidar data.")

        print("Creating RxUnits.")
        units = [unit for unit in shp]
        names = [unit['properties'][unit_name] if unit_name in unit['properties'] else None for unit in self._unit_poly]
        self._units = [RxUnit(unit, self._unit_wkt, self._tao_data, idx, name)
                       for unit, name, idx in zip(units, names, range(len(units))) if unit.area != 0]

# this class is stored in a list of other rxunits.  Each one represents a treatment area.
# most of the communicating with the dll happens here.
class RxUnit:
    def __init__(self, unit_poly, unit_wkt, tao_data, idx, name=None):
        self.idx = idx  # for communicating with the dll which unit this is.
        print("Name: " + str(name) + " " + str(datetime.now()))
        self._unit_poly = self._validate_polygon(unit_poly)
        self._unit_wkt = unit_wkt
        self._tao_data = tao_data  # so we can communicate with the dll.
        self._tao_points = self.get_tao_points_dll()
        self._chm = self.get_canopy_model_dll()
        self._max_height_map = self.get_max_height_map_dll()
        self._basin_map = self.get_basin_map_dll()
        self._mask = self.get_mask_dll()
        self._hillshade = self._calculate_hillshade(self._chm)
        self._clump_sizes = self._get_raw_clumps_dll()
        self._clump_map = self._make_clump_map_dll(self._clump_sizes)
        self._csd = self._get_csd(self._clump_sizes, group=True)
        self._current_structure = self.get_current_structure_dll()
        self._current_structure = self.get_current_structure_dll()
        self._target_structure = copy.deepcopy(self._current_structure)
        self._treat_taos = None
        self.cut_taos = None
        self.treatment_result = None
        self._treat_chm = None
        self._treat_hill = None
        self._treat_basin = None
        self._treat_target_struct = None
        self._treat_best_struct = None
        self._treat_cutoff = None
        self._treat_method = None
        self._treat_clump_sizes = None
        self._treat_clump_map = None
        if name is not None:
            self._name = name
        else:
            self._name = "Unnamed"

    def get_csd(self):
        return self._csd

    def get_clump_sizes(self):
        return self._clump_sizes


    def _get_raw_clumps_dll(self):
        self._tao_data.dll.getRawClumps.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getRawClumps.restype = None
        cl = np.empty(shape=(self._tao_points.shape[0]), dtype=int)
        self._tao_data.dll.getRawClumps(ctypes.c_int(self.idx), cl)
        return cl

    # Convert raw clumps into a record array. 1 row per clump.
    def _get_csd(self, data, clump_thresh_dist=6, nd=None, group=False):
        if data.shape[0] == 0:
            return 0

        if not group and data.shape[1] == 0:
            return 0

        def classify(x, breaks):
            out = np.full_like(x, -1, int)
            for i in range(1, len(breaks)):
                in_bin = (x >= breaks[i - 1]) & (x < breaks[i])
                if len(in_bin) > 0:
                    out[in_bin] = i - 1
            return out
        if not group:
            group = self._get_raw_clumps_dll()
        else:
            group = data
        clump_sizes = self._convert_to_clump_sizes(group)
        bins = classify(clump_sizes, [1, 2, 5, 10, 15, 31, max(clump_sizes) + 1])
        bins = np.array(["Individual", "Small", "Medium", "Large", "Super", "Mega", "Na"])[bins]
        return np.array([(a, b, c) for a, b, c in zip(group, clump_sizes, bins)],
                        dtype=[("clump_id", int), ("clump_size", int), ("clump_bins", str, 10)])

    # return raw clump vector as vector of same length but each clump is the clump size instead of clump id.
    @staticmethod
    def _convert_to_clump_sizes(data):
        all_sizes = dict(zip(*np.unique(data, return_counts=True)))
        clump_sizes = data.copy()
        for key, value in all_sizes.items():
            clump_sizes[data == int(key)] = value
        return clump_sizes


    # return a distribution of BA for trees, by given points.
    def get_ba_dist(self, pts):
        dbh = self._tao_data.get_dbh_from_height_dll(pts[:, 3])
        dbh[dbh < 0] = 1
        return 0.005454 * (dbh / 2.54)**2

    def get_treatment_dll(self, dbhMin, dbhMax):
        self.set_target_structure_dll()

        self._tao_data.dll.doTreatment.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double]
        self._tao_data.dll.doTreatment.restype = None
        self._tao_data.dll.doTreatment(ctypes.c_int(self.idx), ctypes.c_double(dbhMin*2.54), ctypes.c_double(dbhMax*2.54))

        self._treat_chm = self.get_treat_chm_dll()
        self._treat_basin = self.get_treat_basin_dll()
        self._treat_taos = self.get_treat_taos_dll()
        self.cut_taos = self.get_cut_taos_dll()
        self.treatment_result = self.get_treatment_result_dll()
        self._treat_best_struct = self.get_treated_structure_dll()

        self._treat_hill = copy.deepcopy(self._hillshade)

        self._treat_hill[self._treat_basin.values == 1] = 200

        self._treat_clump_sizes = self.get_treat_raw_clumps_dll()
        self._treat_clump_map = self.get_treat_clump_map_dll(self._treat_clump_sizes)

        self._treat_cutoff = dbhMax
        print("treatment")
        return self._treat_chm, self._treat_hill, self._treat_taos, self._treat_basin, self._treat_best_struct

    def get_treatment_result_dll(self):
        self._tao_data.dll.getTreatmentResult.argtypes = [ctypes.c_int]
        self._tao_data.dll.getTreatmentResult.restype = int
        return self._tao_data.dll.getTreatmentResult(ctypes.c_int(self.idx))

    def get_treat_raw_clumps_dll(self):
        self._tao_data.dll.getTreatedRawClumps.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getTreatedRawClumps.restype = None
        cl = np.empty(shape=(self._tao_points.shape[0]), dtype=int)
        self._tao_data.dll.getTreatedRawClumps(ctypes.c_int(self.idx), cl)
        return cl


class LidarDataset:
    def __init__(self, lidar_data_path, dll_path):
        self._dll_path = dll_path
        self._root_path = None
        self._units = None
        self._layout_poly = None
        self._segments_path = None
        self._chm_path = None
        self._allometry = None
        self._wkt = None
        self.dll = None

        self.set_root_path(lidar_data_path)

    @staticmethod
    def _convert_to_shape(shp):
        return [shapely.geometry.shape(feature['geometry']) for feature in shp]

    @staticmethod
    def _validate_polygon(polygon):
        if isinstance(polygon, list):
            if len(polygon) == 1:
                polygon = polygon[0]
            else:
                raise TypeError("Polygon should be a single feature or shape")
        if isinstance(polygon, dict) or isinstance(polygon, collections.OrderedDict):
            return shapely.geometry.shape(polygon['geometry'])
        elif type(polygon) is shapely.geometry.polygon.Polygon or shapely.geometry.multipolygon.MultiPolygon:
            return polygon
        else:
            raise TypeError("Unrecognized data format")

    def get_root_path(self):
        return self._root_path

    def set_root_path(self, path):
        self._root_path = path

        lapis = False
        lidR = False
        if "Layout" in os.listdir(path):
            lapis = True
        if "layout" in os.listdir(path):
            lidR = True

        if lapis:
            f = os.listdir(os.path.join(path, "TreeApproximateObjects"))[0]
            if "Feet" in f:
                self.set_units("feet")
            else:
                self.set_units("meters")
        elif lidR:
            #TODO: UNITS!
            self.set_units("meters")
        else:
            if any("FEET" in f for f in os.listdir(path)):
                self.set_units("feet")
            else:
                self.set_units("meters")

        if lapis:
            self.set_layout_poly(os.path.join(path, "Layout", "TileLayout.shp"))
        elif lidR:
            layout_poly_path = [f for f in os.listdir(os.path.join(path, "layout"))
                                if f.endswith(".shp")][0]
            self.set_layout_poly(os.path.join(path, "layout", layout_poly_path))
        else:
            layout_poly_path = [f for f in os.listdir(os.path.join(path, "Layout_shapefiles"))
                                if f.endswith('ProcessingTiles.shp')][0]
            self.set_layout_poly(os.path.join(path, "Layout_shapefiles", layout_poly_path))

        self.dll = ctypes.CDLL(self._dll_path)
        self.dll.setSeed.restype = None
        self.dll.setSeed(1)

        self.dll.setProjDataDirectory.restype = None
        self.dll.setProjDataDirectory.argtypes = [ctypes.c_char_p]
        b_dll_path = (os.path.dirname(self._dll_path)).encode('utf-8')
        self.dll.setProjDataDirectory(b_dll_path)

        b_root_path = self._root_path.encode('utf-8')
        self.dll.initLidarDataset.argtypes = [ctypes.c_char_p]
        self.dll.initLidarDataset.restype = bool
        success = self.dll.initLidarDataset(b_root_path)
        if not success:
            raise RuntimeError("Unable to read the lidardataset path as a lidar dataset. "
                               "Is it formatted correctly, with complete data?")
        print("init done")

    def doPreProcessing(self, projPoly, nthreads):
        self.dll.queueRx.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        self.dll.queueRx.restype = ctypes.c_bool
        crsWkt = self._wkt
        crsWkt = crsWkt.encode('utf-8')

        result = []
        for shp in projPoly:
            wkt = shapely.wkt.dumps(shp)
            wkt = wkt.encode('utf-8')
            result.append(self.dll.queueRx(wkt, crsWkt))

        self.dll.doPreProcessing.argtypes = [ctypes.c_int]
        self.dll.doPreProcessing.restype = None
        print("Doing preprocessing on " + str(nthreads) + " threads.")
        self.dll.doPreProcessing(nthreads)
        return [shp for shp, res in zip(projPoly, result) if res]

class Allometry:
    def __init__(self):
        self.slope = 0
        self.intercept = 0
        #0 = none
        #1 = sqrt
        #2 = curt
        #3 = log
        #4 = suggest
        self.transform = 0

    def __init__(self, i, s, t):
        transform = ""
        if t == 0:
            transform = "response transform = none"
        elif t == 1:
            transform = "response transform = square"
        elif t == 2:
            transform = "response transform = cube"
        else:
            transform = "transform = log-log"

        print("Allometry: intercept = " + str(round(i, 4)) +
              "; slope = " + str(round(s, 4)) +
              " " + transform)
        self.slope = s
        self.intercept = i
        self.transform = t


class StructureSummary:
    def __init__(self, tpa, ba, mcs, cc):
        self._data = [tpa, ba, mcs, cc]

    def __str__(self):
        return f'''
        TPA = {round(self.tpa, 3)}
        BA = {round(self.ba, 3)}
        MCS = {round(self.mcs, 3)}
        CC = {round(self.cc, 3)}
        '''

    def __eq__(self, other):
        if other is None:
            return False
        else:
            return self.__dict__ == other.__dict__

    def __len__(self):
        return len(self._data)

    def __getitem__(self, i):
        return self._data[i]

    def __setitem__(self, i, value):
        try:
            value = float(value)
        except:
            raise ValueError("Cannot cast %f as type float" % value)
        self._data[i] = value

    @property
    def tpa(self):
        return self._data[0]

    @property
    def ba(self):
        return self._data[1]

    @property
    def mcs(self):
        return self._data[2]

    @property
    def cc(self):
        return self._data[3]