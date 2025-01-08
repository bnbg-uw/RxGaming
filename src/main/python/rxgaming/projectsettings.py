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

    # get the data ready to save to disk.
    def prepToPickle(self):
        self._tao_data.prepToPickle()

    # get the data rady to use when reading from disk.
    def dePickle(self, dll_path):
        self._tao_data.dePickle(dll_path)
        self._tao_data.dll.setRxsSize.argtypes = [ctypes.c_int]
        self._tao_data.dll.setRxsSize.restype = None
        self._tao_data.dll.setRxsSize(len(self._units))

    def add_unit(self, new_rxunit):
        if not isinstance(new_rxunit, RxUnit):
            raise TypeError("Can only add object of class RxUnit to list of units")
        self._units.append(new_rxunit)

    def get_units(self):
        return self._units

    def get_tao_data(self):
        return self._tao_data

    def set_tao_data(self, new_lidar_data):
        if isinstance(new_lidar_data, LidarDataset):
            self._tao_data = new_lidar_data
        elif isinstance(new_lidar_data, str):
            self._tao_data = LidarDataset(new_lidar_data)
        else:
            raise TypeError("new_lidar_data should be an object of type LidarDataset or a string path to a dataset")

    #unused code 10/2/24 - delete?
    def validate_unit_poly(self, prj, unit):
        if not len(self._unit_poly) > 0:
            print("unit len")
            return False

        layout = shapely.ops.unary_union([shapely.geometry.shape(feature['geometry']) for feature in
                                          self._tao_data.get_layout_poly() if feature['geometry'] is not None])

        if all([geom.disjoint(prj) for geom in unit]):
            print("unit disjoint")
            return False
        if any([not geom.is_valid for geom in unit]):
            print("unit invalid")
            return False
        if any([not geom.within(layout) for geom in unit]):
            print("unit layout")
            return False

        return True


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

    # this method is duplicated in lidar dataset.  Make it better.
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

    def set_target_structure(self, target_structure):
        if type(target_structure) is not StructureSummary:
            raise ValueError("To change whole target structure, give a structure summary object.")
        self._target_structure = target_structure

    def set_target_structure_by_idx(self, index, value):
        try:
            value = float(value)
        except:
            raise ValueError("%s cannot be interpreted as type float" % value)

        self._target_structure[index] = value

    def set_target_structure_dll(self):
        self._tao_data.dll.setTargetStructure.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        self._tao_data.dll.setTargetStructure.restype = None
        bah = self._target_structure.ba / 4.356
        tph = self._target_structure.tpa * 2.47105
        cc = self._target_structure.cc / 100
        self._tao_data.dll.setTargetStructure(ctypes.c_int(self.idx),
                                              ctypes.c_double(bah),
                                              ctypes.c_double(tph),
                                              ctypes.c_double(self._target_structure.mcs),
                                              ctypes.c_double(self._target_structure.osi),
                                              ctypes.c_double(cc))

    def get_target_structure(self):
        return self._target_structure

    def get_current_structure(self):
        return self._current_structure

    def get_unit_polygon(self):
        return self._unit_poly

    def set_unit_polygon(self, unit_poly):
        self._unit_poly = self._validate_polygon(unit_poly)
        self.summarize_structure()

    def get_name(self):
        return self._name

    def set_name(self, new_name):
        self._name = new_name

    def get_chm(self):
        return self._chm

    def get_tao_points(self):
        return self._tao_points

    def get_basin_map(self):
        return self._basin_map;

    def get_csd(self):
        return self._csd

    def get_units(self):
        return self._tao_data.get_units()

    def get_hill(self):
        return self._hillshade

    def get_treat_points(self):
        return self._treat_taos

    def get_treat_structure(self):
        return self._treat_best_struct

    def get_treat_chm(self):
        return self._treat_chm

    def get_tao_data(self):
        return self._tao_data

    def get_wkt(self):
        return self._tao_data.get_wkt()

    def get_clump_sizes(self):
        return self._clump_sizes

    def get_clump_map(self):
        return self._clump_map

    def get_treat_clump_map(self):
        return self._treat_clump_map

    def get_max_height_map_dll(self):
        # Getting the raster metadata so we know what the dll will be sending us
        self._tao_data.dll.getMhmMeta.argtypes = [ctypes.c_int,
                                                  ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                                  ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                  ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                  ctypes.c_char_p]
        self._tao_data.dll.getMhmMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self._tao_data.dll.getMhmMeta(ctypes.c_int(self.idx), ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres),
                                      ctypes.byref(yres), ctypes.byref(xmin), ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        # getting the actual raster data.
        self._tao_data.dll.getMhm.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self._tao_data.dll.getMhm.restype = None
        data = np.ndarray((nrow.value, ncol.value), int)
        self._tao_data.dll.getMhm(ctypes.c_int(self.idx), data, -9999999)


        # creating a raster object using our custom python raster library.
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=np.float64)
        return rout

    def get_basin_map_dll(self):
        self._tao_data.dll.getBasinMeta.argtypes = [ctypes.c_int,
                                                    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                                    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                    ctypes.c_char_p]
        self._tao_data.dll.getBasinMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self._tao_data.dll.getBasinMeta(ctypes.c_int(self.idx), ctypes.byref(nrow), ctypes.byref(ncol),
                                        ctypes.byref(xres), ctypes.byref(yres), ctypes.byref(xmin), ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self._tao_data.dll.getBasin.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self._tao_data.dll.getBasin.restype = None
        data = np.zeros((nrow.value, ncol.value), int)
        self._tao_data.dll.getBasin(ctypes.c_int(self.idx), data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=np.float64)
        return rout

    def get_canopy_model_dll(self):
        self._tao_data.dll.getChmMeta.argtypes = [ctypes.c_int,
                                                  ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                                  ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                  ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                  ctypes.c_char_p]
        self._tao_data.dll.getChmMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self._tao_data.dll.getChmMeta(ctypes.c_int(self.idx), ctypes.byref(nrow), ctypes.byref(ncol),
                                      ctypes.byref(xres), ctypes.byref(yres), ctypes.byref(xmin), ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self._tao_data.dll.getChm.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"), ctypes.c_double]
        self._tao_data.dll.getChm.restype = None
        data = np.empty(shape=(nrow.value, ncol.value), dtype=np.float64)
        self._tao_data.dll.getChm(ctypes.c_int(self.idx), data, ctypes.c_double(-9999999))

        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        r = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                          dataType=np.float64)
        return r

    def get_tao_points_dll(self):
        self._tao_data.dll.nTaos.argtypes = [ctypes.c_int]
        self._tao_data.dll.nTaos.restype = int
        x = self._tao_data.dll.nTaos(self.idx)

        self._tao_data.dll.getTaos.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getTaos.restype = None
        taos = np.empty(shape=(x, 5), dtype=np.float64)
        self._tao_data.dll.getTaos(self.idx, taos)
        return taos

    def get_mask_dll(self):
        self._tao_data.dll.getMaskMeta.argtypes = [ctypes.c_int,
                                                    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                                    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                    ctypes.c_char_p]
        self._tao_data.dll.getMaskMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self._tao_data.dll.getMaskMeta(ctypes.c_int(self.idx), ctypes.byref(nrow), ctypes.byref(ncol),
                                        ctypes.byref(xres), ctypes.byref(yres), ctypes.byref(xmin), ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self._tao_data.dll.getMask.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self._tao_data.dll.getMask.restype = None
        data = np.zeros((nrow.value, ncol.value), int)
        self._tao_data.dll.getMask(ctypes.c_int(self.idx), data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=np.float64)
        return rout

    # This is so the trees can be shaded in the display of the rxunit.
    def _calculate_hillshade(self, grid, dx=None, az=315, elev=45):
        if dx is None:
            dx = grid.res

        az_rad, elev_rad = (360 - az + 90) * np.pi / 180, (90 - elev) * np.pi / 180
        sx, sy = self._calculate_slope(grid, dx)

        aspect_rad = np.arctan2(sy, sx)
        s_mag_rad = np.arctan(np.sqrt(sx**2 + sy**2))

        return 255.0 * (np.cos(elev_rad) * np.cos(s_mag_rad) + np.sin(elev_rad) * np.sin(s_mag_rad) *
                        np.cos(az_rad - aspect_rad))

    # TODO: Deal with nodata.
    def _calculate_slope(self, grid, dx):
        zbc = self._assign_bcs(grid)

        sx = (zbc[1:-1, 2:] - zbc[1:-1, :-2])/(2 * dx)
        sy = (zbc[2:, 1:-1] - zbc[:-2, 1:-1])/(2 * dx)

        return sx, sy

    # To be honest I copied this code off stackoverflow years ago and don't know exactly how this works
    # but it creates the hillshade.
    @staticmethod
    def _assign_bcs(grid):
        try:
            grid = grid.values
        except AttributeError:
            pass

        ny, nx = grid.shape
        zbc = np.empty((ny + 2, nx + 2))
        zbc[1:-1, 1:-1] = grid

        zbc[0, 1:-1] = grid[0, :]
        zbc[-1, 1:-1] = grid[-1, :]
        zbc[1:-1, 0] = grid[:, 0]
        zbc[1:-1, -1] = grid[:, -1]

        return zbc

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

    def _make_clump_map_dll(self, group_sizes):
        group_sizes = self._convert_to_clump_sizes(group_sizes)
        self._tao_data.dll.makeClumpMap.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                                    np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self._tao_data.dll.makeClumpMap.restype = None

        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self._tao_data.dll.getBasinMeta(ctypes.c_int(self.idx), ctypes.byref(nrow), ctypes.byref(ncol),
                                        ctypes.byref(xres), ctypes.byref(yres),
                                        ctypes.byref(xmin), ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        data = np.zeros((nrow.value, ncol.value), int)
        self._tao_data.dll.makeClumpMap(ctypes.c_int(self.idx), group_sizes.astype(int), data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=np.float64)
        return rout

    def get_current_structure_dll(self):
        self._tao_data.dll.getCurrentStructure.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double)]
        self._tao_data.dll.getCurrentStructure.restypes = None
        bah = ctypes.c_double(0)
        tph = ctypes.c_double(0)
        mcs = ctypes.c_double(0)
        osi = ctypes.c_double(0)
        cc = ctypes.c_double(0)
        self._tao_data.dll.getCurrentStructure(ctypes.c_int(self.idx), ctypes.byref(bah), ctypes.byref(tph),
                                               ctypes.byref(mcs), ctypes.byref(osi), ctypes.byref(cc))
        tpa = tph.value / 2.47105
        baa = bah.value * 4.356
        cc = cc.value*100
        x = StructureSummary(tpa, baa, mcs.value, osi.value, cc)
        return x

    def get_treated_structure_dll(self):
        self._tao_data.dll.getTreatedStructure.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
        self._tao_data.dll.getTreatedStructure.restypes = None
        bah = ctypes.c_double(0)
        tph = ctypes.c_double(0)
        mcs = ctypes.c_double(0)
        osi = ctypes.c_double(0)
        cc = ctypes.c_double(0)
        self._tao_data.dll.getTreatedStructure(ctypes.c_int(self.idx), ctypes.byref(bah), ctypes.byref(tph),
                                               ctypes.byref(mcs), ctypes.byref(osi), ctypes.byref(cc))
        tpa = tph.value / 2.47105
        baa = bah.value * 4.356
        cc = cc.value * 100
        x = StructureSummary(tpa, baa, mcs.value, osi.value, cc)
        return x

    def get_simulated_structures_dll(self, bb_dbh=30):
        bb_dbh *= 2.54
        self._tao_data.dll.getSimulatedStructures.argtypes = [ctypes.c_int, ctypes.c_double, np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getSimulatedStructures.restype = None
        data = np.empty(31*5, ctypes.c_double)
        self._tao_data.dll.getSimulatedStructures(ctypes.c_int(self.idx), ctypes.c_double(bb_dbh), data)
        out = []
        for i in range(0, 31, 1):
            baa = copy.deepcopy(data[i*5] * 4.356)
            tpa = copy.deepcopy(data[(i*5)+1] / 2.47105)
            mcs = copy.deepcopy(data[(i*5)+2])
            osi = copy.deepcopy(data[(i*5)+3])
            cc = copy.deepcopy(data[(i*5)+4])
            out.append(StructureSummary(tpa, baa, mcs, osi, cc))
        return out


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

    def get_treated_structure(self):
        return self._treat_best_struct

    def get_treat_basin_dll(self):
        self._tao_data.dll.getBasinMeta.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                                                    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
                                                    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                    ctypes.POINTER(ctypes.c_double), ctypes.c_char_p]
        self._tao_data.dll.getBasinMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self._tao_data.dll.getBasinMeta(ctypes.c_int(self.idx), ctypes.byref(nrow), ctypes.byref(ncol),
                                        ctypes.byref(xres), ctypes.byref(yres),
                                        ctypes.byref(xmin), ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self._tao_data.dll.getTreatedBasin.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self._tao_data.dll.getTreatedBasin.restype = None
        data = np.zeros((nrow.value, ncol.value), int)
        self._tao_data.dll.getTreatedBasin(ctypes.c_int(self.idx), data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=np.float64)
        return rout

    def get_treat_chm_dll(self):
        self._tao_data.dll.getChmMeta.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                                                  ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
                                                  ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                  ctypes.POINTER(ctypes.c_double), ctypes.c_char_p]
        self._tao_data.dll.getChmMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self._tao_data.dll.getChmMeta(ctypes.c_int(self.idx), ctypes.byref(nrow), ctypes.byref(ncol),
                                      ctypes.byref(xres), ctypes.byref(yres), ctypes.byref(xmin), ctypes.byref(ymin),
                                      crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self._tao_data.dll.getTreatedChm.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"), ctypes.c_double]
        self._tao_data.dll.getTreatedChm.restype = None
        data = np.empty(shape=(nrow.value, ncol.value), dtype=np.float64)
        self._tao_data.dll.getTreatedChm(ctypes.c_int(self.idx), data, ctypes.c_double(-9999999))

        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        r = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                          dataType=np.float64)
        return r

    def get_treat_structure_dll(self):
        self._tao_data.dll.getTreatedStructure.argtypes = [ctypes.c_int,
                                                           ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double)]
        self._tao_data.dll.getTreatedStructure.restypes = None
        bah = ctypes.c_double(0)
        tph = ctypes.c_double(0)
        mcs = ctypes.c_double(0)
        osi = ctypes.c_double(0)
        cc = ctypes.c_double(0)
        self._tao_data.dll.getTreatedStructure(ctypes.c_int(self.idx), ctypes.byref(bah), ctypes.byref(tph),
                                               ctypes.byref(mcs), ctypes.byref(osi), ctypes.byref(cc))
        tpa = tph.value / 2.47105
        baa = bah.value * 4.356
        cc = cc.value * 100
        x = StructureSummary(tpa, baa, mcs.value, osi.value, cc)
        return x

    def get_treat_taos_dll(self):
        self._tao_data.dll.getNTreatedTaos.argtypes = [ctypes.c_int]
        self._tao_data.dll.getNTreatedTaos.restype = int
        x = self._tao_data.dll.getNTreatedTaos(ctypes.c_int(self.idx))

        self._tao_data.dll.getTreatedTaos.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getTreatedTaos.restype = None
        taos = np.empty(shape=(x, 5), dtype=np.float64)
        self._tao_data.dll.getTreatedTaos(ctypes.c_int(self.idx), taos)
        return taos

    def get_cut_taos_dll(self):
        self._tao_data.dll.getNCutTaos.argtypes = [ctypes.c_int]
        self._tao_data.dll.getNCutTaos.restype = int
        x = self._tao_data.dll.getNCutTaos(ctypes.c_int(self.idx))

        self._tao_data.dll.getCutTaos.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getCutTaos.restype = None
        taos = np.empty(shape=(x, 5), dtype=np.float64)
        self._tao_data.dll.getCutTaos(ctypes.c_int(self.idx), taos)
        return taos

    def get_treat_raw_clumps_dll(self):
        self._tao_data.dll.getTreatedRawClumps.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getTreatedRawClumps.restype = None
        cl = np.empty(shape=(self._tao_points.shape[0]), dtype=int)
        self._tao_data.dll.getTreatedRawClumps(ctypes.c_int(self.idx), cl)
        return cl

    def get_treat_clump_map_dll(self, group_sizes):
        group_sizes = self._convert_to_clump_sizes(group_sizes)
        self._tao_data.dll.getTreatedClumpMap.argtypes = [ctypes.c_int,
                                                          np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                                          ctypes.c_int,
                                                          np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                                          np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                                          ctypes.c_int]
        self._tao_data.dll.getTreatedClumpMap.restype = None

        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self._tao_data.dll.getBasinMeta(ctypes.c_int(self.idx), ctypes.byref(nrow), ctypes.byref(ncol),
                                        ctypes.byref(xres), ctypes.byref(yres),
                                        ctypes.byref(xmin), ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        data = np.zeros((nrow.value, ncol.value), int)
        self._tao_data.dll.getTreatedClumpMap(ctypes.c_int(self.idx), self._treat_basin.values.astype(int), -9999999,
                                              group_sizes.astype(int), data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=np.float64)
        return rout

    def _calculate_treatment_dbh_thin(self, dbh_cutoff, **kwargs):
        dbh = self._tao_data.get_dbh_from_height_dll(self._tao_points[:, 3]) / 2.54
        allowed_to_cut = np.ones_like(dbh, dtype=bool)
        allowed_to_cut[dbh > dbh_cutoff] = False
        self._treatment_cleanup(~allowed_to_cut)
        self._treat_cutoff = dbh_cutoff
        return self._treat_chm, self._treat_hill, self._tao_points[self._treat_mask, :], self._treat_basin, self._treat_best_struct

    def get_treatment(self, method, **kwargs):
        if self._current_structure == self._target_structure and method != 'dbh_thin':
            self._treat_chm = self._chm
            self._treat_hill = self._hillshade
            self._treat_taos = self._tao_points
            self._treat_basin = self._basin_map
            self._treat_best_struct = self._target_structure
            self._treat_cutoff = kwargs['dbh_cutoff']
            return self._treat_chm, self._treat_hill, self._treat_taos, self._treat_basin, self._treat_best_struct
        elif self._target_structure == self._treat_target_struct and self._treat_cutoff == kwargs['dbh_cutoff']:
            return self._treat_chm, self._treat_hill, self._treat_taos, self._treat_basin, self._treat_best_struct
        else:
            self._treat_method = method
            self._treat_target_struct = copy.deepcopy(self._target_structure)
            if method == "add_clumps":
                try:
                    dbh_cutoff = kwargs.pop('dbh_cutoff')
                except KeyError:
                    raise TypeError("Missing keyword arg \'dbh_cutoof\' for method \'add_clumps\'")
                return self.get_treatment_dll(10.16, dbh_cutoff)
            elif method == "dbh_thin":
                try:
                    dbh_cutoff = kwargs.pop('dbh_cutoff')
                except KeyError:
                    raise TypeError("Missing keyword arg \'dbh_cutoff\' for method \'dbh_thin\'.")
                return self._calculate_treatment_dbh_thin(dbh_cutoff, **kwargs)
            else:
                raise ValueError("not implemented")

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
            layout_poly_path = [f for f in os.listdir(os.path.join(path, "Layout"))
                                if f.endswith(".shp")][0]
            self.set_layout_poly(os.path.join(path, "Layout", layout_poly_path))
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
        b_dll_path = (os.path.dirname(self._dll_path)+"/share/proj/").encode('utf-8')
        self.dll.setProjDataDirectory(b_dll_path)

        b_root_path = self._root_path.encode('utf-8')
        self.dll.initLidarDataset.argtypes = [ctypes.c_char_p]
        self.dll.initLidarDataset.restype = bool
        success = self.dll.initLidarDataset(b_root_path)
        if not success:
            raise RuntimeError("Unable to read the lidardataset path as a lidar dataset. "
                               "Is it formatted correctly, with complete data?")
        print("init done")

    def prepToPickle(self):
        self.dll.getConvFactor.argtypes = [ctypes.POINTER(ctypes.c_double)]
        self.dll.getConvFactor.restype = None
        convFactor = ctypes.c_double(0)
        self.dll.getConvFactor(ctypes.byref(convFactor))
        cf = convFactor.value

        self.dll = cf
        print("prepped to pickle")

    def dePickle(self, dll_path):
        if isinstance(self.dll, float):
            cf = self.dll
            self.dll = ctypes.CDLL(dll_path)
            self._dll_path = dll_path
            self.dll.setSeed.restype = None
            self.dll.setSeed(1)

            self.dll.setProjDataDirectory.restype = None
            self.dll.setProjDataDirectory.argtypes = [ctypes.c_char_p]
            b_dll_path = (os.path.dirname(self._dll_path) + "/share/proj/").encode('utf-8')
            self.dll.setProjDataDirectory(b_dll_path)

            b_root_path = self._root_path.encode('utf-8')
            self.dll.initLidarDataset.argtypes = [ctypes.c_char_p]
            self.dll.initLidarDataset.restype = None
            self.dll.initLidarDataset(b_root_path)

            self.dll.setConvFactor.argtypes = [ctypes.c_double]
            self.dll.setConvFactor.restype = None
            self.dll.setConvFactor(cf)

            self.set_allometry(self._allometry)
            print("depickled")

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

    def reprojectPolygon(self, poly, crsWkt):
        self.dll.reprojectPolygon.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
        self.dll.reprojectPolygon.restype = None

        crsWkt = crsWkt.encode('utf-8')
        wkt = shapely.wkt.dumps(poly).encode('utf-8')

        outWkt = ctypes.create_string_buffer(2147483647)

        self.dll.reprojectPolygon(wkt, crsWkt, outWkt)
        outWkt = outWkt.value.decode("utf-8")
        return shapely.wkt.loads(outWkt)
    def exportRxToDll(self, rx):
        self.dll.setMhm.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int,
                                    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                    ctypes.c_int, ctypes.c_char_p]
        self.dll.setMhm.restype = None
        r = rx._max_height_map
        crsWkt = r.projection
        crsWkt = crsWkt.encode('utf-8')
        self.dll.setMhm(rx.idx, r.values.astype(int), r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, -9999999, crsWkt)

        self.dll.setBasin.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int,
                                      ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                      ctypes.c_int, ctypes.c_char_p]
        self.dll.setBasin.restype = None
        r = rx._basin_map
        crsWkt = r.projection
        crsWkt = crsWkt.encode('utf-8')
        self.dll.setBasin(rx.idx, r.values.astype(int), r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, -9999999, crsWkt)

        self.dll.setChm.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int,
                                    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                    ctypes.c_double, ctypes.c_char_p]
        self.dll.setChm.restype = None
        r = rx._chm
        crsWkt = r.projection
        crsWkt = crsWkt.encode('utf-8')
        self.dll.setChm(rx.idx, r.values, r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, ctypes.c_double(-9999999), crsWkt)

        self.dll.setMask.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int,
                                    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                    ctypes.c_double, ctypes.c_char_p]
        self.dll.setMask.restype = None
        r = rx._mask
        crsWkt = r.projection
        crsWkt = crsWkt.encode('utf-8')
        self.dll.setMask(rx.idx, r.values.astype(int), r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, -9999999, crsWkt)

        self.dll.setTaos.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int]
        self.dll.setTaos.restype = None
        self.dll.setTaos(rx.idx, rx._tao_points.astype(float), rx._tao_points.size)

        self.dll.calcCurrentStructure.argtypes = [ctypes.c_int]
        self.dll.calcCurrentStructure.restype = None
        self.dll.calcCurrentStructure(rx.idx)

    def get_units(self):
        return self._units

    def set_units(self, unit):
        self._units = unit.lower()

    def get_layout_poly(self):
        return self._layout_poly

    def set_layout_poly(self, path):
        with fiona.open(path, 'r') as shp:
            self._wkt = shp.meta['crs_wkt']
            self._layout_poly = list(shp)

    def get_wkt(self):
        return self._wkt

    def get_dbh_from_height_dll(self, height):
        self.dll.getDbhFromHeight.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                              np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                              ctypes.c_int]
        self.dll.getDbhFromHeight.restype = None
        dbh = np.empty(shape=height.shape, dtype=np.float64)

        height = height.astype(np.float64)
        self.dll.getDbhFromHeight(height, dbh, height.size)
        return dbh

    def set_allometry(self, a):
        self.dll.setAllometry.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_int]

        if isinstance(a, Allometry):
            self._allometry = a
            self.dll.setAllometry(ctypes.c_double(a.intercept), ctypes.c_double(a.slope),
                                  ctypes.c_int(a.transform))
            return
        try:
            print("Trying allometry from project shape.")
            self.dll.setAllometryFiaPath.argtypes = [ctypes.c_char_p]
            self.dll.setAllometryFiaPath.restype = None
            fiaPath = os.path.dirname(self._dll_path) + "/fia/"
            print("Fia data path is: " + fiaPath)
            fiaPath = fiaPath.encode('utf-8')
            self.dll.setAllometryFiaPath(fiaPath)

            self.dll.setAllometryWkt.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
            self.dll.setAllometryWkt.restype = int
            crsWkt = self.get_wkt()
            wkt = shapely.wkt.dumps(a)
            crsWkt = crsWkt.encode('utf-8')
            wkt = wkt.encode('utf-8')
            nplots = self.dll.setAllometryWkt(wkt, crsWkt)
            print("N plots founds: " + str(nplots))

            self.dll.getAllometry.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_int)]
            self.dll.getAllometry.restype = None
            slope = ctypes.c_double(0)
            intercept = ctypes.c_double(0)
            transform = ctypes.c_int(0)
            self.dll.getAllometry(ctypes.byref(intercept), ctypes.byref(slope), ctypes.byref(transform))
            self._allometry = Allometry(intercept.value, slope.value, transform.value)
        except BaseException as e:
            raise ValueError("Provide an Allometry object or a wkt to project area")

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
    def __init__(self, tpa, ba, mcs, osi, cc):
        self._data = [tpa, ba, mcs, osi, cc]

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
    def osi(self):
        return self._data[3]

    @property
    def cc(self):
        return self._data[4]