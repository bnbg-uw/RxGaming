import fiona
from fiona.transform import transform_geom
import shapely.geometry
import shapely.ops
import shapely.wkt
import collections
import numpy as np
import os
import copy

import raster as raster
import osgeo.osr as osr
from osgeo import gdal
import ctypes
from datetime import datetime
import traceback

class ProjectSettings:
    def __init__(self, name, prj_poly_path, unit_poly_path, ref_db_path, lidar_data_path, unit_name,
                 allometry_coefficients, prop_table_path, dll_path):
        self._name = name
        # repeated lines of code... with set_reference_data()
        try:
            self._ref_data = np.genfromtxt(ref_db_path, delimiter=',', dtype=None, names=True, encoding=None)
        except Exception as exc:
            raise TypeError("Could not read reference data.  Expected a csv") from exc

        try:
            self._mcs_prop = np.genfromtxt(prop_table_path, delimiter=',', dtype=None, names=True, encoding=None)
        except Exception as exc:
            raise TypeError("Could not read MCS proportion data.") from exc

        print(dll_path)
        self.prj_area = ProjectArea(prj_poly_path, unit_poly_path, lidar_data_path, unit_name, allometry_coefficients, dll_path)

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


class ProjectArea:
    def __init__(self, prj_poly_path, unit_poly_path, lidar_data_path, unit_name, allometry_coefficients, dll_path):
        with fiona.open(prj_poly_path, 'r') as shp:
            self._prj_osr = osr.SpatialReference()
            self._prj_osr.ImportFromWkt(shp.meta['crs_wkt'])
            self._prj_poly = list(shp)
        with fiona.open(unit_poly_path, 'r') as shp:
            self._unit_osr = osr.SpatialReference()
            self._unit_osr.ImportFromWkt(shp.meta['crs_wkt'])
            self._unit_poly = list(shp)

        self._tao_data = LidarDataset(lidar_data_path, dll_path)


        shp = [shapely.geometry.shape(feature['geometry']).buffer(0) for feature in self._unit_poly
               if feature['geometry'] is not None]
        shp = [self._tao_data.reprojectPolygon(s, self._unit_osr) for s in shp]
        self._unit_osr = self._tao_data.get_osr()

        prj = shapely.ops.unary_union([shapely.geometry.shape(feature['geometry']) for feature in self._prj_poly
                                       if feature['geometry'] is not None]).buffer(0)
        prj = prj.convex_hull
        prj = self._tao_data.reprojectPolygon(prj, self._prj_osr)
        self._prj_osr = self._tao_data.get_osr()

        if not self.validate_prj_poly(prj) or not self.validate_unit_poly(prj, shp):
            raise ValueError("Invalid unit or project polygon")

        if len(allometry_coefficients) == 0:
            self._tao_data.set_allometry(prj)
        else:
            self._tao_data.set_allometry(allometry_coefficients)
        self._tao_data.doPreProcessing(prj)


        units = [unit.intersection(prj) for unit in shp]
        names = [unit['properties'][unit_name] if unit_name in unit['properties'] else None for unit in self._unit_poly]
        self._units = [RxUnit(unit, self._unit_osr, self._tao_data, name) for unit, name in zip(units, names) if unit.area != 0]

    def prepToPickle(self):
        self._unit_osr = self._unit_osr.ExportToWkt()
        self._prj_osr = self._prj_osr.ExportToWkt()

        for rx in self._units:
            rx._unit_osr = rx._unit_osr.ExportToWkt()

        self._tao_data.prepToPickle()

    def dePickle(self, dll_path):
        tmp = osr.SpatialReference()
        tmp.ImportFromWkt(self._unit_osr)
        self._unit_osr = tmp

        tmp.ImportFromWkt(self._prj_osr)
        self._prj_osr = tmp

        for rx in self._units:
            tmp.ImportFromWkt(rx._unit_osr)
            rx._unit_osr = tmp

        self._tao_data.dePickle(dll_path)

    def add_unit(self, new_rxunit):
        if not isinstance(new_rxunit, RxUnit):
            raise TypeError("Can only add object of class RxUnit to list of units")
        self._units.append(new_rxunit)

    def get_units(self):
        return self._units

    def get_prj_poly(self):
        return self._prj_poly

    def set_prj_poly(self, new_prj_poly):
        if isinstance(new_prj_poly, str):
            with fiona.open(new_prj_poly, 'r') as shp:
                self._prj_poly = list(shp)
        elif isinstance(new_prj_poly, list):
            self._prj_poly = new_prj_poly
        else:
            raise TypeError("new_prj_poly should be a string path to a shapefile or a list of features.")

    def get_tao_data(self):
        return self._tao_data

    def set_tao_data(self, new_lidar_data):
        if isinstance(new_lidar_data, LidarDataset):
            self._tao_data = new_lidar_data
        elif isinstance(new_lidar_data, str):
            self._tao_data = LidarDataset(new_lidar_data)
        else:
            raise TypeError("new_lidar_data should be an object of type LidarDataset or a string path to a dataset")

    # reprojection at some point?
    # logic flow?
    def validate_prj_poly(self, poly):
        x = shapely.ops.unary_union([shapely.geometry.shape(feature['geometry']) for feature in self._tao_data.get_layout_poly()
                                     if feature['geometry'] is not None]).buffer(0)

        if not x.contains(poly) and not poly.equals(x):
            print("prj shape is not within")
            return False

        return True

    #duplicated code?
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


# TODO Let unit_poly be a path to a poly or a shapely shape
class RxUnit:
    def __init__(self, unit_poly, unit_osr, tao_data, name=None):
        print("Name: " + str(name) + " " + str(datetime.now()))
        self._unit_poly = self._validate_polygon(unit_poly)
        self._unit_osr = unit_osr
        self._tao_data = tao_data
        print(1)
        self._tao_data.loadRx(unit_poly, unit_osr)
        print(2)
        self._tao_points = self._tao_data.get_tao_points_dll()
        print(3)
        self._chm = self._tao_data.get_canopy_model_dll()
        print(4)
        self._max_height_map = self._tao_data.get_max_height_map_dll()
        print(5)
        #self._basin_map, self._basin_tiles = self._tao_data.get_basin_map()
        self._basin_map = self._tao_data.get_basin_map_dll()
        print(6)
        self._hillshade = self._calculate_hillshade(self._chm)
        print(7)
        #self._master_nd = neardist.near_dist(self._tao_points[:, 1:3], 6)
        #self._clump_sizes = self._get_raw_clumps(self._tao_points, nd=self._master_nd)
        self._clump_sizes = self._get_raw_clumps_dll()
        print(8)
        self._clump_map = self._make_clump_map_dll(self._clump_sizes)
        print(9)
        self._csd = self._get_csd(self._clump_sizes, group=True)
        self._current_structure = self.get_current_structure_dll()
        self._current_structure = self.get_current_structure_dll()
        self._target_structure = copy.deepcopy(self._current_structure)
        self._treat_taos = None
        self._treat_chm = None
        self._treat_hill = None
        self._treat_basin = None
        self._treat_target_struct = None
        self._treat_best_struct = None
        self._treat_cutoff = None
        self._treat_method = None
        self._treat_clump_sizes = None
        self._treat_clump_map = None
        print(10)
        if name is not None:
            self._name = name
        else:
            self._name = "Unnamed"
        print("Done: " + str(name) + " " + str(datetime.now()))


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
        elif type(polygon) is shapely.geometry.polygon.Polygon:
            return polygon
        elif type(polygon) is shapely.geometry.multipolygon.MultiPolygon:
            return shapely.ops.cascaded_union([
                shapely.geometry.polygon.Polygon(c.exterior).buffer(0.01).buffer(-0.01) for c in polygon
                ])
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
        self._tao_data.dll.setTargetStructure.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        self._tao_data.dll.setTargetStructure.restype = None
        bah = self._target_structure.ba / 4.356
        tph = self._target_structure.tpa * 2.47105
        self._tao_data.dll.setTargetStructure(ctypes.c_double(bah), ctypes.c_double(tph),
                                              ctypes.c_double(self._target_structure.mcs),
                                              ctypes.c_double(self._target_structure.osi),
                                              ctypes.c_double(self._target_structure.cc))

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

    def get_osr(self):
        return self._tao_data.get_osr()

    def get_clump_sizes(self):
        return self._clump_sizes

    def get_clump_map(self):
        return self._clump_map

    def get_treat_clump_map(self):
        return self._treat_clump_map

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
        self._tao_data.dll.getRawClumps.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getRawClumps.restype = None
        cl = np.empty(shape=(self._tao_points.shape[0]), dtype=int)
        self._tao_data.dll.getRawClumps(cl)
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
        self._tao_data.dll.makeClumpMap.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self._tao_data.dll.makeClumpMaprestype = None

        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self._tao_data.dll.getBasinMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres),
                              ctypes.byref(xmin),
                              ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        data = np.zeros((nrow.value, ncol.value), int)
        self._tao_data.dll.makeClumpMap(group_sizes.astype(int), data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=gdal.GDT_Float64)
        return rout

    def get_current_structure_dll(self):
        self._tao_data.dll.getCurrentStructure.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double)]
        self._tao_data.dll.getCurrentStructure.restypes = None
        bah = ctypes.c_double(0)
        tph = ctypes.c_double(0)
        mcs = ctypes.c_double(0)
        osi = ctypes.c_double(0)
        cc = ctypes.c_double(0)
        self._tao_data.dll.getCurrentStructure(ctypes.byref(bah), ctypes.byref(tph), ctypes.byref(mcs), ctypes.byref(osi),
                                               ctypes.byref(cc))
        tpa = tph.value / 2.47105
        baa = bah.value * 4.356
        x = StructureSummary(tpa, baa, mcs.value, osi.value, cc.value)
        return x

    def get_treated_structure_dll(self):
        self._tao_data.dll.getTreatedStructure.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                                           ctypes.POINTER(ctypes.c_double)]
        self._tao_data.dll.getTreatedStructure.restypes = None
        bah = ctypes.c_double(0)
        tph = ctypes.c_double(0)
        mcs = ctypes.c_double(0)
        osi = ctypes.c_double(0)
        cc = ctypes.c_double(0)
        self._tao_data.dll.getTreatedStructure(ctypes.byref(bah), ctypes.byref(tph), ctypes.byref(mcs), ctypes.byref(osi),
                                               ctypes.byref(cc))
        tpa = tph.value / 2.47105
        baa = bah.value * 4.356
        x = StructureSummary(tpa, baa, mcs.value, osi.value, cc.value)
        return x

    def get_simulated_structures_dll(self, bb_dbh=30):
        bb_dbh *= 2.54
        self._tao_data.dll.getSimulatedStructures.argtypes = [ctypes.c_double, np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getSimulatedStructures.restype = None
        data = np.empty(31*5, ctypes.c_double)
        self._tao_data.dll.getSimulatedStructures(ctypes.c_double(bb_dbh), data)
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
        self.set_target_structure_dll();

        self._tao_data.dll.doTreatment.argtypes = [ctypes.c_double, ctypes.c_double]
        self._tao_data.dll.doTreatment.restype = None
        self._tao_data.dll.doTreatment(ctypes.c_double(dbhMin*2.54), ctypes.c_double(dbhMax*2.54))

        self._treat_chm = self.get_treat_chm_dll()
        self._treat_basin = self.get_treat_basin_dll()
        self._treat_taos = self.get_treat_taos_dll()
        self._treat_best_struct = self.get_treated_structure_dll()

        self._treat_hill = copy.deepcopy(self._hillshade)
        self._treat_hill[self._treat_basin == 1] = 200

        self._treat_clump_sizes = self.get_treat_raw_clumps_dll()
        self._treat_clump_map = self.get_treat_clump_map_dll(self._treat_clump_sizes)

        self._treat_cutoff = dbhMax
        print("treatment")
        return self._treat_chm, self._treat_hill, self._treat_taos, self._treat_basin, self._treat_best_struct

    def get_treat_basin_dll(self):
        self._tao_data.dll.getBasinMeta.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
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
        self._tao_data.dll.getBasinMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres),
                              ctypes.byref(xmin),
                              ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self._tao_data.dll.getTreatedBasin.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self._tao_data.dll.getTreatedBasin.restype = None
        data = np.zeros((nrow.value, ncol.value), int)
        self._tao_data.dll.getTreatedBasin(data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=gdal.GDT_Float64)
        return rout

    def get_treat_chm_dll(self):
        self._tao_data.dll.getChmMeta.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
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
        self._tao_data.dll.getChmMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres), ctypes.byref(xmin),
                     ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self._tao_data.dll.getTreatedChm.argtypes = [np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"), ctypes.c_double]
        self._tao_data.dll.getTreatedChm.restype = None
        data = np.empty(shape=(nrow.value, ncol.value), dtype=np.float64)
        self._tao_data.dll.getTreatedChm(data, ctypes.c_double(-9999999))

        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        r = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                          dataType=gdal.GDT_Float64)
        return r

    def get_treat_structure_dll(self):
        self._tao_data.dll.getTreatedStructure.argtypes = [ctypes.POINTER(ctypes.c_double),
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
        self._tao_data.dll.getTreatedStructure(ctypes.byref(bah), ctypes.byref(tph), ctypes.byref(mcs),
                                               ctypes.byref(osi), ctypes.byref(cc))
        tpa = tph.value / 2.47105
        baa = bah.value * 4.356
        x = StructureSummary(tpa, baa, mcs.value, osi.value, cc.value)
        return x

    def get_treat_taos_dll(self):
        self._tao_data.dll.getNTreatedTaos.argtypes = []
        self._tao_data.dll.getNTreatedTaos.restype = int
        x = self._tao_data.dll.getNTreatedTaos()

        self._tao_data.dll.getTreatedTaos.argtypes = [np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getTreatedTaos.restype = None
        taos = np.empty(shape=(x, 5), dtype=np.float64)
        self._tao_data.dll.getTreatedTaos(taos)
        return taos

    def get_treat_raw_clumps_dll(self):
        self._tao_data.dll.getTreatedRawClumps.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        self._tao_data.dll.getTreatedRawClumps.restype = None
        cl = np.empty(shape=(self._tao_points.shape[0]), dtype=int)
        self._tao_data.dll.getTreatedRawClumps(cl)
        return cl

    def get_treat_clump_map_dll(self, group_sizes):
        group_sizes = self._convert_to_clump_sizes(group_sizes)
        print(group_sizes[0:5])
        self._tao_data.dll.getTreatedClumpMap.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
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
        self._tao_data.dll.getBasinMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres),
                              ctypes.byref(xmin),
                              ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        data = np.zeros((nrow.value, ncol.value), int)
        self._tao_data.dll.getTreatedClumpMap(self._treat_basin.values.astype(int), -9999999, group_sizes.astype(int),
                                              data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=gdal.GDT_Float64)
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
        self._osr = None
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
        elif type(polygon) is shapely.geometry.polygon.Polygon:
            return polygon
        elif type(polygon) is shapely.geometry.multipolygon.MultiPolygon:
            return shapely.ops.cascaded_union([
                shapely.geometry.polygon.Polygon(c.exterior).buffer(0.01).buffer(-0.01) for c in polygon
                ])
        else:
            raise TypeError("Unrecognized data format")

    def get_root_path(self):
        return self._root_path

    def set_root_path(self, path):
        self._root_path = path

        lapis = False
        if any("Layout" == f for f in os.listdir(path)):
            lapis = True

        if lapis:
            f = os.listdir(os.path.join(path, "TreeApproximateObjects"))[0]
            if "Feet" in f:
                self.set_units("feet")
            else:
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
        else:
            layout_poly_path = [f for f in os.listdir(os.path.join(path, "Layout_shapefiles"))
                                if f.endswith('ProcessingTiles.shp')][0]
            self.set_layout_poly(os.path.join(path, "Layout_shapefiles", layout_poly_path))

        self.dll = ctypes.CDLL(self._dll_path)
        self.dll.setSeed.restype = None
        self.dll.setSeed(1)

        self.dll.setProjDataDirectory.restype = None
        self.dll.setProjDataDirectory.argtypes = [ctypes.c_char_p]
        #b_dll_path = os.path.dirname(self._dll_path).encode('utf-8')
        b_dll_path = ("C:/OSGeo4W64/share/proj").encode('utf-8')
        #self.dll.setProjDataDirectory(b_dll_path)

        b_root_path = self._root_path.encode('utf-8')
        self.dll.initLidarDataset.argtypes = [ctypes.c_char_p]
        self.dll.initLidarDataset.restype = None
        self.dll.initLidarDataset(b_root_path)
        print("init done")

    def prepToPickle(self):
        out = DllStorage()

        self.dll.getBigMhmMeta.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_char_p]
        self.dll.getBigMhmMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self.dll.getBigMhmMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres),
                            ctypes.byref(xmin),
                            ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self.dll.getBigMhm.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self.dll.getBigMhm.restype = None
        data = np.ndarray((nrow.value, ncol.value), int)
        self.dll.getBigMhm(data, -9999999)

        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        out.bigmhm = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=gdal.GDT_Float64)

        self.dll.getBigChmMeta.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_char_p]
        self.dll.getBigChmMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self.dll.getBigChmMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres),
                            ctypes.byref(xmin),
                            ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self.dll.getBigChm.argtypes = [np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"), ctypes.c_double]
        self.dll.getBigChm.restype = None
        data = np.empty(shape=(nrow.value, ncol.value), dtype=np.float64)
        self.dll.getBigChm(data, ctypes.c_double(-9999999))
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        out.bigchm = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                                   dataType=gdal.GDT_Float64)

        self.dll.getBigBasinMeta.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                          ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                          ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                          ctypes.c_char_p]
        self.dll.getBigBasinMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self.dll.getBigBasinMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres),
                              ctypes.byref(xmin),
                              ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self.dll.getBigBasin.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self.dll.getBigBasin.restype = None
        data = np.zeros((nrow.value, ncol.value), int)
        self.dll.getBigBasin(data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        out.bigbasin = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=gdal.GDT_Float64)

        self.dll.getMaskMeta.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                             ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                             ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                             ctypes.c_char_p]
        self.dll.getMaskMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self.dll.getMaskMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres),
                                 ctypes.byref(xmin),
                                 ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self.dll.getMask.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self.dll.getMask.restype = None
        data = np.zeros((nrow.value, ncol.value), int)
        self.dll.getMask(data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        out.mask = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                                     dataType=gdal.GDT_Float64)

        self.dll.nBigTaos.argtypes = []
        self.dll.nBigTaos.restype = int
        x = self.dll.nBigTaos()

        self.dll.getBigTaos.argtypes = [np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")]
        self.dll.getBigTaos.restype = None
        taos = np.empty(shape=(x, 5), dtype=np.float64)
        self.dll.getBigTaos(taos)
        out.bigtaos = taos

        self.dll.getConvFactor.argtypes = [ctypes.POINTER(ctypes.c_double)]
        self.dll.getConvFactor.restype = None
        convFactor = ctypes.c_double(0)
        self.dll.getConvFactor(ctypes.byref(convFactor))
        out.convFactor = convFactor.value

        self.dll = out
        self._osr = self._osr.ExportToWkt()

    def dePickle(self, dll_path):
        if isinstance(self.dll, DllStorage):
            tmp = self.dll
            self.dll = ctypes.CDLL(dll_path)
            self._dll_path = dll_path
            self.dll.setSeed.restype = None
            self.dll.setSeed(1)

            b_root_path = self._root_path.encode('utf-8')
            self.dll.initLidarDataset.argtypes = [ctypes.c_char_p]
            self.dll.initLidarDataset.restype = None
            self.dll.initLidarDataset(b_root_path)

            self.dll.setBigChm.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int,
                                        ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                        ctypes.c_double, ctypes.c_char_p]
            self.dll.setBigChm.restype = None
            r = tmp.bigchm
            crsWkt = r.projection
            crsWkt = crsWkt.encode('utf-8')
            self.dll.setBigChm(r.values, r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, ctypes.c_double(-9999999), crsWkt)

            self.dll.setBigBasin.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int,
                                          ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                          ctypes.c_int, ctypes.c_char_p]
            self.dll.setBigBasin.restype = None
            r = tmp.bigbasin
            crsWkt = r.projection
            crsWkt = crsWkt.encode('utf-8')
            self.dll.setBigBasin(r.values.astype(int), r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, -9999999, crsWkt)

            self.dll.setBigMhm.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int,
                                             ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                             ctypes.c_double,
                                             ctypes.c_int, ctypes.c_char_p]
            self.dll.setBigMhm.restype = None
            r = tmp.bigmhm
            crsWkt = r.projection
            crsWkt = crsWkt.encode('utf-8')
            self.dll.setBigMhm(r.values.astype(int), r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, -9999999, crsWkt)

            self.dll.setMask.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int,
                                           ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                           ctypes.c_double,
                                           ctypes.c_int, ctypes.c_char_p]
            self.dll.setMask.restype = None
            r = tmp.mask
            crsWkt = r.projection
            crsWkt = crsWkt.encode('utf-8')
            self.dll.setMask(r.values.astype(int), r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, -9999999, crsWkt)

            self.dll.setConvFactor.argtypes = [ctypes.c_double]
            self.dll.setConvFactor.restype = None
            self.dll.setConvFactor(tmp.convFactor)


            self.dll.setBigTaos.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int]
            self.dll.setBigTaos.restype = None
            self.dll.setBigTaos(tmp.bigtaos, tmp.bigtaos.size)

            self.set_allometry(self._allometry)

            tmp = osr.SpatialReference()
            tmp.ImportFromWkt(self._osr)
            self._osr = tmp

    def doPreProcessing(self, projPoly):
        self.dll.doPreProcessing.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
        self.dll.doPreProcessing.restype = None
        crsWkt = self._osr.ExportToWkt()
        print("dopreprocessing")
        wkt = shapely.wkt.dumps(projPoly)
        crsWkt = crsWkt.encode('utf-8')
        wkt = wkt.encode('utf-8')
        self.dll.doPreProcessing(wkt, crsWkt, 1)

    def reprojectPolygon(self, poly, osr):
        self.dll.reprojectPolygon.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
        self.dll.reprojectPolygon.restype = None

        crsWkt = osr.ExportToWkt().encode('utf-8')
        wkt = shapely.wkt.dumps(poly).encode('utf-8')

        outWkt = ctypes.create_string_buffer(2147483647)

        self.dll.reprojectPolygon(wkt, crsWkt, outWkt)
        outWkt = outWkt.value.decode("utf-8")
        return shapely.wkt.loads(outWkt)


    def loadRx(self, poly, osr):
        self.dll.loadRx.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        self.dll.loadRx.restype = None
        crsWkt = osr.ExportToWkt()
        wkt = shapely.wkt.dumps(poly)
        crsWkt = crsWkt.encode('utf-8')
        wkt = wkt.encode('utf-8')
        self.dll.loadRx(wkt, crsWkt)

    def exportRxToDll(self, rx):
        self.dll.setMhm.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int,
                                    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                    ctypes.c_int, ctypes.c_char_p]
        self.dll.setMhm.restype = None
        r = rx._max_height_map
        crsWkt = r.projection
        crsWkt = crsWkt.encode('utf-8')
        self.dll.setMhm(r.values.astype(int), r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, -9999999, crsWkt)
        print(1)
        self.dll.setBasin.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int,
                                      ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                      ctypes.c_int, ctypes.c_char_p]
        self.dll.setBasin.restype = None
        r = rx._basin_map
        crsWkt = r.projection
        crsWkt = crsWkt.encode('utf-8')
        self.dll.setBasin(r.values.astype(int), r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, -9999999, crsWkt)
        print(2)

        self.dll.setChm.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int,
                                    ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                    ctypes.c_double, ctypes.c_char_p]
        self.dll.setChm.restype = None
        r = rx._chm
        crsWkt = r.projection
        crsWkt = crsWkt.encode('utf-8')
        self.dll.setChm(r.values, r.nrow, r.ncol, r.xres, r.yres, r.xmin, r.ymin, ctypes.c_double(-9999999), crsWkt)
        print(3)

        self.dll.importRx.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int,
                                      ctypes.c_double]
        self.dll.importRx.restype = None
        self.dll.importRx(rx._tao_points, rx._tao_points.size, ctypes.c_double(rx._current_structure.osi))
        print(4)

    def get_units(self):
        return self._units

    def set_units(self, unit):
        self._units = unit.lower()

    def get_layout_poly(self):
        return self._layout_poly

    def set_layout_poly(self, path):
        with fiona.open(path, 'r') as shp:
            self._osr = osr.SpatialReference()
            self._osr.ImportFromWkt(shp.meta['crs_wkt'])
            self._layout_poly = list(shp)

    def get_osr(self):
        return self._osr

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
        self.dll.setAllometry.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double]

        if isinstance(a, rxgaming.projectsettings.Allometry):
            self._allometry = a
            print(str(a.intercept) + " " + str(a.slope) + " " + str(a.transform) + " " + str(a.backtransform))
            self.dll.setAllometry(ctypes.c_double(a.intercept), ctypes.c_double(a.slope),
                                  ctypes.c_int(a.transform), ctypes.c_double(a.backtransform))
            return
        try:
            print("Trying allometry from project shape")
            self.dll.setAllometryFiaPath.argtypes = [ctypes.c_char_p]
            self.dll.setAllometryFiaPath.restype = None
            fiaPath = os.path.dirname(self._dll_path) + "/fia/"
            print(fiaPath)
            fiaPath = fiaPath.encode('utf-8')
            self.dll.setAllometryFiaPath(fiaPath)
            print("Set fia path")

            self.dll.setAllometryWkt.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
            self.dll.setAllometryWkt.restype = None
            crsWkt = self._osr.ExportToWkt()
            wkt = shapely.wkt.dumps(a)
            crsWkt = crsWkt.encode('utf-8')
            wkt = wkt.encode('utf-8')
            self.dll.setAllometryWkt(wkt, crsWkt)
            print("Set allometry wkt")

            self.dll.getAllometry.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double)]
            self.dll.getAllometry.restype = None
            slope = ctypes.c_double(0)
            intercept = ctypes.c_double(0)
            transform = ctypes.c_int(0)
            backtransform = ctypes.c_double(0)
            self.dll.getAllometry(ctypes.byref(intercept), ctypes.byref(slope), ctypes.byref(transform), ctypes.byref(backtransform))
            self._allometry = Allometry(intercept.value, slope.value, transform.value, backtransform.value)
            print("get allometry")

        except BaseException as e:
            print(traceback.format_exc())
            raise ValueError("Provide an Allometry object or a wkt to project area")


    def get_max_height_map_dll(self):
        self.dll.getMhmMeta.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.c_char_p]
        self.dll.getMhmMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self.dll.getMhmMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres), ctypes.byref(xmin),
                     ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self.dll.getMhm.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self.dll.getMhm.restype = None
        data = np.ndarray((nrow.value, ncol.value), int)
        self.dll.getMhm(data, -9999999)

        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=gdal.GDT_Float64)
        return rout

    def get_basin_map_dll(self):
        self.dll.getBasinMeta.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.c_char_p]
        self.dll.getBasinMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self.dll.getBasinMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres), ctypes.byref(xmin),
                     ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self.dll.getBasin.argtypes = [np.ctypeslib.ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), ctypes.c_int]
        self.dll.getBasin.restype = None
        data = np.zeros((nrow.value, ncol.value), int)
        self.dll.getBasin(data, -9999999)
        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        rout = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                             dataType=gdal.GDT_Float64)
        return rout

    def get_canopy_model_dll(self):
        self.dll.getChmMeta.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.c_char_p]
        self.dll.getChmMeta.restype = None
        nrow = ctypes.c_int(0)
        ncol = ctypes.c_int(0)
        xres = ctypes.c_double(0)
        yres = ctypes.c_double(0)
        xmin = ctypes.c_double(0)
        ymin = ctypes.c_double(0)
        crsWkt = ctypes.create_string_buffer(9999)
        self.dll.getChmMeta(ctypes.byref(nrow), ctypes.byref(ncol), ctypes.byref(xres), ctypes.byref(yres), ctypes.byref(xmin),
                     ctypes.byref(ymin), crsWkt)
        crsWkt = crsWkt.value.decode("utf-8")

        self.dll.getChm.argtypes = [np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"), ctypes.c_double]
        self.dll.getChm.restype = None
        data = np.empty(shape=(nrow.value, ncol.value), dtype=np.float64)
        self.dll.getChm(data, ctypes.c_double(-9999999))

        e = raster.Extent(xmin.value, xmin.value + xres.value * ncol.value, ymin.value,
                          ymin.value + yres.value * nrow.value)
        r = raster.Raster(e, [xres.value, yres.value], ncol.value, nrow.value, crsWkt, data=data,
                          dataType=gdal.GDT_Float64)
        return r

    def get_tao_points_dll(self):
        self.dll.nCurrentTaos.argtypes = []
        self.dll.nCurrentTaos.restype = int
        x = self.dll.nCurrentTaos()

        self.dll.getCurrentTaos.argtypes = [np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")]
        self.dll.getCurrentTaos.restype = None
        taos = np.empty(shape=(x, 5), dtype=np.float64)
        self.dll.getCurrentTaos(taos)
        return taos


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
        self.backtransform = 0

    def __init__(self, i, s, t, b):
        self.slope = s
        self.intercept = i
        self.transform = t
        self.backtransform = b


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

class DllStorage:
    def __init__(self):
        self.lidarDatasetPath = None
        self.allometry = None
        self.bigmhm = None
        self.bigchm = None
        self.bigosinum = None
        self.bigosiden = None
        self.bigbasin = None
        self.bigtaos = None
        self.mask = None
        self.convFactor = None
