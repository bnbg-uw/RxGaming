__author__ = 'jontkane'
from osgeo import gdal
from osgeo.gdalconst import *
import numbers
import struct
import copy
import numpy
import math
import operator
import os
import rasterio.features
import affine

GDTtoStruct = {gdal.GDT_Float32: 'f', gdal.GDT_Float64: 'd', gdal.GDT_Int16: 'h', gdal.GDT_Int32: 'i'}
DTMtoGDT = {0: gdal.GDT_Int16, 1: gdal.GDT_Int32, 2: gdal.GDT_Float32, 3: gdal.GDT_Float64}
DTMheader = {'xmin': ('d', 86, 94), 'ymin': ('d', 94, 102), 'xres': ('d', 126, 134), 'yres': ('d', 134, 142),
             'ncol': ('i', 142, 146), 'nrow': ('i', 146, 150), 'dataType': ('h', 154, 156)}
DTMstart = 200
GDTtoNP = {gdal.GDT_Float32: numpy.float32, gdal.GDT_Float64: numpy.float64, gdal.GDT_Int16: numpy.int16,
           gdal.GDT_Int32: numpy.int32}


def rasterFromFile(f, band=1, dataType=gdal.GDT_Float32, empty=False):
    if os.path.splitext(f)[1] == ".dtm":
        with open(f, mode='rb') as file:
            content = file.read(200)

        def getHeaderInfo(name):
            x = DTMheader[name]
            s = struct.Struct(x[0])
            return s.unpack(content[x[1]:x[2]])[0]

        ymin = getHeaderInfo('ymin')
        yres = getHeaderInfo('yres')
        ymin = ymin - yres / 2
        nrow = getHeaderInfo('nrow')
        ymax = ymin + yres * nrow
        ncol = getHeaderInfo('ncol')
        dataType = DTMtoGDT[getHeaderInfo('dataType')]
        xmin = getHeaderInfo('xmin')
        xres = getHeaderInfo('xres')
        xmin = xmin - xres / 2
        geo = (xmin, xres, 0, ymax, 0, -1 * yres)
        if empty:
            return rasterFromGeo(geo, nrow, ncol, "", True)
        else:
            with open(f, mode='rb') as file:
                content = file.read()
            data = numpy.rot90(numpy.frombuffer(content[200:], GDTtoNP[dataType]).reshape(ncol, nrow))
            return rasterFromGeo(geo, nrow, ncol, "", False, data, dataType, -1)

    r = gdal.Open(f, GA_ReadOnly)
    geo = r.GetGeoTransform()
    nrow = r.RasterYSize
    ncol = r.RasterXSize
    projection = r.GetProjection()
    if empty:
        return rasterFromGeo(geo, nrow, ncol, projection, True)
    else:
        data = r.GetRasterBand(band).ReadRaster(0, 0, ncol, nrow, ncol, nrow, dataType)
        data = numpy.frombuffer(data, dtype=GDTtoNP[dataType])
        NAvalue = r.GetRasterBand(band).GetNoDataValue()
        if not numpy.min(data) == NAvalue and numpy.min(data) < -1e10:
            NAvalue = numpy.min(data)
        return rasterFromGeo(geo, nrow, ncol, projection, False, data, dataType, NAvalue)


def rasterWithNewValues(r, data, NAvalue=None, dataType=gdal.GDT_Float32):
    geo = (r.xmin, r.xres, 0, r.ymax, 0, -1 * r.yres)
    if NAvalue is None:
        return rasterFromGeo(geo, r.nrow, r.ncol, r.projection, False, data, dataType, r.NAvalue)
    else:
        return rasterFromGeo(geo, r.nrow, r.ncol, r.projection, False, data, dataType, NAvalue)


# A function that accepts gdal-ish input and makes a raster from it, as per Raster's old __init__.
# You shouldn't ever need to use it under normal circumstances, but it exists for all previous calls to Raster() to point to this instead.
def rasterFromGeo(geo, nrow, ncol, projection, empty=False, data=None, dataType=None, NAvalue=None):
    nrow = int(nrow)
    ncol = int(ncol)
    xres = numpy.float32(abs(geo[1]))
    yres = numpy.float32(abs(geo[5]))
    e = Extent(min(geo[0], geo[0] + geo[0] * ncol),
               max(geo[0], geo[0] + geo[1] * ncol),
               min(geo[3], geo[3] + geo[5] * nrow),
               max(geo[3], geo[3] + geo[5] * nrow))
    return Raster(e, (xres, yres), ncol, nrow, projection, empty, data, dataType, NAvalue)


# origin and cellsize are always respected. Extent is expanded to match. specify exactly one of cellsize, nrow, or ncol to make it work. data=None results in an empty raster
# If origin isn't specified, it will align the center with (0,0)
def rasterFromExtent(e, origin=None, res=None, ncol=None, nrow=None, projection="", data=None,
                     dataType=gdal.GDT_Float32, NAvalue=None):
    if not ((res is None) + (ncol is None) + (nrow is None)) == 2:
        raise ValueError("Must specify exactly one of res, ncol, and nrow)")
    if ncol is not None:
        res = (e.xmax - e.xmin) // ncol
    if nrow is not None:
        res = (e.ymax - e.ymin) // nrow
    if origin is None:
        origin = (res / 2, res / 2)
    e = Extent(e.xmin, e.xmax, e.ymin, e.ymax)
    e.xmin = ((e.xmin - origin[0]) // res) * res + origin[0]
    e.xmax = (((e.xmax - origin[0]) // res) + 1) * res + origin[0]
    e.ymin = ((e.ymin - origin[1]) // res) * res + origin[1]
    e.ymax = ((e.ymax - origin[1]) // res + 1) * res + origin[1]
    geo = (e.xmin, res, 0, e.ymax, 0, -1 * res)
    nrow = int((e.ymax - e.ymin) // res)
    ncol = int((e.xmax - e.xmin) // res)
    empty = data is None
    return rasterFromGeo(geo, nrow, ncol, projection, empty, data, dataType, NAvalue)


def _myround(x, base, origin):
    return base * round(float(x - origin) / base) + origin


# x is always self, which is 2-dimensional. y can be a raster, a scalar, a vector, a 2-dimensional matrix, or a rasterindex
def _saferastermathloop(x, y, fun, right=False):
    funequiv = {operator.add: numpy.add, operator.sub: numpy.subtract, operator.mul: numpy.multiply,
                operator.floordiv: numpy.floor_divide, operator.mod: numpy.mod, operator.truediv: numpy.true_divide}
    if fun not in funequiv and fun not in funequiv.values():
        raise ValueError("not a supported function for raster math")
    if fun in funequiv:
        fun = funequiv[fun]

    isdiv = (fun == numpy.true_divide or fun == numpy.floor_divide or fun == numpy.mod)

    v = copy.copy(x.values)
    xNA = numpy.isnan(v)
    if right and isdiv:
        xNA = numpy.logical_or(xNA, v == 0)

    if numpy.isscalar(y):
        nxNA = numpy.logical_not(xNA)
        if right:
            v[nxNA] = fun(y, v[nxNA])
            v[xNA] = numpy.nan
        else:
            v[nxNA] = fun(v[nxNA], y)
            v[xNA] = numpy.nan

    elif isinstance(y, RasterIndex):
        # loop over non-NA cells
        # get best remaining raster
        # if res is equal or worse, move on
        # if res is better, subtract it from all cells that it's better for and record that res, then remove the raster
        # important--if dem is NA for the cell in question, don't remove it from consideration. (Still subtract for cells its good for while it's loaded and then remove it)
        y = RasterIndex(y.data, scale=y.scale, shift=y.shift, projection=y.projection)
        for f in y.filelist:
            if not x.extent.overlaps(y.getempty(f)):
                y.removeRaster(f)

        v = v.flatten()
        newv = numpy.full_like(v, numpy.nan)
        xNA = xNA.flatten()
        res = numpy.array(numpy.full_like(v, numpy.inf))

        def applyRaster(r, newv, res):
            cells = x.cellsFromExtent(r.extent)
            wres = res[cells] > r.res
            rv = r.extract(x.xFromCell(cells, True), x.yFromCell(cells, True))
            wrnotNA = numpy.logical_not(numpy.isnan(rv))
            if not right and isdiv:
                wrnotNA = numpy.logical_and(wrnotNA, numpy.logical_not(rv == 0))
            wuse = numpy.logical_and(wres, wrnotNA)
            if right:
                newv[cells[wuse]] = fun(rv[wuse], v[cells[wuse]])
            else:
                newv[cells[wuse]] = fun(v[cells[wuse]], rv[wuse])
            res[cells[wuse]] = r.res
            return newv, res

        for c in range(x.ncell):
            if xNA[c]:
                continue
            foundr = False
            coords = x.xyFromCell(c)
            while not foundr:
                f, r = y.bestRasterFromXY(coords[0], coords[1], True)
                if r is None or r.res >= res[c]:
                    break
                f, r = y.bestRasterFromXY(coords[0], coords[1], False)
                newv, res = applyRaster(r, newv, res)
                y.removeRaster(f)
                foundr = r.extract(coords[0], coords[1]) is not None
                if len(y.filelist) == 0:
                    break
            if len(y.filelist) == 0:
                break
        return newv
    else:
        if isinstance(y, Raster):
            vy = y.values
            yNA = numpy.isnan(vy)
            if not right and isdiv:
                yNA = numpy.logical_or(yNA, vy == 0)
        else:
            vy = numpy.array(y)
            if not right and isdiv:
                yNA = vy == 0
                yNA = yNA.reshape(v.shape)
            else:
                yNA = numpy.zeros_like(v)

        if vy.ndim == 1:
            if not vy.shape[0] == x.ncell:
                raise ValueError("length of second object is not equal to ncell of raster")
            vy = vy.reshape(v.shape)
        else:
            if not vy.shape == v.shape:
                raise ValueError("shape of second object not same as shape of raster")
        xyNA = numpy.logical_or(yNA, xNA)
        xynNA = numpy.logical_not(xyNA)
        if right:
            v[xynNA] = fun(vy[xynNA], v[xynNA])
            v[xyNA] = numpy.nan
        else:
            v[xynNA] = fun(v[xynNA], vy[xynNA])
            v[xyNA] = numpy.nan
    return v


def merge(*r, method="first", mask=None, NAvalue=None, origin=None, res=None, projection=None,safe=False):
    if method.lower() not in ["first", "last", "inside"]:
        raise ValueError("First, last, and inside are the only supported methods")
    if isinstance(r[0], RasterIndex):
        if len(r) > 1:
            raise ValueError("Only one RasterIndex at a time; please merge them separately and then merge the results")
        r = list(r[0].filelist)
    if method.lower() == "last":
        r = list(r)
        r.reverse()
        method = "first"
    if mask is None:
        origin = None
        res = None
        projection = None
        e = None
        if NAvalue is None:
            NAvalue = -9999
    else:
        if origin is None:
            origin = mask.origin
        else:
            if not numpy.allclose(origin, mask.origin, rtol=0.05):
                raise ValueError("Inputs and mask must all have the same origin")
        if res is None:
            res = mask.res
        else:
            if not numpy.allclose(res, mask.res):
                raise ValueError("Inputs and mask must all have same resolution")
        if projection is None:
            projection = mask.projection
        else:
            if not projection == mask.projection or mask.projection == "":
                raise ValueError("Inputs and mask must all have same or undefined projection")
        e = mask.extent
        if NAvalue is None:
            NAvalue = mask.NAvalue
    for rcurrent in r:
        if isinstance(rcurrent, str):
            rcurrent = rasterFromFile(rcurrent, empty=True)
        if origin is None:
            origin = rcurrent.origin
        if res is None:
            res = rcurrent.res
        if projection is None:
            projection = rcurrent.projection
        if not safe and not numpy.allclose(origin, rcurrent.origin, rtol=0.05):
            print(origin)
            print(rcurrent.origin)
            raise ValueError("Inputs and mask must all have the same origin")
        if not safe and not numpy.allclose(res, rcurrent.res):
            raise ValueError("Inputs and mask must all have same resolution")
        if not safe and not projection == rcurrent.projection or rcurrent.projection == "":
            raise ValueError("Inputs and mask must all have same or undefined projection")
        if mask is None:
            if e is None:
                e = Extent(rcurrent.xmin, rcurrent.xmax, rcurrent.ymin, rcurrent.ymax)
            if rcurrent.xmin < e.xmin:
                e.xmin = rcurrent.xmin
            if rcurrent.xmax > e.xmax:
                e.xmax = rcurrent.xmax
            if rcurrent.ymin < e.ymin:
                e.ymin = rcurrent.ymin
            if rcurrent.ymax > e.ymax:
                e.ymax = rcurrent.ymax
    if mask is None:
        mask = rasterFromExtent(e, origin=origin, res=res, projection=projection)
        mask = rasterWithNewValues(mask, numpy.full((mask.nrow, mask.ncol), NAvalue), NAvalue=NAvalue)
    else:
        mask = rasterWithNewValues(mask, numpy.full((mask.nrow, mask.ncol), NAvalue), NAvalue=NAvalue)
    if method == "inside":
        ydist = numpy.full(mask.ncell, numpy.NINF)
        xdist = numpy.full(mask.ncell, numpy.NINF)
    mv = mask.valuesList()
    for rcurrent in r:
        if isinstance(rcurrent, str):
            rcurrent = rasterFromFile(rcurrent)
        c = mask.cellsFromExtent(rcurrent)
        rv = rcurrent.valuesList()
        if len(c) != rcurrent.ncell:
            raise NotImplementedError(
                "Something went wrong with merge implementation; work with Jonathan to solve this corner case")
        if method == "first":
            use = numpy.logical_not(numpy.isnan(rv))
            use = numpy.logical_and(use, numpy.isnan(mv[c]))
            mv[c[use]] = rv[use]
        elif method == "inside":
            x = mask.xFromCell(c)
            y = mask.yFromCell(c)
            currentxdist = numpy.full(mask.ncell, numpy.inf)
            currentydist = numpy.full(mask.ncell, numpy.inf)
            currentxdist[c] = numpy.minimum(currentxdist[c], x - rcurrent.xmin)
            currentydist[c] = numpy.minimum(currentydist[c], rcurrent.xmax - x)
            currentxdist[c] = numpy.minimum(currentxdist[c], y - rcurrent.ymin)
            currentydist[c] = numpy.minimum(currentydist[c], rcurrent.ymax - y)

            mincurrentdist = numpy.minimum(currentxdist[c], currentydist[c])
            maxcurrentdist = numpy.maximum(currentxdist[c], currentydist[c])
            mindist = numpy.minimum(xdist[c], ydist[c])
            maxdist = numpy.maximum(xdist[c], ydist[c])
            mindisthigher = mincurrentdist > mindist
            mindistequal = mincurrentdist == mindist
            maxdisthigher = maxcurrentdist > maxdist

            use = numpy.logical_and(mindistequal, maxdisthigher)
            use = numpy.logical_or(use, mindisthigher)
            use = numpy.logical_and(use, numpy.logical_not(numpy.isnan(rv)))

            mv[c[use]] = rv[use]
            xdist[c[use]] = currentxdist[c[use]]
            ydist[c[use]] = currentydist[c[use]]
    mask.values = mv
    return mask


def crop(r, e, snap="out"):
    r = copy.copy(r)
    r.crop(e, snap=snap)
    return r


class Extent:
    def __init__(self, xmin, xmax, ymin, ymax):
        if xmin > xmax or ymin > ymax:
            print([xmin, xmax, ymin, ymax])
            raise ValueError("mins must be less than maxes")
        self._xmin = xmin
        self._xmax = xmax
        self._ymin = ymin
        self._ymax = ymax

    def intersect(self, r):
        xmin = max(self.xmin, r.xmin)
        xmax = min(self.xmax, r.xmax)
        ymin = max(self.ymin, r.ymin)
        ymax = min(self.ymax, r.ymax)
        if xmin >= xmax or ymin >= ymax:
            return None
        else:
            return Extent(xmin, xmax, ymin, ymax)

    def contains(self, x, y):
        if numpy.isscalar(x) and numpy.isscalar(y):
            return (x >= self.xmin and x <= self.xmax and y >= self.ymin and y <= self.ymax)
        else:
            xin = numpy.logical_and(x >= self.xmin, x <= self.xmax)
            yin = numpy.logical_and(y >= self.ymin, y <= self.ymax)
            return numpy.logical_and(xin, yin)

    def overlaps(self, e):
        if (e.xmin <= self.xmax and e.xmin >= self.xmin) or (e.xmax <= self.xmin and e.xmax >= self.xmin) or (
                self.xmin <= e.xmax and self.xmin >= e.xmin) or (self.xmax <= e.xmax and self.xmax >= e.xmin):
            if (e.ymin <= self.ymax and e.ymin >= self.ymin) or (e.ymax <= self.ymin and e.ymax >= self.ymin) or (
                    self.ymin <= e.ymax and self.ymin >= e.ymin) or (self.ymax <= e.ymax and self.ymax >= e.ymin):
                return True
        return False

    def __str__(self):
        return "xmin: " + str(self.xmin) + " xmax: " + str(self.xmax) + " ymin: " + str(self.ymin) + " ymax: " + str(
            self.ymax)

    @property
    def xmin(self):
        return self._xmin

    @xmin.setter
    def xmin(self, x):
        if x >= self.xmax:
            raise ValueError("xmin must be strictly less than xmax")
        self._xmin = x

    @property
    def xmax(self):
        return self._xmax

    @xmax.setter
    def xmax(self, x):
        if x <= self.xmin:
            raise ValueError("xmax must be strictly greater than xmin")
        self._xmax = x

    @property
    def ymin(self):
        return self._ymin

    @ymin.setter
    def ymin(self, y):
        if y >= self.ymax:
            raise ValueError("ymin must be strictly less than ymax")
        self._ymin = y

    @property
    def ymax(self):
        return self._ymax

    @ymax.setter
    def ymax(self, y):
        if y <= self.ymin:
            raise ValueError("ymax must be strictly greater than ymin")
        self._ymax = y


class Alignment(Extent):

    def __init__(self, e, res, ncol, nrow, projection):
        super().__init__(e.xmin, e.xmax, e.ymin, e.ymax)
        if numpy.isscalar(res):
            self._xres = res
            self._yres = res
        else:
            self._xres = res[0]
            self._yres = res[1]
        self._ncol = ncol
        self._nrow = nrow
        self._projection = projection

    # returns boolean 2-tuple; first value is whether we're dealing with scalars
    # 2nd is whether the coords are inbounds
    def _checkcoords(self, row=None, col=None, x=None, y=None, cell=None):
        rowcol = x is None and y is None and cell is None
        xy = row is None and col is None and cell is None
        cellused = row is None and col is None and x is None and y is None
        scalar = True
        inbounds = True
        if rowcol + xy + cellused > 1:
            raise ValueError("multiple of rowcol, xy, and cell passed to function")
        if not (rowcol or xy) and cell is None:
            raise ValueError("x, y, row, col and cell all None")
        if rowcol:
            if row is not None and col is not None:
                if numpy.isscalar(row) ^ numpy.isscalar(col):
                    raise ValueError("one of rowcol is scalar, other is not")
            if row is not None:
                if numpy.isscalar(row):
                    scalar = True
                    if row < 0 or row >= self.nrow:
                        inbounds = False
                else:
                    scalar = False
                    if numpy.min(row) < 0 or numpy.max(row) >= self.nrow:
                        inbounds = False
            if col is not None:
                if numpy.isscalar(col):
                    scalar = True
                    if col < 0 or col >= self.ncol:
                        inbounds = False
                else:
                    scalar = False
                    if numpy.min(col) < 0 or numpy.max(col) >= self.ncol:
                        inbounds = False
        if xy:
            if x is not None and y is not None:
                if numpy.isscalar(x) ^ numpy.isscalar(y):
                    raise ValueError("one of xy is scalar, other is not")
            if x is not None:
                if numpy.isscalar(x):
                    scalar = True
                    if x < self.xmin or x > self.xmax:
                        inbounds = False
                else:
                    scalar = False
                    if numpy.min(x) < self.xmin or numpy.max(x) > self.xmax:
                        inbounds = False
            if y is not None:
                if numpy.isscalar(y):
                    scalar = True
                    if y < self.ymin or y > self.ymax:
                        inbounds = False
                else:
                    scalar = False
                    if numpy.min(y) < self.ymin or numpy.max(y) > self.ymax:
                        inbounds = False
        if cellused:
            if numpy.isscalar(cell):
                scalar = True
                if cell < 0 or cell >= self.ncell:
                    inbounds = False
            else:
                scalar = False
                if numpy.min(cell) < 0 or numpy.max(cell) >= self.ncell:
                    inbounds = False
        return scalar, inbounds

    @staticmethod
    def _intconvert(n):
        if numpy.isscalar(n):
            return int(n)
        return n.astype(numpy.int32)

    def cellFromRowCol(self, row, col, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(row=row, col=col)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("rows or columns fall outside of raster")
        return col + self.ncol * row

    def colFromX(self, x, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(x=x)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("x falls outside of raster")
        return self._intconvert((x - self.xmin) / self.xres)

    def rowFromY(self, y, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(y=y)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("y falls outside of raster")
        return self._intconvert((self.ymax - y) / self.yres)

    def cellFromXY(self, x, y, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(x=x, y=y)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("x or y falls outside of raster")
        return self.cellFromRowCol(self.rowFromY(y, True), self.colFromX(x, True), safe)

    def cellsFromExtent(self, e):
        cells = numpy.arange(self.ncell)
        x, y = self.xyFromCell(cells, True)
        wxin = numpy.logical_and(x >= e.xmin, x <= e.xmax)
        wyin = numpy.logical_and(y >= e.ymin, y <= e.ymax)
        win = numpy.logical_and(wyin, wxin)
        return cells[win]

    def colFromCell(self, cell, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(cell=cell)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("cells fall outside of raster")
        return cell % self.ncol

    def rowFromCell(self, cell, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(cell=cell)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("cells fall outside of raster")
        return cell // self.ncol

    def rowColFromCell(self, cell, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(cell=cell)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("cells fall outside of raster")
        return (self.rowFromCell(cell, True), self.colFromCell(cell, True))

    def xFromCol(self, col, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(col=col)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("cols fall outside of raster")
        return self.xmin + self.xres * col + self.xres / 2

    def yFromRow(self, row, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(row=row)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("rows fall outside of raster")
        return self.ymax - self.yres * row - self.yres / 2

    def xyFromCell(self, cell, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(cell=cell)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("cells fall outside of raster")
        return (self.xFromCol(self.colFromCell(cell, True), True), self.yFromRow(self.rowFromCell(cell, True), True))

    def xFromCell(self, cell, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(cell=cell)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("cells fall outside of raster")
        return self.xFromCol(self.colFromCell(cell, True), True)

    def yFromCell(self, cell, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(cell=cell)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("cells fall outside of raster")
        return self.yFromRow(self.rowFromCell(cell, True), True)

    def alignExtent(self, e, snap="near"):
        if snap == "near":
            xmin = round((e.xmin - self.origin[0]) / self.xres) * self.xres + self.origin[0]
            xmax = round((e.xmax - self.origin[0]) / self.xres) * self.xres + self.origin[0]
            ymin = round((e.ymin - self.origin[1]) / self.yres) * self.yres + self.origin[1]
            ymax = round((e.ymax - self.origin[1]) / self.xres) * self.yres + self.origin[1]
        elif snap == "out":
            xmin = math.floor((e.xmin - self.origin[0]) / self.xres) * self.xres + self.origin[0]
            xmax = math.ceil((e.xmax - self.origin[0]) / self.xres) * self.xres + self.origin[0]
            ymin = math.floor((e.ymin - self.origin[1]) / self.yres) * self.yres + self.origin[1]
            ymax = math.ceil((e.ymax - self.origin[1]) / self.yres) * self.yres + self.origin[1]
        elif snap == "in":
            xmin = math.ceil((e.xmin - self.origin[0]) / self.xres) * self.xres + self.origin[0]
            xmax = math.floor((e.xmax - self.origin[0]) / self.xres) * self.xres + self.origin[0]
            ymin = math.ceil((e.ymin - self.origin[1]) / self.yres) * self.yres + self.origin[1]
            ymax = math.floor((e.ymax - self.origin[1]) / self.yres) * self.yres + self.origin[1]
        elif snap == "ll":
            xmin = math.floor((e.xmin - self.origin[0]) / self.xres) * self.xres + self.origin[0]
            xmax = math.floor((e.xmax - self.origin[0]) / self.xres) * self.xres + self.origin[0]
            ymin = math.floor((e.ymin - self.origin[1]) / self.yres) * self.yres + self.origin[1]
            ymax = math.floor((e.ymax - self.origin[1]) / self.yres) * self.yres + self.origin[1]
        if xmin == xmax:
            if xmin < e.xmin:
                xmax = xmax + self.xres
            else:
                xmin = xmin - self.xres
        if ymin == ymax:
            if ymin < e.ymin:
                ymax = ymax + self.yres
            else:
                ymin = ymin - self.yres
        return Extent(xmin, xmax, ymin, ymax)

    def indexExtent(self, e):
        if e.ymax > self.ymax:
            minrow = 0
        else:
            minrow = self.rowFromY(e.ymax)
        if minrow < 0 or not e.contains((e.xmin + e.xmax) / 2, self.yFromRow(minrow)):
            minrow += 1
        if e.ymin < self.ymin:
            maxrow = self.nrow - 1
        else:
            maxrow = self.rowFromY(e.ymin)
        if maxrow >= self.nrow or not e.contains((e.xmin + e.xmax) / 2, self.yFromRow(maxrow)):
            maxrow = maxrow - 1
        if e.xmin < self.xmin:
            mincol = 0
        else:
            mincol = self.colFromX(e.xmin)
        if mincol < 0 or not e.contains(self.xFromCol(mincol), (e.ymin + e.ymax) / 2):
            mincol += 1
        if e.xmax > self.xmax:
            maxcol = self.ncol - 1
        else:
            maxcol = self.colFromX(e.xmax)
        if maxcol >= self.ncol or not e.contains(self.xFromCol(maxcol), (e.ymin + e.ymax) / 2):
            maxcol = maxcol - 1
        return (minrow, maxrow, mincol, maxcol)

    @property
    def extent(self):
        return Extent(self._xmin, self._xmax, self._ymin, self._ymax)

    @property
    def res(self):
        if self._xres == self._yres:
            return self._xres
        else:
            return None

    @property
    def xres(self):
        return self._xres

    @property
    def yres(self):
        return self._yres

    @property
    def ncol(self):
        return self._ncol

    @property
    def nrow(self):
        return self._nrow

    @property
    def projection(self):
        return self._projection

    @projection.setter
    def projection(self, prj):
        self._projection = prj

    @property
    def xmin(self):
        return self._xmin

    @xmin.setter
    def xmin(self, x):
        raise AttributeError("Use functions such as crop or extend to change the properties of an Alignment or Raster")

    @property
    def xmax(self):
        return self._xmax

    @xmax.setter
    def xmax(self, x):
        raise AttributeError("Use functions such as crop or extend to change the properties of an Alignment or Raster")

    @property
    def ymin(self):
        return self._ymin

    @ymin.setter
    def ymin(self, y):
        raise AttributeError("Use functions such as crop or extend to change the properties of an Alignment or Raster")

    @property
    def ymax(self):
        return self._ymax

    @ymax.setter
    def ymax(self, y):
        raise AttributeError("Use functions such as crop or extend to change the properties of an Alignment or Raster")

    @property
    def ncell(self):
        return self._ncol * self._nrow

    @property
    def origin(self):
        return (self.xmin % self.xres, self.ymin % self.yres)

    def crop(self, e, snap="out"):
        if self.intersect(e) is None:
            raise ValueError("Extents do not overlap")
        e = self.intersect(e)
        e = self.alignExtent(e, snap=snap)
        (rowmin, rowmax, colmin, colmax) = self.indexExtent(Extent(e.xmin, e.xmax, e.ymin, e.ymax))
        self._xmin = e.xmin
        self._xmax = e.xmax
        self._ymin = e.ymin
        self._ymax = e.ymax
        self._nrow = rowmax - rowmin + 1
        self._ncol = colmax - colmin + 1

    @property
    def empty(self):
        return True

    def extend(self, e):
        if e.xmin >= self.xmin and e.xmax <= self.xmax and e.ymin >= self.ymin and e.ymax <= self.ymax:
            return
        olde = self.extent
        precols = 0
        postcols = 0
        prerows = 0
        postrows = 0
        if e.xmin < self.xmin:
            precols = -1 * self.colFromX(e.xmin, True)
            self._xmin = e.xmin
        if e.xmax > self.xmax:
            postcols = self.colFromX(e.xmax, True) - self.ncol + 1
            self._xmax = e.xmax
        if e.ymin < self.ymin:
            postrows = self.rowFromY(e.ymin, True) - self.nrow + 1
            self._ymin = e.ymin
        if e.ymax > self.ymax:
            prerows = -1 * self.rowFromY(e.ymax, True)
            self._ymax = e.ymax
        self._nrow = self.nrow + postrows + prerows
        self._ncol = self.ncol + precols + postcols


class Raster(Alignment):
    def __init__(self, e, res, ncol, nrow, projection, empty=False, data=None, dataType=None, NAvalue=-9999999):
        super().__init__(e, res, ncol, nrow, projection)
        self._dataType = dataType
        self._empty = empty
        if not empty:
            self._NAvalue = NAvalue
            self._dataType = dataType
            self.values = data

    def __getitem__(self, item):
        return self.values[item]

    # if both the key and value are scalars, the key is interpreted as a cell number
    # otherwise, the value needs to be something numpy-ish and the key needs to correspond properly
    # in either case,
    def __setitem__(self, key, value):
        if numpy.isscalar(key) and numpy.isscalar(value):
            row = self.rowFromCell(key)
            col = self.colFromCell(key)
            if numpy.isinf(value) or value == self.NAvalue:
                value = numpy.nan
            self._values[row][col] = value
        elif numpy.isscalar(value):
            if numpy.isinf(value) or value == self.NAvalue:
                value = numpy.nan
            self._values[key] = value
        else:
            value[numpy.isinf(value)] = numpy.nan
            value[value == self.NAvalue] = numpy.nan
            self._values[key] = value

    def __add__(self, other):
        return rasterWithNewValues(self, _saferastermathloop(self, other, operator.add), self.dataType)

    def __sub__(self, other):
        return rasterWithNewValues(self, _saferastermathloop(self, other, operator.sub), self.dataType)

    def __mul__(self, other):
        return rasterWithNewValues(self, _saferastermathloop(self, other, operator.mul), self.dataType)

    def __floordiv__(self, other):
        return rasterWithNewValues(self, _saferastermathloop(self, other, operator.floordiv),
                                   self.dataType)

    def __mod__(self, other):
        return rasterWithNewValues(self, _saferastermathloop(self, other, operator.mod), self.dataType)

    def __neg__(self):
        return rasterWithNewValues(self, numpy.mul(self.values, -1), self.dataType)

    def __truediv__(self, other):
        return rasterWithNewValues(self, _saferastermathloop(self, other, operator.truediv),
                                   self.dataType)

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return rasterWithNewValues(self, _saferastermathloop(self, other, operator.add, True),
                                   self.dataType)

    def __rmul__(self, other):
        return self * other

    def __rtruediv__(self, other):
        return rasterWithNewValues(self, _saferastermathloop(self, other, operator.truediv, True),
                                   self.dataType)

    def __rfloordiv__(self, other):
        return rasterWithNewValues(self, _saferastermathloop(self, other, operator.floordiv, True),
                                   self.dataType)

    def __iadd__(self, other):
        self._values = _saferastermathloop(self, other, operator.add)
        return self

    def __ifloordiv__(self, other):
        self._values = _saferastermathloop(self, other, operator.floordiv)

    def __imod__(self, other):
        self._values = _saferastermathloop(self, other, operator.mod)
        return self

    def __imul__(self, other):
        self._values = _saferastermathloop(self, other, operator.mul)
        return self

    def __isub__(self, other):
        self._values = _saferastermathloop(self, other, operator.sub)

    def __itruediv__(self, other):
        self._values = _saferastermathloop(self, other, operator.truediv)
        return self

    def __lt__(self, other):
        return self.values < other

    def __le__(self, other):
        return self.values <= other

    def __eq__(self, other):
        return self.values == other

    def __ne__(self, other):
        return self.values != other

    def __gt__(self, other):
        return self.values > other

    def __ge__(self, other):
        return self.values >= other

    def __len__(self):
        if self.empty:
            return 0
        else:
            return self.ncell

    def crop(self, e, snap="out"):
        if self.intersect(e) is None:
            raise ValueError("Extents do not overlap")
        e = self.intersect(e)
        e = self.alignExtent(e, snap=snap)
        usecells = self.cellsFromExtent(e)
        if not self.empty:
            v = self.valuesList()
            usevalues = v[usecells]
        super().crop(e)
        if not self.empty:
            self.values = usevalues

    @property
    def dataType(self):
        if self._empty:
            return None
        return self._dataType

    def extract(self, x, y, safe=False):
        if not safe:
            scalar, inbounds = self._checkcoords(x=x, y=y)
            if not inbounds:
                if scalar:
                    return None
                else:
                    raise ValueError("x or y falls outside of raster")
        return self.values[self.rowFromY(y, True), self.colFromX(x, True)]

    @property
    def values(self):
        if self._empty:
            return None
        return self._values

    @values.setter
    def values(self, v):
        if self._empty:
            raise RuntimeError("Empty rasters can't have values; use rasterWithNewValues()")
        if self.dataType == gdal.GDT_Float64:
            data = numpy.array(v, dtype=numpy.float64)
        else:
            data = numpy.array(v, dtype=numpy.float32)
        if len(data.shape) == 1 and data.shape[0] == self.ncell:
            self._values = data.reshape(self.nrow, self.ncol)
            self._values[self._values == self.NAvalue] = numpy.nan
        elif len(data.shape) == 2 and data.shape[0] == self.nrow and data.shape[1] == self.ncol:
            self._values = data
            self._values[self._values == self.NAvalue] = numpy.nan
        else:
            raise ValueError("Values for raster do not have the correct shape")

    def valuesList(self):
        return self._values.flatten()

    @property
    def NAvalue(self):
        if self._empty:
            return None
        return self._NAvalue

    @NAvalue.setter
    def NAvalue(self, n):
        if self._empty:
            raise RuntimeError("Empty rasters can't have values; use rasterWithNewValues()")
        if isinstance(n, numbers.Number):
            self.values.NAvalue = n
            self._NAvalue = n
        else:
            raise TypeError("NA Value must be a number")

    def writeRaster(self, f, dataType=None, driver="HFA"):
        if dataType is None:
            dataType = self.dataType
        d = gdal.GetDriverByName(driver)
        outr = d.Create(f, self.ncol, self.nrow, 1, dataType)
        outr.SetGeoTransform((self.xmin, self.xres, 0, self.ymax, 0, -1 * self.yres))
        outr.SetProjection(self.projection)
        if self.NAvalue is None:
            outr.GetRasterBand(1).SetNoDataValue(-9999)
        else:
            outr.GetRasterBand(1).SetNoDataValue(float(self.NAvalue))
        noNAs = numpy.array(self.values)
        noNAs[numpy.isnan(noNAs)] = self.NAvalue
        noNAs[numpy.isinf(noNAs)] = self.NAvalue
        outr.GetRasterBand(1).WriteArray(noNAs)
        del noNAs
        del outr

    def extend(self, e):
        if e.xmin >= self.xmin and e.xmax <= self.xmax and e.ymin >= self.ymin and e.ymax <= self.ymax:
            return
        olde = self.extent
        super().extend(e)
        if not self.empty:
            v = self.values
            newv = numpy.full(self.nrow, numpy.nan)
            newv[self.cellsFromExtent(olde, True)] = v
            self.values = newv

    @property
    def empty(self):
        return self._empty

    def trim(self):
        NAs = numpy.logical_not(numpy.isnan(self.values))
        cols = numpy.apply_along_axis(numpy.sum, 0, NAs)
        rows = numpy.apply_along_axis(numpy.sum, 1, NAs)
        wcol = numpy.where(cols > 0)[0]
        wrow = numpy.where(rows > 0)[0]
        mincol = wcol[0]
        maxcol = wcol[-1]
        minrow = wrow[0]
        maxrow = wrow[-1]
        e = Extent(_myround(self.xFromCol(mincol) - self.xres / 2, self.xres, self.origin[0]),
                   _myround(self.xFromCol(maxcol) + self.xres / 2, self.xres, self.origin[0]),
                   _myround(self.yFromRow(maxrow) - self.yres / 2, self.yres, self.origin[1]),
                   _myround(self.yFromRow(minrow) + self.yres / 2, self.yres, self.origin[1]))
        self.crop(e)

    def maskByRaster(self, m):
        if not (self.xres == m.xres and self.yres == m.yres and self.nrow == m.nrow and
                self.ncol == m.ncol and numpy.allclose(self.origin, m.origin,
                                                       rtol=0.05) and self.projection == m.projection):
            raise ValueError("Mask must exactly match this raster")
        self[numpy.isnan(m.values)] = numpy.nan

    def maskByPoly(self, m):
        try:
            polyraster = rasterio.features.rasterize([m['geometry']], out_shape=(self.nrow, self.ncol),
                                                     transform=self.affine,
                                                     fill=1, default_value=0).astype(numpy.bool_)
        except TypeError:
            polyraster = rasterio.features.rasterize([m], out_shape=(self.nrow, self.ncol), transform=self.affine,
                                                     fill=1, default_value=0).astype(numpy.bool_)
        self[polyraster] = numpy.nan

    # fun needs to take a vector and preferably treat nan properly
    def zonal(self, r, fun):
        d = {}
        e = r.extent
        c = self.cellsFromExtent(e)
        v = self.valuesList()[c]
        rv = r.extract(self.xFromCell(c, True), self.yFromCell(c, True))
        dtemp = {}
        for i in range(len(c)):
            if rv[i] is None:
                continue
            if numpy.isnan(rv[i]):
                continue
            if rv[i] not in dtemp:
                dtemp[rv[i]] = [v[i]]
            else:
                dtemp[rv[i]].append(v[i])
        for value in dtemp:
            d[value] = fun(dtemp[value])
        del dtemp
        return d

    def isna(self, v):
        return numpy.isinf(v) or numpy.isnan(v) or v == self.NAValue

    @property
    def geo(self):
        return self.xmin, self.xres, 0, self.ymax, 0, -self.yres

    @property
    def affine(self):
        return affine.Affine.from_gdal(*self.geo)


class RasterIndex:
    def __init__(self, flist, recursive=False, scale=1, shift=0, projection="", extent=None, verbose=False):

        self._customextent = False
        if isinstance(flist, type({})):
            self._data = copy.deepcopy(flist)
            self._projection = projection
            self._scale = scale
            self._shift = shift
            if extent is not None:
                self._extent = extent
            else:
                self._extent = self.__getExtent()
        else:
            self._data = {}
            self._projection = projection

            def rastersindir(d):
                for f in os.listdir(d):
                    if os.path.isdir(os.path.join(d, f)):
                        if recursive:
                            rastersindir(os.path.join(d, f))
                    try:
                        addraster(os.path.join(d, f))
                    except AttributeError:
                        pass

            def addraster(f):
                r = rasterFromFile(f, 1, empty=True)
                if verbose:
                    print(f)
                if self.projection == "":
                    if not r.projection == "":
                        self.projection = r.projection
                else:
                    if not r.projection == "" and not r.projection == self.projection:
                        raise ValueError("Rasters do not all have the same projection")
                self._data[f] = r

            if isinstance(flist, str):
                flist = [flist]
            for f in flist:
                if not os.access(f, os.F_OK):
                    raise FileNotFoundError(f)
                if os.path.isdir(f):
                    rastersindir(f)
                else:
                    addraster(f)

            self._scale = scale
            self._shift = shift
            self._extent = self.__getExtent()

    @property
    def filelist(self):
        return [f for f in self._data.keys()]

    @property
    def rasterlist(self):
        return self._data.values()

    @property
    def data(self):
        return copy.copy(self._data)

    @property
    def projection(self):
        return self._projection

    @projection.setter
    def projection(self, prj):
        self._projection = prj

    def __getitem__(self, item):
        if item in self._data:
            r = rasterFromFile(item)
            if not self._scale == 1:
                r *= self._scale
            if not self._shift == 0:
                r += self._shift
            return r
        else:
            raise IndexError()

    def getempty(self, item):
        return self._data[item]

    def __add__(self, other):
        if numpy.isscalar(other):
            return RasterIndex(self._data, scale=self._scale, shift=self._shift + other, projection=self.projection,
                               extent=self._extent)
        if isinstance(other, type(self)):
            raise TypeError("math between raster indices not supported")
        return other.__radd__(self)

    def __sub__(self, other):
        if numpy.isscalar(other):
            return RasterIndex(self._data, scale=self._scale, shift=self._shift - other, projection=self.projection,
                               extent=self._extent)
        if isinstance(other, type(self)):
            raise TypeError("math between raster indices not supported")
        else:
            return -other.__add__(self)

    def __mul__(self, other):
        if numpy.isscalar(other):
            return RasterIndex(self._data, scale=self._scale * other, shift=self._shift * other,
                               projection=self.projection, extent=self._extent)
        if isinstance(other, type(self)):
            raise TypeError("math between raster indices not supported")
        else:
            return other.__rmul__(self)

    def __neg__(self):
        return RasterIndex(self._data, scale=-self._scale, shift=-self._shift, projection=self.projection,
                           extent=self._extent)

    def __truediv__(self, other):
        if numpy.isscalar(other):
            return RasterIndex(self._data, scale=self._scale / other, shift=self._shift / other,
                               projection=self.projection, extent=self._extent)
        if isinstance(other, type(self)):
            raise TypeError("math between raster indices not supported")
        else:
            return other.rtruediv(self)

    def __radd__(self, other):
        if numpy.isscalar(other):
            return RasterIndex(self._data, scale=self._scale, shift=self._shift + other,
                               projection=self.projection, extent=self._extent)
        if isinstance(other, type(self)):
            raise TypeError("math between raster indices not supported")
        else:
            return other.__add__(self)

    def __rsub__(self, other):
        if numpy.isscalar(other):
            return RasterIndex(self._data, scale=-self._scale, shift=-self._shift + other,
                               projection=self.projection, extent=self._extent)
        if isinstance(other, type(self)):
            raise TypeError("math between raster indices not supported")
        else:
            return other.__sub__(self)

    def __rmul__(self, other):
        if numpy.isscalar(other):
            return RasterIndex(self._data, scale=self._scale * other, shift=self._shift * other,
                               projection=self.projection, extent=self._extent)
        if isinstance(other, type(self)):
            raise TypeError("math between raster indices not supported")
        else:
            return other.__mul__(self)

    def __rtruediv__(self, other):
        if numpy.isscalar(other):
            return RasterIndex(self._data, scale=self._scale / other, shift=self._shift / other,
                               projection=self.projection, extent=self._extent)
        if isinstance(other, type(self)):
            raise TypeError("math between raster indices not supported")
        else:
            return other.__div__(self)

    def __iadd__(self, other):
        if numpy.isscalar(other):
            self._shift += other
            return self
        else:
            raise TypeError()

    def __imul__(self, other):
        if numpy.isscalar(other):
            self._shift *= other
            self._scale *= other
            return self
        else:
            raise TypeError()

    def __isub__(self, other):
        if numpy.isscalar(other):
            self._shift -= other
            return self
        else:
            raise TypeError()

    def __itruediv__(self, other):
        if numpy.isscalar(other):
            self._shift /= other
            self._scale /= other
            return self
        else:
            raise TypeError()

    def multiply(self, other):
        if numpy.isscalar(other):
            self._shift *= other
            self._scale *= other
        else:
            raise TypeError

    def add(self, other):
        if numpy.isscalar(other):
            self._shift *= other
            self._scale *= other
        else:
            raise TypeError

    @property
    def extent(self):
        return self._extent

    @extent.setter
    def extent(self, e):
        self._extent = Extent(e.xmin, e.xmax, e.ymin, e.ymax)
        self._customextent = True

    def removeRaster(self, f):
        r = self._data.pop(f)
        if not self._customextent:
            self._extent = self.__getExtent()

    def __getExtent(self):
        xmin = numpy.inf
        ymin = numpy.inf
        xmax = numpy.NINF
        ymax = numpy.NINF
        for r in self._data.values():
            if r.xmin < xmin:
                xmin = r.xmin
            if r.xmax > xmax:
                xmax = r.xmax
            if r.ymin < ymin:
                ymin = r.ymin
            if r.ymax > ymax:
                ymax = r.ymax
        self._customextent = True
        return Extent(xmin, xmax, ymin, ymax)

    def rasterListFromXY(self, x, y, empty=True):
        if not self.extent.contains(x, y):
            return None
        out = []
        for f in self.filelist:
            if self.getempty(f).extent.contains(x, y):
                if empty:
                    out.append(self.getempty(f))
                else:
                    r = self[f]
                    if r.extract(x, y) is not None:
                        out.append(self[f])
        return out

    def bestRasterFromXY(self, x, y, empty=True):
        if not self.extent.contains(x, y):
            return None
        res = numpy.inf
        out = None
        outf = None
        for f in self.filelist:
            emp = self.getempty(f)
            if emp.res < res and emp.extent.contains(x, y):
                if not empty:
                    r = self[f]
                    if r.extract(x, y) is not None:
                        out = r
                        outf = f
                        res = r.res
                else:
                    out = emp
                    res = emp.res
                    outf = f
        return outf, out

    def extract(self, x, y):
        _, r = self.bestRasterFromXY(x, y)
        if r is None:
            return None
        return r.extract(x, y)

    @property
    def xmin(self):
        return self.extent.xmin

    @xmin.setter
    def xmin(self, x):
        self.extent = Extent(x, self.xmax, self.ymin, self.ymax)

    @property
    def xmax(self):
        return self.extent.xmax

    @xmax.setter
    def xmax(self, x):
        self.extent = Extent(self.xmin, x, self.ymin, self.ymax)

    @property
    def ymin(self):
        return self.extent.ymin

    @ymin.setter
    def ymin(self, y):
        self.extent = Extent(self.xmin, self.xmax, y, self.ymax)

    @property
    def ymax(self):
        return self.extent.ymax

    @ymax.setter
    def ymax(self, y):
        self.extent = Extent(self.xmin, self.xmax, self.ymin, y)

    @property
    def scale(self):
        return self._scale

    @scale.setter
    def scale(self, s):
        self._scale = s

    @property
    def shift(self):
        return self._shift

    @shift.setter
    def shift(self, s):
        self._shift = s
