import rasterio
import numpy as np
import math

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
        if np.isscalar(x) and np.isscalar(y):
            return (x >= self.xmin and x <= self.xmax and y >= self.ymin and y <= self.ymax)
        else:
            xin = np.logical_and(x >= self.xmin, x <= self.xmax)
            yin = np.logical_and(y >= self.ymin, y <= self.ymax)
            return np.logical_and(xin, yin)

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
        if np.isscalar(res):
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
                if np.isscalar(row) ^ np.isscalar(col):
                    raise ValueError("one of rowcol is scalar, other is not")
            if row is not None:
                if np.isscalar(row):
                    scalar = True
                    if row < 0 or row >= self.nrow:
                        inbounds = False
                else:
                    scalar = False
                    if np.min(row) < 0 or np.max(row) >= self.nrow:
                        inbounds = False
            if col is not None:
                if np.isscalar(col):
                    scalar = True
                    if col < 0 or col >= self.ncol:
                        inbounds = False
                else:
                    scalar = False
                    if np.min(col) < 0 or np.max(col) >= self.ncol:
                        inbounds = False
        if xy:
            if x is not None and y is not None:
                if np.isscalar(x) ^ np.isscalar(y):
                    raise ValueError("one of xy is scalar, other is not")
            if x is not None:
                if np.isscalar(x):
                    scalar = True
                    if x < self.xmin or x > self.xmax:
                        inbounds = False
                else:
                    scalar = False
                    if np.min(x) < self.xmin or np.max(x) > self.xmax:
                        inbounds = False
            if y is not None:
                if np.isscalar(y):
                    scalar = True
                    if y < self.ymin or y > self.ymax:
                        inbounds = False
                else:
                    scalar = False
                    if np.min(y) < self.ymin or np.max(y) > self.ymax:
                        inbounds = False
        if cellused:
            if np.isscalar(cell):
                scalar = True
                if cell < 0 or cell >= self.ncell:
                    inbounds = False
            else:
                scalar = False
                if np.min(cell) < 0 or np.max(cell) >= self.ncell:
                    inbounds = False
        return scalar, inbounds

    @staticmethod
    def _intconvert(n):
        if np.isscalar(n):
            return int(n)
        return n.astype(np.int32)

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
        cells = np.arange(self.ncell)
        x, y = self.xyFromCell(cells, True)
        wxin = np.logical_and(x >= e.xmin, x <= e.xmax)
        wyin = np.logical_and(y >= e.ymin, y <= e.ymax)
        win = np.logical_and(wyin, wxin)
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
    def __init__(self, e, res, ncol, nrow, projection, data, dataType, NAvalue=-9999999):
        super().__init__(e, res, ncol, nrow, projection)
        self.NAvalue = NAvalue
        self.values = data
        self.dtype = dataType

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, v):
        self._values = np.array(v, dtype=np.float64)
        self._values[self._values == self.NAvalue] = np.nan

    def writeRaster(self, filename):
        if len(self.values.shape) > 2:
            bands = self.values.shape[0]
        else:
            bands = 1
        transform = rasterio.Affine(self.xres, 0, self.xmin, 0, -self.yres, self.ymax)
        with rasterio.open(
                filename,
                'w',
                driver='GTiff',
                height=self.nrow,
                width=self.ncol,
                count=bands,
                dtype=self.dtype,
                crs=rasterio.CRS.from_wkt(self.projection),
                transform=transform,
        ) as dst:
            if bands == 1:
                dst.write(self.values, 1)
            else:
                dst.write(self.values)
