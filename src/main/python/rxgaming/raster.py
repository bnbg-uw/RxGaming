import rasterio
import numpy as np
class Extent:
    def __init__(self, xmin, xmax, ymin, ymax):
        if xmin > xmax or ymin > ymax:
            print([xmin, xmax, ymin, ymax])
            raise ValueError("mins must be less than maxes")
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def __str__(self):
        return "xmin: " + str(self.xmin) + " xmax: " + str(self.xmax) + " ymin: " + str(self.ymin) + " ymax: " + str(
            self.ymax)

class Raster:
    def __init__(self, e, res, ncol, nrow, projection, data, dataType, NAvalue=-9999999):
        self.NAvalue = NAvalue
        self.values = data
        self.extent = e
        if len(res) == 1:
            self.res = [res, res]
        else:
            self.res = res
        self.ncol = ncol
        self.nrow = nrow
        self.projection = projection
        self.dtype = dataType

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, v):
        self._values = np.array(v, dtype=np.float64)
        self._values[self._values == self.NAvalue] = np.nan

    @property
    def xres(self):
        return self.res[0]

    @property
    def yres(self):
        return self.res[1]

    @property
    def xmin(self):
        return self.extent.xmin

    @property
    def xmax(self):
        return self.extent.xmax

    @property
    def ymin(self):
        return self.extent.ymin

    @property
    def ymax(self):
        return self.extent.ymax
    def writeRaster(self, filename):
        transform = rasterio.Affine(self.res[0], 0, self.extent.xmin, 0, -self.res[1], self.extent.ymax)
        with rasterio.open(
                filename,
            'w',
            driver='GTiff',
            height=self.nrow,
            width=self.ncol,
            count=1,
            dtype=self.dtype,
            crs=rasterio.CRS.from_wkt(self.projection),
            transform=transform,
        ) as dst:
            dst.write(self.values, 1)
