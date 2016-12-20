import numpy as np
import netCDF4 as nc
import os
from pyproj import Proj

def getnemo(filename,fillmissing=False):
    # Loads bathy from a NEMO netcdf file
    sspath='/home/mdunphy/MEOPAR/NEMO-forcing/grid'
    bnc = nc.Dataset(filename, 'r')
    bathy = bnc.variables['Bathymetry'][:]
    bnc.close()
    if fillmissing is True:
        bathy = bathy.filled(0)
    return bathy


def getchs(filename):
    # Load CHS data, cache in a .npz file for faster loading
    cache=filename+".npz"
    if os.path.exists(cache):
        npzfile = np.load(cache)
        x, y, z = npzfile['x'], npzfile['y'], npzfile['z']
    else:
        a = np.loadtxt(filename, delimiter=',')
        x, y, z = np.float32(a[:,0]), np.float32(a[:,1]), a[:,2]
        np.savez(cache, x=x, y=y, z=z)
    p = Proj("+proj=utm +zone=10U")    # UTM to lon, lat using zone 10U
    return x,y,z,p


def getcascadia(filename):
    # Inspired by: https://pymorton.wordpress.com/2016/02/26/plotting-prism-bil-arrays-without-using-gdal/
    nrows, ncols = 6435, 5951
    z = np.fromfile(filename, dtype=np.int16).byteswap()   # load data, 16 bit big-endian integers
    z = z.reshape(nrows,ncols)                 # reshape
    z = np.flipud(z )                          # flip because the data is stored from north to south
    mask = (z == 0) | (z >= 10000)             # mask for nonexistant points and land points
    z -= 10000                                 # remove offset
    z *= -1                                    # make depths positive
    z[mask] = 0                                # set masked values to zero
    # Construct Cascadia coordinates
    xmin, xmax, dx = -738044.062, 749705.938, 250
    ymin, ymax, dy = 101590.289, 1710340.289, 250
    x = xmin + dx*np.arange(0, z.shape[1]) + dx/2
    y = ymin + dy*np.arange(0, z.shape[0]) + dy/2
    p = Proj(r'+proj=lcc +lat_1=41.5 +lat_2=50.5 +lat_0=38 +lon_0=-124.5 +x_0=0 +y_0=0 +ellps=clrk66 +no_defs')
    return x,y,z,p


def getbc3(filename):
    # Load BC3 data, cache depths in a .npz file for faster loading
    cache=filename+".npz"
    if os.path.exists(cache):
        npzfile = np.load(cache)
        z = npzfile['z']
    else:
        z = np.loadtxt(filename, dtype='float32', skiprows=6)
        np.savez(cache, z=z)
    z = np.flipud(z)              # flip because the data is stored from north to south
    mask = (z == -9999) | (z > 0) # mask for nonexistant points and land points
    z *= -1                       # make depth positive
    z[mask] = 0                   # set masked values to zero
    # Documentation: https://www.ngdc.noaa.gov/dem/squareCellGrid/getReport/4956
    # As of Dec 19, 2016, there is a typo in the metadata, which claims
    # that ymin is 48.50 N, but ymin is actually 48.05 N
    # The headers in the .asc file are:
    #    ncols         18301
    #    nrows         7381
    #    xllcorner     -137.45041666368
    #    yllcorner     48.049583333335
    #    cellsize      0.000833333333
    #    NODATA_value  -9999
    # We recover (xllcorner,yllcorner) by computing (xmin-d/2, ymin-d/2)
    xmin, ymin, d = -137.45, 48.05, 3/3600
    x = xmin + d*np.arange(0, z.shape[1])
    y = ymin + d*np.arange(0, z.shape[0])
    p = None
    return x,y,z,p

