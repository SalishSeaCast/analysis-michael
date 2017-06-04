import netCDF4 as nc
import numpy as np
import sys

def fix_bathy(infile, mindep):
    """
    Simple script to enforce minimum depth and fix the longitudes
    """
    with nc.Dataset(infile, 'r+') as f:
        # Enforce minimum bathymetry
        bm = f.variables['Bathymetry'][:]
        idx = (bm > 0) & (bm < mindep)
        if np.any(idx):
            md = np.min(bm[idx])
            print("Min depth {:3f} m, resetting to {:3f} m".format(md, mindep))
            bm[idx] = mindep
            f.variables['Bathymetry'][:] = bm
            
        # Enforce nav_lon to be in [-180,180] and not [0,360]
        lon = f.variables['nav_lon'][:]
        if np.any(lon > 180):
            lon[lon>180] -= 360
        f.variables['nav_lon'][:] = lon
        f.variables['nav_lon'].valid_min = np.min(lon)
        f.variables['nav_lon'].valid_max = np.max(lon)
    
if __name__ == "__main__":
    fix_bathy(sys.argv[1], 4)
