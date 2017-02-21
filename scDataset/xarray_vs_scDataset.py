import os,sys,resource,fnmatch,time
from scDataset import scDataset
from IPython import embed
import xarray as xr
import numpy as np

# Run this as:
#    python xarray_vs_scDataset.py scDataset
# or
#    python xarray_vs_scDataset.py xarray

def gethourlyfiles(prefix, days, grid):
    files = []    
    for day in days:
        dirname = os.path.join(prefix, day)
        for item in os.listdir(dirname):
            if fnmatch.fnmatchcase(item, "SalishSea_1h_*"+grid+"*.nc"):
                files += [os.path.join(dirname,item)]
    return files

# Mar 17 is a netcdf discontinuity in nowcast-green
prefix = '/results/SalishSea/nowcast-green'
days = ['{:02d}mar16'.format(x) for x in range(5,12)]
files = gethourlyfiles(prefix, days, 'grid_T')
print("Found {} files".format(len(files)))

if sys.argv[1] == "scDataset":
    t0 = time.time()
    ds = scDataset(files)
    t  = ds.variables['votemper'][:,:,500,250]
    ds.close()
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print("scDataset        , {:.2f}s, Peak Memory Usage {:.2f} MB".format(time.time()-t0, mem/1024.0))

if sys.argv[1] == "xarray":
    t0 = time.time()
    ds = xr.open_mfdataset(files)
    t  = ds.variables['votemper'][:,:,500,250].values
    ds.close()
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print("xr.open_mfdataset, {:.2f}s, Peak Memory Usage {:.2f} MB".format(time.time()-t0, mem/1024.0))

# Benchmark results at salish:
#
# Found 7 files
# scDataset        , 27.68s, Peak Memory Usage 227.39 MB
#
# Found 7 files
# xr.open_mfdataset, 39.98s, Peak Memory Usage 19828.78 MB
