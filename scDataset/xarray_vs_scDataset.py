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

t0 = time.time()

if sys.argv[1] == "scDataset":
    ds = scDataset(files)
    t  = ds.variables['votemper'][:,:,500,250]
    ds.close()
if sys.argv[1] == "xarray":
    ds = xr.open_mfdataset(files)
    t  = ds.variables['votemper'][:,:,500,250].values
    ds.close()

mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print("{:10s}, {:.2f}s, Peak Memory Usage {:.2f} MB".format(sys.argv[1], time.time()-t0, mem/1024.0))

# Benchmark results at salish:
#
# Found 7 files
# scDataset , 27.58s, Peak Memory Usage 227.46 MB
# Found 7 files
# xarray    , 38.81s, Peak Memory Usage 19833.56 MB
