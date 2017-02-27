import os, fnmatch
import netCDF4 as nc
import numpy as np
from time import time
from scDataset import scDataset
#from IPython import embed


def gethourlyfiles(prefix, days, grid):
    files = []
    for day in days:
        dirname = os.path.join(prefix, day)
        for item in os.listdir(dirname):
            if fnmatch.fnmatchcase(item, "SalishSea_1h_*" + grid + "*.nc"):
                files += [os.path.join(dirname, item)]
    return files


t0 = time()
prefix = '/results/SalishSea/nowcast-blue'
days = ['{:02d}dec16'.format(x) for x in range(1, 4 + 1)]
filesT = gethourlyfiles(prefix, days, 'grid_T')
filesU = gethourlyfiles(prefix, days, 'grid_U')
filesV = gethourlyfiles(prefix, days, 'grid_V')
print("Found {} files [{:.3f} s]".format(len(filesT), time() - t0))


print("\n\nTest 1 (load entire time_counter):")
t0 = time()
tc1 = np.array([])
for f in filesT:
    ds = nc.Dataset(f, 'r')
    tc1 = np.append(tc1, ds.variables['time_counter'][:])
    ds.close()
print("Loop on nc.Dataset [{:.3f} s]".format(time() - t0))

t0 = time()
ds = scDataset(filesT)
tc2 = ds.variables['time_counter'][:]
ds.close()
print("scDataset (no context manager), max diff {:g}, [{:.3f} s]".format(np.max(np.abs(tc1 - tc2)), time() - t0))

t0 = time()
with scDataset(filesT) as ds:
    tc3 = ds.variables['time_counter'][:]
print("scDataset, max diff {:g}, [{:.3f} s]".format(np.max(np.abs(tc1 - tc3)), time() - t0))


print("\n\nTest 2 (temperature at all t, z)")
t0 = time()
tz1 = None
for f in filesT:
    ds = nc.Dataset(f, 'r')
    data = ds.variables['votemper'][:, :, 300, 200]
    if tz1 is None:
        tz1 = data
    else:
        tz1 = np.append(tz1, data, 0)
    ds.close()
print("Loop on nc.Dataset [{:.3f} s]".format(time() - t0))

t0 = time()
with scDataset(filesT) as ds:
    tz2 = ds.variables['votemper'][:, :, 300, 200]
print("scDataset, max diff {:g}, [{:.3f} s]".format(np.max(np.abs(tz1 - tz2)), time() - t0))


print("\n\nTest 3 (max surface salinity at all t)")
t0 = time()
mx1 = np.array([])
for f in filesT:
    ds = nc.Dataset(f, 'r')
    for ti in range(ds.variables['vosaline'].shape[0]):
        sss = ds.variables['vosaline'][ti, 0, :, :]
        mx1 = np.append(mx1, np.max(sss))
    ds.close()
print("Loop on nc.Dataset [{:.3f} s]".format(time() - t0))

t0 = time()
with scDataset(filesT) as ds:
    mx2 = np.array([])
    for ti in range(ds.variables['vosaline'].shape[0]):
        sss = ds.variables['vosaline'][ti, 0, :, :]
        mx2 = np.append(mx2, np.max(sss))
print("scDataset, max diff {:g}, [{:.3f} s]".format(np.max(np.abs(mx1 - mx2)), time() - t0))


print("\n\nTest 4 (slicing test)")
with scDataset(filesT) as ds:
    tx = ds.variables['votemper'][20:30, 0:4, -20:, ::2]
    expectation = (10, 4, 20, 199)
    print("Expected shape:", expectation, "\nActual shape:  ", tx.shape)


print("\n\nTest 5 (open two context managers)")
t0 = time()
with scDataset(filesU) as dsu, scDataset(filesV) as dsv:
    tu = dsu.variables['time_counter'][:]
    tv = dsv.variables['time_counter'][:]
print("[{:.3f} s]".format(np.max(np.abs(time() - t0))))


print("\n\nTest 6 (we expect to get an exception for mismatched dimensions)")
t0 = time()
with scDataset(filesT) as ds:
    try:
        t = ds.variables['time_counter'][:, :, :]
    except ValueError as err:
        print("Caught ValueError: {}".format(err))
print("[{:.3f} s]".format(np.max(np.abs(time() - t0))))
