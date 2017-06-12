# This script runs a bunch of indexing tests against scDataset to
# make sure it's doing the Right Thing. The expected output is that the size
# matches on every test with maxerr = 0.

import netCDF4 as nc
import numpy as np

from salishsea_tools.nc_tools import scDataset
#import sys; sys.path.append('/home/mdunphy/MEOPAR/analysis-michael/scDataset')
#from scDataset import scDataset

NT,NZ,NY,NX = 22,24,26,28
f,f1,f2 = '/tmp/f.nc', '/tmp/f1.nc', '/tmp/f2.nc'

def opn(fx):
    dx = nc.Dataset(fx, 'w', clobber=True)
    dx.createDimension('x', NX)
    dx.createDimension('y', NY)
    dx.createDimension('z', NZ)
    dx.createDimension('t', None)
    dx.createVariable('x', 'f', ('x'))
    dx.createVariable('y', 'f', ('y'))
    dx.createVariable('z', 'f', ('z'))
    dx.createVariable('time', 'f', ('t'))
    dx.createVariable('salt', 'f', ('t','z','y','x'))
    dx.variables['x'] = [q for q in range(NX)]
    dx.variables['y'] = [q for q in range(NY)]
    dx.variables['z'] = [q for q in range(NZ)]
    dx.description = "test file: " + fx
    return dx

with opn(f) as d, opn(f1) as d1, opn(f2) as d2:
    for ti in range(NT):
        vals = np.random.rand(NZ,NY,NX)
        # Combined file
        d.variables['salt'][ti,...] = vals
        d.variables['time'][ti] = ti
        if ti < NT/2:
            # Part 1
            d1.variables['salt'][ti % (NT/2),...] = vals
            d1.variables['time'][ti % (NT/2)] = ti
        else:
            # Part 2
            d2.variables['salt'][ti % (NT/2),...] = vals
            d2.variables['time'][ti % (NT/2)] = ti

# Helper class for verifying the indexing
class checker(object):
    def __init__(self, a, b):
        self.a, self.b = a, b
    def __getitem__(self, slices):
        v1, v2 = self.a[slices], self.b[slices]
        if v1.shape == v2.shape:
            maxerr = np.max(np.abs(v1-v2))
            print("Size matches, max err = {}".format(maxerr))
        else:
            print("Size mismatch for indexing {}".format(slices))

# Open both datasets
with nc.Dataset(f) as ncds, scDataset([f1,f2]) as scds:
    c = checker(ncds.variables['salt'], scds.variables['salt'])
    
    # Finally, we run a slew of indexing tests
    c[:]
    c[...]
    c[:,...,:]
    c[...,:]
    c[:,...]
    c[0,0,0,0]
    c[0,0,0]
    c[0,0,0,:]
    c[0,0,0,...]
    c[0,0,:,...]
    c[0,...,:]
    c[0,:,...]
    c[:,:,:,:]
    c[4:15,2:8,3:7,4:12]
    c[::-1,::-2,::-3,::-4]
    c[::-1,...,::-3,:]
    c[:,4,5,6]
    c[:,4:6,3:9,::-2]
    c[...,:,:,:]
    c[...,4,:,8]
    c[6,6,6,6,...]
    c[...,6,6,6,6]
    c[1,2,:,...,:]
