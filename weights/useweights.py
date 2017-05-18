from salishsea_tools import grid_tools

from IPython import embed
import netCDF4 as nc
import scipy.interpolate as spi
import scipy.sparse as sp
import numpy as np

# Weights
weightsfile = '/home/mdunphy/MEOPAR/NEMO-forcing/grid/weights-gem2.5-ops.nc'
with nc.Dataset(weightsfile) as f:
    s1 = f.variables['src01'][:]-1  # minus one for fortran-to-python indexing
    s2 = f.variables['src02'][:]-1
    s3 = f.variables['src03'][:]-1
    s4 = f.variables['src04'][:]-1
    w1 = f.variables['wgt01'][:]
    w2 = f.variables['wgt02'][:]
    w3 = f.variables['wgt03'][:]
    w4 = f.variables['wgt04'][:]

# Operational data
opsfile='/results/forcing/atmospheric/GEM2.5/operational/ops_y2017m04d29.nc'
with nc.Dataset(opsfile) as f:
    odata = f.variables['tair'][0,...]   # Load a 2D field


NO = odata.size   # number of operational grid points
NN = s1.size      # number of NEMO grid points

# Build matrix
n = np.array([x for x in range(0,NN)])
M1 = sp.csr_matrix((w1.flatten(), (n, s1.flatten())), (NN,NO))
M2 = sp.csr_matrix((w2.flatten(), (n, s2.flatten())), (NN,NO))
M3 = sp.csr_matrix((w3.flatten(), (n, s3.flatten())), (NN,NO))
M4 = sp.csr_matrix((w4.flatten(), (n, s4.flatten())), (NN,NO))
M = M1+M2+M3+M4

# Interpolate by matrix multiply - quite fast
ndata = M*odata.flatten()

# Reshape to NEMO shaped array
ndata=ndata.reshape(s1.shape)



# We can /approximately/ verify by griddata; NEMO's weights file uses spherical
# coordinates for the remapping, while griddata does not, so expect some differences.

# Get the NEMO T grid points
coordsfile = '/home/mdunphy/MEOPAR/NEMO-forcing/grid/coordinates_seagrid_SalishSea.nc'
with nc.Dataset(coordsfile) as f:
    nlon = f.variables['glamt'][:]
    nlat = f.variables['gphit'][:]

# Also the operational grid points 
with nc.Dataset(opsfile) as f:
    olon = f.variables['nav_lon'][:]-360
    olat = f.variables['nav_lat'][:]

# Call griddata
gdata=spi.griddata((olon.flatten(),olat.flatten()),odata.flatten(),(nlon.flatten(),nlat.flatten()))
gdata=gdata.reshape(s1.shape)



# Check that the grid_tools version is doing the same thing
matrix, nemosize = grid_tools.build_matrix(weightsfile, opsfile)
ndata2 = grid_tools.use_matrix(opsfile, matrix, nemosize, 'tair', 0)
print(np.max(np.abs(ndata-ndata2)))

embed()
