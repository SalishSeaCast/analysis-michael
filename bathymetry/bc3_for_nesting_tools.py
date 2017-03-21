import numpy as np
import netCDF4 as nc
from bathy_helpers import *
from bathy_readers import *

with nc.Dataset('/home/mdunphy/MEOPAR/NEMO-forcing/grid/coordinates_seagrid_SalishSea201702.nc', 'r') as cnc:
    glamf = cnc.variables['glamf'][0,...]; gphif = cnc.variables['gphif'][0,...]
    glamt = cnc.variables['glamt'][0,...]; gphit = cnc.variables['gphit'][0,...]
glamfe, gphife = expandf(glamf, gphif)

# Create external elevation database for AGRIF nesting tools
outfile = "/home/mdunphy/MEOPAR/WORK/Bathy-201702/BC3/BC3_For_Nesting_Tools.nc"

# Source data
bc3file = '/home/mdunphy/MEOPAR/WORK/Bathy-201702/BC3/british_columbia_3sec.asc'
x,y,z,p = getbc3(bc3file)
    
# There are too many points, so we filter based on minimum and maximum x,y
idx1 = (x >= np.min(glamfe)) & (x <= np.max(glamfe))
idx2 = (y >= np.min(gphife)) & (y <= np.max(gphife))
x,y = x[idx1], y[idx2]
z = z[idx2,:]
z = z[:,idx1]
X,Y = np.meshgrid(x, y, sparse=False, indexing='xy')

# Now create NC file
with nc.Dataset(outfile, 'w') as fout:
    fout.createDimension('x', x.size)
    fout.createDimension('y', y.size)
    # Write each variable
    def writevar(name,val):
        fout.createVariable(name, 'f8', ('y','x'), zlib=True, complevel=4)
        fout.variables[name][:] = val
    writevar('nav_lon', X);
    writevar('nav_lat', Y)
    writevar('Bathymetry', z)
