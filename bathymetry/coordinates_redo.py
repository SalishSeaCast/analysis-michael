from IPython import embed
from coordinates_helpers import *
import scipy.io as sio
import netCDF4 as nc
import numpy as np
import datetime

# This script rebuilds the coordinates file for SalishSea based on the
# original seagrid output found in a .mat file in the NEMO_PREPARATION archive.
#
# We use the same grid points as in the original coordinates file, but
# recompute the scaling factors because they were not correct in the
# original file.
#
# For reference, the Arakawa-C grid looks like this:
#           f---v---f---v---f  
#           |       |       |
#           u   T   u   T   u
#           |       |       |
#           f---v---f---v---f
#

# Load JP's original seagrid output, which are our T points
# After transposing, the dimensions are 901x401 while SalishSea is 898x398.
# We drop the last three in both dimensions when writing which gives us the
# same grid as in the original coordinates file. However we keep the extra
# points during the computations so we can get scaling factors at the right/top.
mfile = sio.loadmat('seagrid_west_coast_1km_900x400_rot_new.mat')
glamt = np.transpose( mfile['s'][0,0]['geographic_grids'][0,0] )
gphit = np.transpose( mfile['s'][0,0]['geographic_grids'][0,1] )

# Compute the rest of the grid points
glamu,gphiu = t2u(glamt,gphit)
glamv,gphiv = t2v(glamt,gphit)
glamf,gphif = t2f(glamt,gphit)

# Compute scaling factors (with extrapolation for the left/bottom most scaling factor)
#
e1t = gete1(glamu,gphiu,expandleft=True)   # Need a left u point
e1u = gete1(glamt,gphit)
e1v = gete1(glamf,gphif,expandleft=True)   # Need a left f point
e1f = gete1(glamv,gphiv)
#
e2t = gete2(glamv,gphiv,expanddown=True)   # Need a lower v point
e2u = gete2(glamf,gphif,expanddown=True)   # Need a lower f point
e2v = gete2(glamt,gphit)
e2f = gete2(glamu,gphiu)

# Output slices
J,I = slice(0,898), slice(0,398)

#
filename = "coordinates_seagrid_SalishSea3.nc"
writecoords(filename,
            glamt[J,I],glamu[J,I],glamv[J,I],glamf[J,I],
            gphit[J,I],gphiu[J,I],gphiv[J,I],gphif[J,I],
            e1t[J,I],e1u[J,I],e1v[J,I],e1f[J,I],
            e2t[J,I],e2u[J,I],e2v[J,I],e2f[J,I])

# Add note to history
cnc = nc.Dataset(filename, 'r+')
note ='[{}] Rebuilt with correct scaling factors (e1* and e2*)'
cnc.setncattr('history', note.format(datetime.datetime.today().strftime('%Y-%m-%d')))
cnc.close()
