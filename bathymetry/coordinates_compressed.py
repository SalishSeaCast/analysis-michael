import scipy.io as sio
import netCDF4 as nc
import numpy as np
import datetime
from IPython import embed
from bathy_helpers import *
from coordinates_helpers import *

# SalishSea data
cnc = nc.Dataset('/home/mdunphy/MEOPAR/NEMO-forcing/grid/coordinates_seagrid_SalishSea.nc', 'r')
glamf = cnc.variables['glamf'][0,...]; gphif = cnc.variables['gphif'][0,...]
glamt = cnc.variables['glamt'][0,...]; gphit = cnc.variables['gphit'][0,...]
cnc.close()

# Box coordinates to extract from SalishSea grid for compression
i0,i1 = 324,388
j0,j1 = 402,446
ni,nj = i1-i0+1, j1-j0+1
print ("Dims for seagrid, first: {}, second: {}".format(nj,ni))


# Get f-point corners
crt=[(i0,j1),(i0,j0),(i1,j0),(i1,j1)]
crf=[(i0-1,j1),(i0-1,j0-1),(i1,j0-1),(i1,j1)]
xc=[glamf[j,i] for i,j in crf]
yc=[gphif[j,i] for i,j in crf]

# Write corners (counter clockwise) for SeaGrid
with open('/home/mdunphy/seagrid/FraserCorners.txt','w') as f:
    for xk,yk in zip(xc,yc):
        f.write('{:.10g} {:.10g} 1\n'.format(xk,yk))

# Now we go to seagrid in MATLAB
# 1) Load the coastline file PNWrivers-ll.mat
# 2) Load the boundary file FraserCorners.txt
# 3) Toggle MenuBar so you can zoom the Fraser region
# 3) View->Setup to set number of grid points (45 first, 65 second, as printed above)
# 4) Adjust north and south bounding line
# 5) Save as FraserCompressed.mat file

# Import the seagrid result (J,I indexing, flip in j dimension needed)
mfile = sio.loadmat('/home/mdunphy/seagrid/FraserCompressed.mat')
nlonf = np.flipud( mfile['s'][0,0]['geographic_grids'][0,0] )
nlatf = np.flipud( mfile['s'][0,0]['geographic_grids'][0,1] )


# SeaGrid compresses the grid cells too much along the curved boundaries
# So here we equi-space them in j direction while aligning with the original i grid
# This will reduce the orthogonality of the grid in the compressed region, however
# we deem this unimportant as it enables is to resolve the Fraser river arms
mlonf = np.copy(nlonf)
mlatf = np.copy(nlatf)
for i in range(mlonf.shape[1]):
    ys,yn=nlatf[0,i],nlatf[-1,i]  # y endpoints on curved boundaries
    xs=np.interp(ys, gphif[:,i+i0-1], glamf[:,i+i0-1]) # interp in i to remain
    xn=np.interp(yn, gphif[:,i+i0-1], glamf[:,i+i0-1]) # aligned with old grid
    mlonf[:,i] = np.linspace(xs,xn,nj+1)
    mlatf[:,i] = np.linspace(ys,yn,nj+1)
nlonf=mlonf
nlatf=mlatf    
# Move to t points
nlont = 0.25*(nlonf[0:-1,0:-1] + nlonf[1:,0:-1] + nlonf[0:-1,1:] + nlonf[1:,1:])
nlatt = 0.25*(nlatf[0:-1,0:-1] + nlatf[1:,0:-1] + nlatf[0:-1,1:] + nlatf[1:,1:])


# Now we follow the coordinates synthesis of coordinates_redo.py, except that we
# substitute our new compressed Fraser river t-points, and compute the rest of
# the grid points and scaling factors in the exact same fashion
mfile = sio.loadmat('seagrid_west_coast_1km_900x400_rot_new.mat')
glamt = np.transpose( mfile['s'][0,0]['geographic_grids'][0,0] )
gphit = np.transpose( mfile['s'][0,0]['geographic_grids'][0,1] )

# Insert compressed region t points
glamt[j0:j1+1,i0:i1+1] = nlont
gphit[j0:j1+1,i0:i1+1] = nlatt

# Compute the rest of the grid points
glamu,gphiu = t2u(glamt,gphit)
glamv,gphiv = t2v(glamt,gphit)
glamf,gphif = t2f(glamt,gphit)

# Compute scaling factors (with extrapolation for the left/bottom most scaling factor)
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
filename = "coordinates_seagrid_SalishSea4.nc"
writecoords(filename,
            glamt[J,I],glamu[J,I],glamv[J,I],glamf[J,I],
            gphit[J,I],gphiu[J,I],gphiv[J,I],gphif[J,I],
            e1t[J,I],e1u[J,I],e1v[J,I],e1f[J,I],
            e2t[J,I],e2u[J,I],e2v[J,I],e2f[J,I])

# Add note to history
cnc = nc.Dataset(filename, 'r+')
note ='[{}] Compressed Fraser river region'
cnc.setncattr('history', note.format(datetime.datetime.today().strftime('%Y-%m-%d')))
cnc.close()
