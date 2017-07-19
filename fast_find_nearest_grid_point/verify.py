import netCDF4 as nc
from fast_ll2ij_SalishSea201702 import fast_ll2ij_SalishSea201702
from salishsea_tools.geo_tools import find_closest_model_point

with nc.Dataset('/home/mdunphy/MEOPAR/NEMO-forcing/grid/coordinates_seagrid_SalishSea201702.nc') as nbl:
    glamt = nbl.variables['glamt'][0,...]
    gphit = nbl.variables['gphit'][0,...]

#@profile
def dothetest():
    for j in range(glamt.shape[0]):
        for i in range(glamt.shape[1]):
            jfast, ifast = fast_ll2ij_SalishSea201702(glamt[j,i], gphit[j,i], glamt, gphit)
            jgeo, igeo = find_closest_model_point(glamt[j,i], gphit[j,i], glamt, gphit)
            if jfast != jgeo or ifast != igeo:
                print("Mismatch:", ifast,jfast,igeo,jgeo)
dothetest()

# kernprof reveals that the inversion approach is ~7x faster than the geo_tools function

#Total time: 702.584 s
#File: verify.py
#Function: dothetest at line 9
#
#Line #      Hits         Time  Per Hit   % Time  Line Contents
#==============================================================
#     9                                           @profile
#    10                                           def dothetest():
#    11       899          443      0.5      0.0      for j in range(glamt.shape[0]):
#    12    358302       189278      0.5      0.0          for i in range(glamt.shape[1]):
#    13    357404     86849607    243.0     12.4              jfast, ifast = fast_ll2ij_SalishSea201702(glamt[j,i], gphit[j,i], glamt, gphit)
#    14    357404    615135587   1721.1     87.6              jgeo, igeo = find_closest_model_point(glamt[j,i], gphit[j,i], glamt, gphit)
#    15    357404       408594      1.1      0.1              if jfast != jgeo or ifast != igeo:
#    16                                                           print("Mismatch:", ifast,jfast,igeo,jgeo)
