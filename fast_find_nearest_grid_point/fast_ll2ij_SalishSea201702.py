import numpy as np
from salishsea_tools.geo_tools import haversine

# Magic numbers from GridPlaneFit order 2
Ci = np.array([ 1.26921890e+04,   1.59279223e+02,   1.39650621e+02,
                3.21769836e-01,   1.37597009e+00,   1.58488524e+00])
Cj = np.array([-1.26067300e+04,   6.53817493e+00,   1.76514195e+02,
                4.90478919e-01,   9.06976896e-01,   1.33030515e+00])

def fast_ll2ij_SalishSea201702(lon,lat,model_lons,model_lats):
    def useC(X,Y,C):
        return C[0] + C[1]*X + C[2]*Y + C[3]*X**2 + C[4]*X*Y + C[5]*Y**2
    # Use coefficients to get an estimate of i and j
    ie = int(useC(lon,lat,Ci))
    je = int(useC(lon,lat,Cj))

    # Define search box, +/- 4 and + /- 19 is based on max error during fit
    i1 = max(0, ie - 4)
    i2 = min(ie + 4, model_lons.shape[1])
    j1 = max(0, je - 19)
    j2 = min(je + 19, model_lons.shape[0])

    # Prepare lists of test points in search box and their indices
    lons = model_lons[j1:j2,i1:i2].flatten()
    lats = model_lats[j1:j2,i1:i2].flatten()
    i_list = np.array([x for x in range(i1,i2)]*(j2-j1))
    j_list = j1 + np.array([x for x in range(0,(j2-j1)*(i2-i1))]) // (i2-i1)

    # Test all of these points, return closest point
    dists = haversine(
        np.array([lon] * i_list.size), np.array([lat] * j_list.size),
        lons, lats)
    n = dists.argmin()
    j, i = map(np.asscalar, (j_list[n], i_list[n]))
    return j, i
