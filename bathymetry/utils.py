"""
Some utility functions
"""
import numpy as np

def slidingmean(t,d,w):
    """
    Sliding mean, computed as:  sm(t) = (1/w) * \int_{t-w}^{t} d(t') dt'
    - t and d should be length N vectors, w is the averaging window width
    """
    from slidingmean_c import slidingmean_c as smean
    tt, dd = t.copy(), d.copy()  # Interfacing with C routines requires care...
    ww = np.array([float(w)])
    return smean(tt,dd,ww)
