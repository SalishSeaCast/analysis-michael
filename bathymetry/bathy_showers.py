import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import cmocean
from bathy_helpers import *

# Plotting helper funtions

def dezero(x):
    return x.flatten()[x.flatten()>0]

def mm(x):
    return ma.masked_array(x, mask=(x == 0 | np.isnan(x)))

def cmp(a,b,c,d,e):  # compare a and b, mask from c, colorbar limits d, title e
    im=plt.pcolormesh(c*0 + a -b); cb=plt.colorbar(im);
    im.set_cmap('seismic'); im.set_clim([-d,d]); cb.set_clim(-d,d); plt.title(e)
def shw(a,b,c,d):  # plot a, mask b, colorbar limits c, title d
    if b is None: b=mm(a)
    im=plt.pcolormesh(b*0+a); cb=plt.colorbar(im);
    im.set_cmap(cmocean.cm.haline); im.set_clim([0,c]); cb.set_clim(0,c); plt.title(d)


def sl(x):
    imin = 350; imax = 505; jmin = 280; jmax = 398
    islice = slice(imin,imax); jslice = slice(jmin, jmax)
    if x is not None: return x[islice,jslice]

# These two are same as cmp/shw, except with a Frasher river mouth zoom
def cmpz(a,b,c,d,e):  # compare a and b, mask from c, colorbar limits d, title e
    cmp(sl(a), sl(b), sl(c), d, e)
def shwz(a,b,c,d):  # plot a, mask b, colorbar limits c, title d
    shw(sl(a), sl(b), c, d)


# Also with Fraser river zoom, but with lat,lon, such that we can easily draw rivers
def shwll(xf,yf,a,b,c,d,z=True):  # plot a, mask b, colorbar limits c, title d, zoom
    if b is None: b=mm(a)
    xfe, yfe = expandf(xf, yf)
    im=plt.pcolormesh(xfe[:-1,:-1],yfe[:-1,:-1],b*0+a); cb=plt.colorbar(im);
    im.set_cmap(cmocean.cm.haline); im.set_clim([0,c]); cb.set_clim(0,c); plt.title(d)
    if z:
        plt.gca().set_xlim(-123.4,-122.6)
        plt.gca().set_ylim(48.9,49.5)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
