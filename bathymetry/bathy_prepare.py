import os
from bathy_helpers import *
from bathy_readers import *
from matplotlib import path
import numpy as np
from pyproj import Proj
from scipy import interpolate

"""
Routines to re-grid the bathy datasets on the SalishSea grid
"""
def prepare_cascadia(cascadiafile,glamt,gphit,glamf,gphif,glamfe,gphife):
    cache=cascadiafile+".results.npz"
    if os.path.exists(cache):
        return cache
    # Work with Cascadia data
    x,y,z,p = getcascadia(cascadiafile)

    # Convert NEMO coords to Cascadia projected coordinates
    NY, NX = glamt.shape[0], glamt.shape[1]
    Xt, Yt = p(glamt, gphit)
    Xf, Yf = p(glamf, gphif)
    Xfe, Yfe = p(glamfe, gphife)

    # There are too many points, so we filter based on minimum and maximum x,y
    idx1 = (x >= np.min(Xfe)) & (x <= np.max(Xfe))
    idx2 = (y >= np.min(Yfe)) & (y <= np.max(Yfe))
    x,y = x[idx1], y[idx2]
    z = z[idx2,:]
    z = z[:,idx1]
    X,Y = np.meshgrid(x, y, sparse=False, indexing='xy')
    X,Y,z = X.flatten(), Y.flatten(), z.flatten()

    # Further filter points not within the SalishSea domain
    poly0 = path.Path(makebox(Xfe,Yfe,0,NX,0,NY), closed=True)
    tmp = np.zeros([X.size,2]); tmp[:,0]=X; tmp[:,1]=Y;
    idx = poly0.contains_points(tmp)
    X, Y, z = X[idx], Y[idx], z[idx]

    # Construct new bathy using interpolation
    points = (X,Y)
    xi = (Xt.flatten(), Yt.flatten())
    casnearest = np.reshape(interpolate.griddata(points, z, xi, method='nearest'), Xt.shape)
    caslinear = np.reshape(interpolate.griddata(points, z, xi, method='linear'), Xt.shape)

    Xo,Yo = X-125,Y-125
    points = (Xo,Yo)
    casnearestjp = np.reshape(interpolate.griddata(points, z, xi, method='nearest'), Xt.shape)
    
    # Get the bin indices, apply bin methods
    boxi, boxj = getboxij(Xfe,Yfe,X,Y,cache=cascadiafile+".boxij.npz",searchmore=False)
    casmin,casmax,casmean,casmedian,cascount = binstobathy(boxi,boxj,X,Y,z,NX,NY)

    np.savez(cache, casmin=casmin, casmax=casmax, casmean=casmean, casmedian=casmedian,
                    cascount=cascount, casnearest=casnearest, caslinear=caslinear, casnearestjp=casnearestjp)
    return cache
#1239279 points in 407.3525502681732 s
#hit0 report: 466518/1239279=37.644307698266495%
#hit1 report: 768443/1239279=62.00726390102632%
#hit2 report: 0/1239279=0.0%
#hit3 report: 0/1239279=0.0%
#hit4 report: 0/1239279=0.0%
#no hit report: 4318/1239279=0.3484284007071854%
#points not in domain: 0

def prepare_chs(chsfile,glamt,glamf,gphif,gphit,glamfe,gphife):
    cache=chsfile+".results.npz"
    if os.path.exists(cache):
        return cache
    # Work with CHS data
    x,y,z,p=getchs(chsfile)

    # Convert NEMO coords to UTM Zone 10 coordinates
    NY, NX = glamt.shape[0], glamt.shape[1]
    Xt, Yt = p(glamt, gphit)
    Xf, Yf = p(glamf, gphif)
    Xfe, Yfe = p(glamfe, gphife)

    # Construct new bathy using interpolation
    points = (x,y)
    xi = (Xt.flatten(), Yt.flatten())
    chsnearest = np.reshape(interpolate.griddata(points, z, xi, method='nearest'), Xt.shape)
    chslinear = np.reshape(interpolate.griddata(points, z, xi, method='linear'), Xt.shape)

    # Get the bin indices, apply bin methods
    boxi, boxj = getboxij(Xfe,Yfe,x,y,cache=chsfile+".boxij.npz",searchmore=True)
    chsmin,chsmax,chsmean,chsmedian,chscount = binstobathy(boxi,boxj,x,y,z,NX,NY)
    np.savez(cache, chsmin=chsmin, chsmax=chsmax, chsmean=chsmean, chsmedian=chsmedian,
                    chscount=chscount, chsnearest=chsnearest, chslinear=chslinear)
    return cache
#13108913 points in 908.1155421733856 s
#hit0 report: 12013979/13108913=91.64740814131576%
#hit1 report: 969689/13108913=7.3971732057417725%
#hit2 report: 12159/13108913=0.09275368598449009%
#hit3 report: 39251/13108913=0.29942223279687646%
#hit4 report: 20135/13108913=0.15359778495745605%
#no hit report: 236/13108913=0.0018003018251780297%
#points not in domain: 53464                    
                    
             
             
def prepare_chs2(chs2path,glamt,gphit,glamf,gphif,glamfe,gphife):
    cache=chs2path+"results.npz"
    if os.path.exists(cache):
        return cache
    # Work with CHS2 data
    chs2files=['Area-A_MB-10m_CGVD28.txt','Area-A_SB_CGVD28.txt','Area-B_MB-10m_CGVD28.txt']
    chs2files+=['Area-B_SB_CGVD28.txt','Fraser_NAD83_CGVD28_all.txt','S-of-G_Data_10M_CGVD28.txt']
    x,y,z=np.zeros([0]),np.zeros([0]),np.zeros([0])
    for fn in chs2files:
        xn,yn,zn,p=getchs2(chs2path + fn)
        x=np.append(x,xn)
        y=np.append(y,yn)
        z=np.append(z,zn)

    NY, NX = glamt.shape[0], glamt.shape[1]
    Xt, Yt = glamt, gphit
    Xf, Yf = glamf, gphif
    Xfe, Yfe = glamfe, gphife

    # Construct new bathy using interpolation
    points = (x,y)
    xi = (Xt.flatten(), Yt.flatten())
    chs2nearest = np.reshape(interpolate.griddata(points, z, xi, method='nearest'), Xt.shape)
    chs2linear = np.reshape(interpolate.griddata(points, z, xi, method='linear'), Xt.shape)

    # Get the bin indices, apply bin methods
    boxi, boxj = getboxij(Xfe,Yfe,x,y,cache=chs2path+"boxij.npz",searchmore=True)
    chs2min,chs2max,chs2mean,chs2median,chs2count = binstobathy(boxi,boxj,x,y,z,NX,NY)
    np.savez(cache, chs2min=chs2min, chs2max=chs2max, chs2mean=chs2mean, chs2median=chs2median,
             chs2count=chs2count, chs2nearest=chs2nearest, chs2linear=chs2linear)
    return cache
#26676210 points in 1197.0264418125153 s
#hit0 report: 25376536/26676210=95.12796607913943%
#hit1 report: 981883/26676210=3.6807440037396617%
#hit2 report: 49917/26676210=0.18712178379162558%
#hit3 report: 91056/26676210=0.341337843719179%
#hit4 report: 2525/26676210=0.009465362583365479%
#no hit report: 2525/26676210=0.009465362583365479%
#points not in domain: 171768
             
def prepare_bc3(bc3file,glamt,gphit,glamf,gphif,glamfe,gphife):
    cache=bc3file+".results.npz"
    if os.path.exists(cache):
        return cache
    # Work with BC3 data
    x,y,z,p = getbc3(bc3file)

    NY, NX = glamt.shape[0], glamt.shape[1]

    # There are too many points, so we filter based on minimum and maximum x,y
    idx1 = (x >= np.min(glamfe)) & (x <= np.max(glamfe))
    idx2 = (y >= np.min(gphife)) & (y <= np.max(gphife))
    x,y = x[idx1], y[idx2]
    z = z[idx2,:]
    z = z[:,idx1]
    X,Y = np.meshgrid(x, y, sparse=False, indexing='xy')
    X,Y,z = X.flatten(), Y.flatten(), z.flatten()

    # Further filter points not within the SalishSea domain
    poly0 = path.Path(makebox(glamfe,gphife,0,NX,0,NY), closed=True)
    tmp = np.zeros([X.size,2]); tmp[:,0]=X; tmp[:,1]=Y;
    idx = poly0.contains_points(tmp)
    X, Y, z = X[idx], Y[idx], z[idx]

    # Construct new bathy using interpolation
    points = (X,Y)
    xi = (glamt.flatten(), gphit.flatten())
    bc3nearest = np.reshape(interpolate.griddata(points, z, xi, method='nearest'), glamt.shape)
    bc3linear = np.reshape(interpolate.griddata(points, z, xi, method='linear'), glamt.shape)

    # Get the bin indices, apply bin methods
    boxi, boxj = getboxij(glamfe,gphife,X,Y,cache=bc3file+".boxij.npz",searchmore=False)   # takes 15-20 mins
    bc3min,bc3max,bc3mean,bc3median,bc3count = binstobathy(boxi,boxj,X,Y,z,NX,NY)

    np.savez(cache, bc3min=bc3min, bc3max=bc3max, bc3mean=bc3mean, bc3median=bc3median,
                    bc3count=bc3count, bc3nearest=bc3nearest, bc3linear=bc3linear)
    return cache
#10309989 points in 1063.4851875305176 s
#hit0 report: 8509929/10309989=82.54062152733626%
#hit1 report: 1790466/10309989=17.366323087250628%
#hit2 report: 0/10309989=0.0%
#hit3 report: 0/10309989=0.0%
#hit4 report: 0/10309989=0.0%
#no hit report: 9594/10309989=0.09305538541311732%
#points not in domain: 0
    
    
    
    
## Feb 14 rerun with coordinates4
#cascadiafile = '/home/mdunphy/MEOPAR/Bathy/Cascadia/cascadia.bil'
#cache=prepare_cascadia(cascadiafile,glamt,gphit,glamf,gphif,glamfe,gphife)
#locals().update(np.load(cache))
#
#1239279 points in 406.56440019607544 s
#hit0 report: 466528/1239279=37.64511461906479%
#hit1 report: 768433/1239279=62.00645698022802%
#hit2 report: 0/1239279=0.0%
#hit3 report: 0/1239279=0.0%
#hit4 report: 0/1239279=0.0%
#no hit report: 4318/1239279=0.3484284007071854%
#points not in domain: 0
#
#chsfile = '/home/mdunphy/MEOPAR/Bathy/CHS/Salish Sea 25m Grid.txt'
#cache=prepare_chs(chsfile,glamt,gphit,glamf,gphif,glamfe,gphife)
#locals().update(np.load(cache))
#
#13108913 points in 834.6945581436157 s
#hit0 report: 12013971/13108913=91.64734711413524%
#hit1 report: 969681/13108913=7.397112178561258%
#hit2 report: 12161/13108913=0.09276894277961872%
#hit3 report: 39265/13108913=0.29952903036277684%
#hit4 report: 20135/13108913=0.15359778495745605%
#no hit report: 236/13108913=0.0018003018251780297%
#points not in domain: 53464
#
#chs2path='/home/mdunphy/MEOPAR/Bathy/Data_From_Mitchell_CHS/'
#cache=prepare_chs2(chs2path,glamt,gphit,glamf,gphif,glamfe,gphife)
#locals().update(np.load(cache))
#
#26676210 points in 1270.055364370346 s
#hit0 report: 25370822/26676210=95.10654624476265%
#hit1 report: 988041/26676210=3.703828242467727%
#hit2 report: 49749/26676210=0.1864920091722175%
#hit3 report: 90778/26676210=0.3402957166703966%
#hit4 report: 2527/26676210=0.009472859900263194%
#no hit report: 2525/26676210=0.009465362583365479%
#points not in domain: 171768
#
#bc3file = '/home/mdunphy/MEOPAR/Bathy/BC3/british_columbia_3sec.asc'
#cache=prepare_bc3(bc3file,glamt,gphit,glamf,gphif,glamfe,gphife)
#locals().update(np.load(cache))
#
#10309989 points in 1094.212818622589 s
#hit0 report: 8509821/10309989=82.5395739995455%
#hit1 report: 1790574/10309989=17.367370615041395%
#hit2 report: 0/10309989=0.0%
#hit3 report: 0/10309989=0.0%
#hit4 report: 0/10309989=0.0%
#no hit report: 9594/10309989=0.09305538541311732%
#points not in domain: 0
