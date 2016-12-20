from matplotlib import path
import numpy as np
import time,os
from IPython import embed

def expandf(glamf, gphif):
    # Expand the f points grid so the f points form complete boxes around the t points
    # This is needed because the coordinates file truncates the f points by one.
    NY, NX = glamf.shape[0], glamf.shape[1]
    glamfe = np.zeros([NY+1, NX+1])
    gphife = np.zeros([NY+1, NX+1])
    # Long
    glamfe[1:,1:] = glamf
    glamfe[0,1:] = glamf[0,:] - (glamf[1,:] - glamf[0,:])     # extraoplation
    glamfe[:,0] = glamfe[:,1] - (glamfe[:,2] - glamfe[:,1])   # extraoplation
    # Lat
    gphife[1:,1:] = gphif
    gphife[0,1:] = gphif[0,:] - (gphif[1,:] - gphif[0,:])     # extraoplation
    gphife[:,0] = gphife[:,1] - (gphife[:,2] - gphife[:,1])   # extraoplation
    return glamfe, gphife

# Helper functions
def makebox(glamfe,gphife,il,ir,jl,jr):
    # Builds a bounding polygon enclosing grid boxes il..ir by jl..jr
    n, L = 0, 2*(ir-il+1) + 2*(jr-jl-1) + 1
    px, py = np.zeros([L,1]), np.zeros([L,1])
    for i in range(il,ir+1):
        px[n] = glamfe[jl,i]
        py[n] = gphife[jl,i]
        n+=1
    for j in range(jl+1,jr+1):
        px[n] = glamfe[j,ir]
        py[n] = gphife[j,ir]
        n+=1
    for i in range(ir-1,il-1,-1):
        px[n] = glamfe[jr,i]
        py[n] = gphife[jr,i]
        n+=1
    for j in range(jr-1,jl-1,-1):
        px[n] = glamfe[j,il]
        py[n] = gphife[j,il]
        n+=1
    if n < L:
        print("n,L:",n,l)
    return np.hstack((px,py))

def search(glamfe,gphife,cx,cy,il,ir,jl,jr):
    # Conducts a binary search to find which grid box contains point cx,cy
    while (ir-il > 1) or (jr-jl > 1):
        # Do a refinement in i
        if ir-il>1:
            irr = il + (ir-il)//2
            p = makebox(glamfe,gphife,il,irr,jl,jr)
            poly = path.Path(p, closed=True)
            test = poly.contains_points([(cx,cy)])
            if test: ir = irr
            else: il = irr
        # Do a refinement in j
        if jr-jl > 1:
            jrr = jl + (jr-jl)//2
            p = makebox(glamfe,gphife,il,ir,jl,jrr)
            poly = path.Path(p, closed=True)
            test = poly.contains_points([(cx,cy)])
            if test: jr = jrr
            else: jl = jrr
    p = makebox(glamfe,gphife,il,ir,jl,jr)
    poly = path.Path(p, closed=True)
    return il,jl,poly

def expand(glamfe,gphife,cx,cy,il,jl,e):
    # expands the search box, see comments below in getboxij
    NY, NX = glamfe.shape[0]-1, glamfe.shape[1]-1
    i1, i2 = max(0,il-e), min(NX-1,il+e)
    j1, j2 = max(0,jl-e), min(NY-1,jl+e)
    p1 = makebox(glamfe,gphife,i1,i2,j1,j2)
    poly1 = path.Path(p1, closed=True)
    test = poly1.contains_points([(cx,cy)])
    return i1,i2,j1,j2,test

def getboxij(glamfe,gphife,lon,lat,nmax=-1,searchmore=False,cache=None):
    if cache is not None and os.path.exists(cache):
        npzfile = np.load(cache)
        boxi, boxj = npzfile['boxi'], npzfile['boxj']
    else:
        # The goal here is to find the model grid box that each input point falls inside
        # We record this information as the i,j coordinate in arrays boxi,boxj
        t0 = time.time()
        hit0, hit1, hit2, hit3, hit4, nout = 0, 0, 0, 0, 0, 0
        NY, NX = glamfe.shape[0], glamfe.shape[1]

        p0 = makebox(glamfe,gphife,0,NX-1,0,NY-1)
        poly0 = path.Path(p0, closed=True)

        boxi, boxj = -1*np.ones(lon.shape), -1*np.ones(lon.shape)
        if nmax>0: N = nmax
        else: N=len(lon)

        for i in range(N):
            cx,cy = lon[i],lat[i]

            if i>0 and 'poly' in locals(): # If we're on the second point or later ...

                # ... then it's quite likely that the current point is in the same box as the previous point.
                # So we check if it's in the same box first, because it's faster that running the complete binary search
                test = poly.contains_points([(cx,cy)])
                if test:
                    hit0 += 1
                    boxi[i], boxj[i] = il, jl
                    continue

                # ... and if we're not in the same box, it's still quite likely that we're in an adjacent box
                # and it's still faster to check nearby boxes than to do the complete search
                i1,i2,j1,j2,test = expand(glamfe,gphife,cx,cy,il,jl,2)
                if test:
                    hit1+=1
                    il,jl,poly = search(glamfe,gphife,cx,cy,i1,i2,j1,j2)
                    boxi[i], boxj[i] = il, jl
                    continue

                if searchmore:
                    # as before with an expanded search box
                    i1,i2,j1,j2,test = expand(glamfe,gphife,cx,cy,il,jl,4)
                    if test:
                        hit2+=1
                        il,jl,poly = search(glamfe,gphife,cx,cy,i1,i2,j1,j2)
                        boxi[i], boxj[i] = il, jl
                        continue

                    # as before with an expanded search box
                    i1,i2,j1,j2,test = expand(glamfe,gphife,cx,cy,il,jl,20)
                    if test:
                        hit3+=1
                        il,jl,poly = search(glamfe,gphife,cx,cy,i1,i2,j1,j2)
                        boxi[i], boxj[i] = il, jl
                        continue

                    # as before with an expanded search box
                    i1,i2,j1,j2,test = expand(glamfe,gphife,cx,cy,il,jl,50)
                    if test:
                        hit4+=1
                        il,jl,poly = search(glamfe,gphife,cx,cy,i1,i2,j1,j2)
                        boxi[i], boxj[i] = il, jl
                        continue


            # If we get here, then this is either the first point, or we're not very close to the previous point,
            # so we start the search from scratch
            test = poly0.contains_points([(cx,cy)])  # Ensure the point is inside the domain
            if test:
                # We're in the domain, so start the binary search from the complete domain
                il,jl,poly = search(glamfe,gphife,cx,cy,0,NX-1,0,NY-1)
                boxi[i], boxj[i] = il, jl
            else:
                nout += 1  # Found a point outside of the domain, increment counter but otherwise do nothing
                continue

        # Summary of how often each search was used (hit0 is fastest, followed by hit1, ...)
        i += 1
        t = time.time()-t0
        nohit = i-hit0-hit1-hit2-hit3-hit4-nout
        print("{} points in {} s".format(i,t))
        if nmax > 0:
            print("ETA: {:.3} min".format( (len(lon)/nmax)*t/60 ))
        print("hit0 report: {}/{}={}%".format(hit0,i,100*hit0/i))
        print("hit1 report: {}/{}={}%".format(hit1,i,100*hit1/i))
        print("hit2 report: {}/{}={}%".format(hit2,i,100*hit2/i))
        print("hit3 report: {}/{}={}%".format(hit3,i,100*hit3/i))
        print("hit4 report: {}/{}={}%".format(hit4,i,100*hit4/i))
        print("no hit report: {}/{}={}%".format(nohit,i,100*nohit/i))
        print("points not in domain: {}".format(nout))
        np.savez(cache, boxi=boxi, boxj=boxj)

    return boxi.astype(np.int32), boxj.astype(np.int32)

def binstobathy(boxi,boxj,x,y,z,NX,NY):
    # Now we use the boxi,boxj information to assemble the binned input bathy on NEMO grid
    bmin, bmax, bmean = np.zeros([NY,NX]), np.zeros([NY,NX]), np.zeros([NY,NX])
    bmedian, bcount = np.zeros([NY,NX]), np.zeros([NY,NX])
    xx, yy, zz = x.flatten(), y.flatten(), z.flatten()
    for j in range(NY):
        idx = (boxj == j)
        if np.any(idx):   # j coordinate matches
            itmp = boxi[idx]
            xtmp = xx[idx]
            ytmp = yy[idx]
            ztmp = zz[idx]
            for i in range(NX): # i coordinate matches
                idx = (itmp == i)
                if np.any(idx):
                    # We've arrived at grid box i,j and we have at least 1 data point
                    d            = ztmp[idx]
                    bmin[j,i]    = np.min(d)
                    bmax[j,i]    = np.max(d)
                    bmean[j,i]   = np.mean(d)
                    bmedian[j,i] = np.median(d)
                    bcount[j,i]  = len(d)
    return bmin,bmax,bmean,bmedian,bcount.astype(np.int32)
