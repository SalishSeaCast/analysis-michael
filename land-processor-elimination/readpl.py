"""
    Script to process the results of mpp_optimiz_zoom_nc
"""
from __future__ import division    # For python2 compatibility
from IPython import embed
import numpy as np

# Get the results from processor.layout
jpni,jpnj,jpi,jpj,nw,nl=[],[],[],[],[],[]
with open('processor.layout','r') as f:
    for line in f:
        if "jpni=" in line:
            jpni.append(int(line.split()[1]))
            jpnj.append(int(line.split()[3]))
        elif "jpi=" in line:
            jpi.append(int(line.split()[1]))
            jpj.append(int(line.split()[3]))
        elif "processeurs mer" in line:
            nw.append(int(line.split()[4]))
        elif "processeurs terre" in line:
            nl.append(int(line.split()[4]))
        elif "choix optimum" in line:
            break

# Convert lists to np arrays
jpni = np.array(jpni)
jpnj = np.array(jpnj)
jpi = np.array(jpi)
jpj = np.array(jpj)
nw = np.array(nw)
nl = np.array(nl)

# Sanity check
if (jpni*jpnj == nw+nl).all() is False:
    print("Something strange happened ... !")

n  = nw+nl   # number of processors
r  = nw/n    # ratio of water processors to total
ar = np.fmax(jpi/jpj, jpj/jpi) # tile aspect ratio

def filt(idx):
    return jpni[idx],jpnj[idx],jpi[idx],jpj[idx],nw[idx],nl[idx],n[idx],r[idx],ar[idx]
    
# Filter out useless configurations
jpni,jpnj,jpi,jpj,nw,nl,n,r,ar = filt(jpni > 2)
jpni,jpnj,jpi,jpj,nw,nl,n,r,ar = filt(jpnj > 2)
jpni,jpnj,jpi,jpj,nw,nl,n,r,ar = filt(nw <= 384)
jpni,jpnj,jpi,jpj,nw,nl,n,r,ar = filt(r < 1)

# Produce the lookup table for salishsea command
LUT = np.vstack((jpni,jpnj,nw)).T
np.savetxt('salish.csv', LUT, fmt='%g', delimiter=',')

head = ("MPI breakdown  Water  Land  r      Tile size        Tile aspect\n")
line = ("=============  =====  ====  =====  ===============  ===========\n")

# Produce the complete table
with open('LPE-SalishSea-complete.rst', 'w') as f:
    title="All decompositions"
    stars='*'*len(title)
    f.write(stars+'\n'+title+'\n'+stars+'\n\n')
    f.write(line+head+line)
    jpni,jpnj,jpi,jpj,nw,nl,n,r,ar = filt(np.argsort(nw))  # Sort by nw
    for i, _ in enumerate(n):
        f.write("{:>3d}x{:<3d} = {:3d}   {:3d}   {:3d}   {:.3f}  {:>3d}x{:<3d} = {:5d}   {:.3f}\n"
               .format(jpni[i],jpnj[i],n[i],nw[i],nl[i],r[i],jpi[i],jpj[i],jpi[i]*jpj[i],ar[i]))
    f.write(line)

def writebest(f,k,jpni,jpnj,jpi,jpj,nw,nl,n,r,ar):
    def filt(idx):
        return jpni[idx],jpnj[idx],jpi[idx],jpj[idx],nw[idx],nl[idx],n[idx],r[idx],ar[idx]
    # Extract the configurations with k water processors
    jpni,jpnj,jpi,jpj,nw,nl,n,r,ar = filt(nw == k)
    if len(n)==0: return
    # Filter out cases with poor aspect ratio
    jpni,jpnj,jpi,jpj,nw,nl,n,r,ar = filt(ar < 1.15)
    if len(n)==0: return
    # Sort remaining by best land elimination
    jpni,jpnj,jpi,jpj,nw,nl,n,r,ar = filt(np.argsort(r))  # Sort by r
    # Select all entries with same elimination ratio
    jpni,jpnj,jpi,jpj,nw,nl,n,r,ar = filt(r == r[0])
    for i in range(len(n)):
        f.write("{:>3d}x{:<3d} = {:3d}   {:3d}   {:3d}   {:.3f}  {:>3d}x{:<3d} = {:5d}   {:.3f}\n"
               .format(jpni[i],jpnj[i],n[i],nw[i],nl[i],r[i],jpi[i],jpj[i],jpi[i]*jpj[i],ar[i]))

# Produce the table of preferred decompositions
with open('LPE-SalishSea-preferred.rst', 'w') as f:
    title="Preferred decompositions"
    stars='*'*len(title)
    f.write(stars+'\n'+title+'\n'+stars+'\n\n')
    f.write(line+head+line)
    for k in range(min(nw), 256+1):
        writebest(f,k,jpni,jpnj,jpi,jpj,nw,nl,n,r,ar)
    f.write(line)
