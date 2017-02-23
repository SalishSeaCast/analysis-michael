#River comparator
from IPython import embed
import netCDF4 as nc
import numpy as np
import cmocean.cm as cm
import matplotlib.pyplot as plt
import sys
import numpy.ma as ma

from salishsea_tools import river_201702 as wsnew
from salishsea_tools import rivertools as wsold

old_rivers_file = nc.Dataset('/home/mdunphy/MEOPAR/NEMO-forcing/rivers/rivers_month.nc') #_allArms.nc
old_runoff = old_rivers_file.variables['rorunoff'][0]
old_rivers_file.close()

bathy_file = nc.Dataset('/home/mdunphy/MEOPAR/NEMO-forcing/grid/bathymetry_201702.nc')
bathy = bathy_file.variables['Bathymetry'][:]
bathy_file.close()
mmfile = nc.Dataset('/home/mdunphy/MEOPAR/NEMO-forcing/grid/mesh_mask201702.nc')
tmask = mmfile.variables['tmask'][0,0,...]
mmfile.close()


oldbathy_file = nc.Dataset('/home/mdunphy/MEOPAR/NEMO-forcing/grid/bathy_downonegrid2.nc')
oldbathy = oldbathy_file.variables['Bathymetry'][:]
oldbathy_file.close()

oldmmfile = nc.Dataset('/home/mdunphy/MEOPAR/NEMO-forcing/grid/mesh_mask_downbyone2.nc')
oldtmask = oldmmfile.variables['tmask'][0,0,...]
olde12t = oldmmfile.variables['e1t'][0,...]*oldmmfile.variables['e2t'][0,...]
oldmmfile.close()

oldbathy=np.clip(oldbathy,0,80)
bathy=np.clip(bathy,0,80)

# Mask by tmask from the mesh mask file
bathy = ma.masked_array(bathy.data, mask=1-tmask)
oldbathy = ma.masked_array(oldbathy.data, mask=1-oldtmask)

#np.sum(old_runoff*olde12t)
#6337869.4791815793
#np.sum(old_runoff*olde12t*oldtmask)
#6112259.3660572302


def getChar(): # Found this getChar function online
    import tty, termios # raises ImportError if unsupported
    fd = sys.stdin.fileno()
    oldSettings = termios.tcgetattr(fd)
    try:
        tty.setraw(fd)
        answer = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, oldSettings)
    return answer
 
def getaction():
    print("Press arrow keys to move river, enter to continue, or q to quit: ")
    c = getChar()
    if c is 'q' or c is 'Q': return "quit"
    if c is 'e' or c is 'E': return "embed"
    if c == '\x1b': # Assume it's the first byte of a three-byte arrow key
        k=''.join([c,getChar(),getChar()])  # Merge it with the other two bytes
        if k=='\x1b[A': return "up"
        if k=='\x1b[B': return "down"
        if k=='\x1b[C': return "right"
        if k=='\x1b[D': return "left"
    return "next"

fig = plt.figure(figsize=(12, 6))
plt.clf(); plt.ion(); plt.show(); plt.pause(0.1)
cmap = cm.deep
cmap.set_bad('darkgreen')

def plotrvr(i,j,i0,j0):
    plt.clf(); axs = plt.subplot(1,2,1), plt.subplot(1,2,2)

    NY,NX = oldbathy.shape
    x=[q-0.5 for q in range(NX)]
    y=[q-0.5 for q in range(NY)]

    axs[0].pcolormesh(x,y,oldbathy, cmap=cmap)
    axs[0].set_title("Old bathy")
    axs[1].pcolormesh(x,y,bathy, cmap=cmap)
    axs[1].set_title("New bathy")

#    jj,ii = np.where(oldbathy>0)
#    axs[0].plot(ii,jj, 'b.')
#    jj,ii = np.where(oldbathy==0)
#    axs[0].plot(ii,jj, 'k.')

    jj,ii = np.where(old_runoff>0)
    axs[0].scatter(ii,jj, facecolors='none', edgecolors='r', label='old runoff cells')
    axs[1].plot(i,j,'r*',markeredgecolor='r',label='new locn')
    
    if i0 is not None:
        axs[0].plot(i0,j0,'k*',label='old locn')
        axs[1].plot(i0,j0,'k*',label='old locn')
        axs[1].plot([i0,i],[j0,j],'b',label='adjustment')


    d=10   # viewing window half width
    xlo,xhi = max(0 ,i-d), min(NX,i+d)
    ylo,yhi = max(0 ,j-d), min(NY,j+d)
    for i in range(2):
        axs[i].set_xlim((xlo,xhi))
        axs[i].set_ylim((ylo,yhi))
        axs[i].legend(loc='best')
    plt.pause(0.1)

# OK now we get to work
watersheds = ['skagit', 'fraser', 'evi_n', 'howe', 'bute', 'puget', 'jdf', 'evi_s', 'jervis', 'toba']
for watershed in watersheds:
    old = wsold.get_watershed_prop_dict(watershed)
    new = wsnew.prop_dict[watershed]
    for rname in new.keys():
        print("Watershed {}, river {}".format(watershed,rname))
        # Caution: i and j are swapped in the prop dicts
        i1,j1 = new[rname]['j'],new[rname]['i']
        try:
            i0,j0 = old[rname]['j'],old[rname]['i']
        except:
            print("Watershed {}, river {}: not found in old prop dict, continuing ...".format(watershed,rname))
            i0,j0=None,None
            continue
        print("Old position {}".format((i0,j0)))
        print("New position {}".format((i1,j1)))
        di,dj = 0,0
        while False:
            plotrvr(i1,j1,i0,j0)
            action = getaction()
            if action is "quit": exit()
            if action is "embed": embed();exit()
            if action is "next": break
            if action is "left":  di-=1
            if action is "right": di+=1
            if action is "down":  dj-=1
            if action is "up":    dj+=1
            print("Moving {}".format(action))
            # Caution: i and j are swapped in the prop dicts
            i1,j1 = new[rname]['j']+di, new[rname]['i']+dj
        
        # OK we're done with this river
        if i0==i1 and j0==j1:
            print("Finished with river {}, no adjustment".format(rname))
        else:
            msg="Watershed {}: moved river {:20s} from {} to {}".format(watershed,rname,(j0,i0),(j1,i1))
            if di!=0 or dj!=0:
                # We moved the river, so log the edits that we need to make in river_201702.py
                with open("river_updates_for_201702.txt",'a') as f:
                    f.write(msg+"\n")
                print(msg+" [LOGGED]")
            else:
                print(msg)
        print("\n")
