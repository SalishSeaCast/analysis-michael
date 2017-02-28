import scipy.io as sio
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import cmocean
import matplotlib.pyplot as plt
from IPython import embed
from bathy_helpers import *
from coordinates_helpers import *

pre='/home/mdunphy/MEOPAR/NEMO-forcing/grid/'
cnc = nc.Dataset(pre+'coordinates_seagrid_SalishSea3.nc', 'r')
glamf = cnc.variables['glamf'][0,...]; gphif = cnc.variables['gphif'][0,...]
glamt = cnc.variables['glamt'][0,...]; gphit = cnc.variables['gphit'][0,...]
cnc.close()
glamfe, gphife = expandf(glamf, gphif)

# Load M0 bathy that has erroneous land cells in Puget Sound region
with nc.Dataset('/home/mdunphy/MEOPAR/Bathy/bathy_meter_SalishSeaM0.nc', 'r') as bnc:
    bathy = bnc.variables['Bathymetry'][:]
    bathyM0 = np.copy(bathy)

# Rivers
mfile = sio.loadmat('/ocean/rich/more/mmapbase/bcgeo/PNWrivers.mat')
ncstr = mfile['ncst']
# Coastlines
mfile = sio.loadmat('/ocean/rich/more/mmapbase/bcgeo/PNW.mat')
ncst = mfile['ncst']

#Plotting functions
def gridlines(x,y,c):
    for j in range(x.shape[0]):
        plt.plot(x[j,:],y[j,:],c,linewidth=0.25)
    for i in range(x.shape[1]):
        plt.plot(x[:,i],y[:,i],c,linewidth=0.25)

def drawrivers():
    pb=np.zeros([5,2])
    pb[:,0] = [x0,x0,x1,x1,x0]
    pb[:,1] = [y0,y1,y1,y0,y0]
    poly = path.Path(pb, closed=True)
    idx=np.where(np.isnan(ncst[:,0]))[0]
    for i,val in enumerate(idx[:-1]):
        x = ncstr[idx[i]+1:idx[i+1],0]
        y = ncstr[idx[i]+1:idx[i+1],1]
        pts = np.zeros([x.shape[0],2])
        pts[:,0] = x
        pts[:,1] = y
        if poly.contains_points(pts).any():
            plt.plot(x,y,'r')
    plt.plot(ncst[:,0],ncst[:,1],'r')

# Locations
coords, cid = [], []
def onclick(event):
    global coords
    global fig
    if event.inaxes is None:
        fig.canvas.mpl_disconnect(cid)
    else:
        coords.append((event.xdata, event.ydata))
        shw()
    return coords

def shw():
    # Ghastly function to print coords as python code
    tups=[]
    for c in coords:
        cx,cy = c
        i,j,xxxxx = search(glamfe,gphife,cx,cy,0,NX-1,0,NY-1)
        plt.plot(glamt[j,i], gphit[j,i], 'k*')
        tups+=['({:d},{:d}),'.format(i,j)]
    print('\n')
    s=''.join(tups)
    print('Fill = ['+s+']')

 
def mkfig():
    plt.clf(); #plt.ion(); #plt.show()
    plt.gca().set_xlim(x0,x1)
    plt.gca().set_ylim(y0,y1)

    # Add elements
    cbmax=150
    im=plt.pcolormesh(glamfe[:-1,:-1],gphife[:-1,:-1],bathy)
    cb=plt.colorbar(im);
    #im.set_cmap(cmocean.cm.haline);
    im.set_cmap('winter_r');
    im.set_clim([0,cbmax]); cb.set_clim(0,cbmax);

    gridlines(glamfe,gphife,'k')
    drawrivers()
    ax = plt.gca()
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    plt.pause(0.1)

def mkfig_click():
    mkfig()
    global cid
    cid = fig.canvas.mpl_connect('button_press_event', onclick)



# Grid limits
fig=plt.figure(figsize=(24,12))
x0,x1=-123.2, -122.15
y0,y1=47,48.1
NY, NX = glamfe.shape[0], glamfe.shape[1]



# Snohomish river mouth
Fill = [(309,145),(309,144),(310,148),(310,147),(310,146),(310,145),(310,144),(310,143),(310,142),(311,148),(311,147),(311,146),(311,145),(311,144),(311,143),(311,142),(311,141),(312,149),(312,148),(312,147),(312,146),(312,145),(312,144),(312,143),(312,142),(312,141),(312,140),(312,138),(312,137),(313,149),(313,148),(313,147),(313,146),(313,145),(313,144),(313,143),(313,142),(313,141),(313,140),(313,138),(313,137),(314,146),(314,145),(314,144),(314,143),(314,142),(314,141),(314,139),(315,144),(315,143),(315,142),(315,140),(315,139),(316,144),(316,143),(316,141),(316,140),(317,144),(317,143),(317,142),(317,141),(313,139),(314,140),(315,141),(316,142),(318,145),(318,144),(318,143),(318,142),(319,145),(319,144),(319,143),(319,142),(320,145),(320,144),(320,143),(321,143),(318,141)]
#Dabob bay 
Fill += [(198,145),(198,146),(197,146),(199,146),(199,147),(198,147),(197,147),(199,148),(200,148),(198,144),(208,149),(209,149),(207,149),(236,111),(237,111),(180,127),(180,128),(181,128)]
# Near Brinnon and Seabeck
Fill += [(168,119),(167,119),(177,124),(178,124),(179,124),(179,125),(180,125),(180,126),(179,126),(181,112)]
# North Bay
Fill += [(150,59),(150,60),(151,60),(151,59),(152,61),(152,62),(151,61),(153,61),(153,62),(153,63),(153,64),(154,64)]

# Simple fill in with 4m depth here
for i,j in Fill:
    bathy[j,i]=4


# Down towards Skokomish
Fill2 = [(135,100),(134,98),(133,98),(132,97),(131,97),(131,96),(130,96),(129,96),(127,96),(127,95),(126,95),(125,94),(124,94),(124,93),(123,93),(122,93),(122,92),(121,92),(120,92),(119,91),(118,90),(117,90),(116,89),(115,88),(114,87),(113,87),(112,86),(111,85),(107,77),(109,76),(108,76),(109,75),(106,76),(105,76),(104,76),(103,76),(102,75),(103,75),(104,75),(105,75),(104,74),(103,74),(103,73),(112,72),(114,71),(116,70),(121,70),(122,70),(148,73),(149,73),(150,73),(151,73),(152,73),(153,73),(149,72),(148,71),(150,72),(151,72),(136,73),(107,76)]
# Southern most edge of Puget Sound (Budd Inlet, Oyster Bay, etc)
Fill2 += [(90,13),(90,14),(89,14),(88,14),(88,13),(89,13),(87,13),(86,13),(86,12),(87,14),(85,11),(85,12),(87,23),(87,24),(88,24),(86,24),(86,25),(85,25),(85,24),(84,24),(84,25),(83,25),(83,24),(82,24),(81,24),(80,24),(119,13),(119,12),(119,11),(118,11),(118,10),(100,3),(101,3),(102,4),(102,5),(131,1),(132,1),(132,2),(133,1),(132,0),(133,0),(134,0),(131,0),(154,12),(135,37),(136,37),(136,38),(137,38),(138,38),(118,37),(118,36),(116,38),(118,38)]

# These ones we'll fill with max(4, mean of nearby values)
def chk(bathy,ijlist):
    for i,j in ijlist:
        if bathy[j,i] is ma.masked: return True
    return False
while(chk(bathy,Fill2)):
    print("Filling pass")
    bathy0=ma.masked_array(np.copy(bathy),mask=bathy.mask)
    ny,nx=bathy.shape
    for i,j in Fill2:
        i1,i2=max(0,i-1),min(nx-1,i+1)
        j1,j2=max(0,j-1),min(ny-1,j+1)
        tmp = bathy0[j1:j2+1,i1:i2+1]
        nzvals=tmp[~tmp.mask].data
        bathy[j,i]=max(4, np.mean(nzvals))


if False:
    # Interactive mode, click boxes to fill and accumulate a list of coordinates
    mkfig_click()

else:
    # Plot the result of filling
    mkfig()
    for i,j in Fill+Fill2:
        c='k'
        lw=1
        #plt.plot(glamt[j,i],gphit[j,i],'k*')  # mark with stars
        # Box them in
        plt.plot(glamfe[j,i:i+2],gphife[j,i:i+2],c,linewidth=lw)
        plt.plot(glamfe[j+1,i:i+2],gphife[j+1,i:i+2],c,linewidth=lw)
        plt.plot(glamfe[j:j+2,i],gphife[j:j+2,i],c,linewidth=lw)
        plt.plot(glamfe[j:j+2,i+1],gphife[j:j+2,i+1],c,linewidth=lw)
    plt.tight_layout(pad=0.5)
    plt.savefig('fix_puget.png')

embed()
