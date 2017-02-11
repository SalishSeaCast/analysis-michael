from IPython import embed
import numpy as np
import netCDF4 as nc

def t2u(lont,latt):
    # We don't have a value for the rightmost u points
    lonu = np.zeros(lont.shape)
    latu = np.zeros(lont.shape)
    latu[:,0:-1]=0.5*(latt[:,0:-1] + latt[:,1:])
    lonu[:,0:-1]=0.5*(lont[:,0:-1] + lont[:,1:])
    return lonu,latu


def t2v(lont,latt):
    # We don't have a value for the northmost v points
    lonv = np.zeros(lont.shape)
    latv = np.zeros(lont.shape)
    latv[0:-1,:]=0.5*(latt[0:-1,:] + latt[1:,:])
    lonv[0:-1,:]=0.5*(lont[0:-1,:] + lont[1:,:])
    return lonv,latv


def t2f(lont,latt):
    # We don't have values in the rightmost u pt or northmost v pt
    lonf = np.zeros(lont.shape)
    latf = np.zeros(lont.shape)
    lonf[0:-1,0:-1] = 0.25*(lont[0:-1,0:-1] + lont[1:,0:-1] + lont[0:-1,1:] + lont[1:,1:])
    latf[0:-1,0:-1] = 0.25*(latt[0:-1,0:-1] + latt[1:,0:-1] + latt[0:-1,1:] + latt[1:,1:])
    return lonf,latf
    
    
def haversine(lon1, lat1, lon2, lat2):
    """ This is copied from salishsea_tools """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c # from ss tools
    #km = 6371 * c  # from earthdist.m
    return km*1000


def gete1(lon, lat, expandleft=False):
    if expandleft:
        lon,lat = expandi(lon, lat)
    dx = np.zeros(lon.shape)
    lon1 = lon[0:-1, 0:-1]
    lat1 = lat[0:-1, 0:-1]
    lon2 = lon[0:-1, 1:  ]
    lat2 = lat[0:-1, 1:  ]
    dx[:-1,:-1] = haversine(lon1,lat1,lon2,lat2)
    return dx


def gete2(lon, lat, expanddown=False):
    if expanddown:
        lon,lat = expandj(lon, lat)
    dy = np.zeros(lon.shape)
    lon1 = lon[0:-1, 0:-1]
    lat1 = lat[0:-1, 0:-1]
    lon2 = lon[1:  , 0:-1]
    lat2 = lat[1:  , 0:-1]
    dy[:-1,:-1] = haversine(lon1,lat1,lon2,lat2)
    return dy


def expandi(lon, lat):
    # Expands the grid one to the left by linear extrapolation
    NY, NX = lon.shape
    lon2, lat2 = np.zeros([NY, NX+1]), np.zeros([NY, NX+1])
    # Lon
    lon2[:,1:] = lon
    lon2[:,0]  = lon[:,0] - (lon[:,1] - lon[:,0])
    # Lat
    lat2[:,1:] = lat
    lat2[:,0]  = lat[:,0] - (lat[:,1] - lat[:,0])
    return lon2, lat2


def expandj(lon, lat):
    # Expands the grid one down by linear extrapolation
    NY, NX = lon.shape
    lon2, lat2 = np.zeros([NY+1, NX]), np.zeros([NY+1, NX])
    # Long
    lon2[1:,:] = lon
    lon2[0,:]  = lon[0,:] - (lon[1,:] - lon[0,:])
    # Lat
    lat2[1:,:] = lat
    lat2[0,:]  = lat[0,:] - (lat[1,:] - lat[0,:])
    return lon2, lat2


def writecoords(fname,
                glamt,glamu,glamv,glamf,
                gphit,gphiu,gphiv,gphif,
                e1t,e1u,e1v,e1f,
                e2t,e2u,e2v,e2f):
    # Build a NEMO format coordinates file

    cnc = nc.Dataset(fname, 'w', clobber=True)
    NY,NX = glamt.shape
    
    # Create the dimensions
    cnc.createDimension('x', NX)
    cnc.createDimension('y', NY)
    cnc.createDimension('time', None)
   
    # Create the float variables
    cnc.createVariable('nav_lon', 'f', ('y','x'), zlib=True, complevel=4)
    cnc.variables['nav_lon'].setncattr('units','degrees_east')
    cnc.variables['nav_lon'].setncattr('comment','at t points')
    
    cnc.createVariable('nav_lat', 'f', ('y','x'), zlib=True, complevel=4)
    cnc.variables['nav_lat'].setncattr('units','degrees_north')
    cnc.variables['nav_lat'].setncattr('comment','at t points')

    cnc.createVariable('time', 'f', ('time'), zlib=True, complevel=4)
    cnc.variables['time'].setncattr('units','seconds since 0001-01-01 00:00:00')
    cnc.variables['time'].setncattr('time_origin','0000-JAN-01 00:00:00')
    cnc.variables['time'].setncattr('calendar','gregorian')

    # Create the double variables
    varlist = ['glamt','glamu','glamv','glamf','gphit','gphiu','gphiv','gphif']
    varlist+= ['e1t','e1u','e1v','e1f','e2t','e2u','e2v','e2f']
    for v in varlist:
        cnc.createVariable(v, 'd', ('time','y','x'), zlib=True, complevel=4)
        cnc.variables[v].setncattr('missing_value',1e+20)

    # Write the data
    cnc.variables['nav_lon'][...] = glamt
    cnc.variables['nav_lat'][...] = gphit
    #
    cnc.variables['glamt'][0,...] = glamt
    cnc.variables['glamu'][0,...] = glamu
    cnc.variables['glamv'][0,...] = glamv
    cnc.variables['glamf'][0,...] = glamf
    #
    cnc.variables['gphit'][0,...] = gphit
    cnc.variables['gphiu'][0,...] = gphiu
    cnc.variables['gphiv'][0,...] = gphiv
    cnc.variables['gphif'][0,...] = gphif
    #
    cnc.variables['e1t'][0,...] = e1t
    cnc.variables['e1u'][0,...] = e1u
    cnc.variables['e1v'][0,...] = e1v
    cnc.variables['e1f'][0,...] = e1f
    #
    cnc.variables['e2t'][0,...] = e2t
    cnc.variables['e2u'][0,...] = e2u
    cnc.variables['e2v'][0,...] = e2v
    cnc.variables['e2f'][0,...] = e2f
                    
    cnc.close()
