# This is a script version of this notebook:
#  https://nbviewer.jupyter.org/urls/bitbucket.org/salishsea/tools/raw/tip/I_ForcingFiles/Atmos/ImproveWeightsFile.ipynb

from datetime import datetime
import os
import netCDF4 as nc
import numpy as np
from salishsea_tools import nc_tools

met_gem_weight = 'met_gem_weight.nc'           # Output from the get_weight_nemo tool
netcdf4_weight = 'wcvi-weights-gem2.5-ops.nc'  # Filename for the improved weights file
atmos_grid_name = 'GEM 2.5km Operational'

history = ('[{}] Converted to netCDF4 zlib=True dataset, and added CF-1.6 metadata.'
            .format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

mtime = datetime.fromtimestamp(os.path.getmtime(met_gem_weight))
history += '\n [{}] Created netCDF3 classic dataset with get_weight_nemo.'.format(mtime)


with nc.Dataset(met_gem_weight,'r') as src, nc.Dataset(netcdf4_weight, 'w') as weights:

    y_size, x_size = src.variables['src01'].shape
    weights.createDimension('x', x_size)
    weights.createDimension('y', y_size)
    lat_size, lon_size = src.variables['nav_lon'].shape
    weights.createDimension('lon', lon_size)
    weights.createDimension('lat', lat_size)
    weights.createDimension('numwgt', 4)

    lats = weights.createVariable('nav_lat', float, ('lat', 'lon'), zlib=True)
    lons = weights.createVariable('nav_lon', float, ('lat', 'lon'), zlib=True)
    src01 = weights.createVariable('src01', int, ('y', 'x'), zlib=True)
    wgt01 = weights.createVariable('wgt01', float, ('y', 'x'), zlib=True)
    src02 = weights.createVariable('src02', int, ('y', 'x'), zlib=True)
    wgt02 = weights.createVariable('wgt02', float, ('y', 'x'), zlib=True)
    src03 = weights.createVariable('src03', int, ('y', 'x'), zlib=True)
    wgt03 = weights.createVariable('wgt03', float, ('y', 'x'), zlib=True)
    src04 = weights.createVariable('src04', int, ('y', 'x'), zlib=True)
    wgt04 = weights.createVariable('wgt04', float, ('y', 'x'), zlib=True)

    lats[:] = src.variables['nav_lat'][:]
    lats.units = 'degrees_north'
    lats.long_name = 'Latitude'
    lats.valid_range = (np.min(lats[:]), np.max(lats[:]))
        
    lons[:] = src.variables['nav_lon'][:]
    lons.units = 'degrees_east'
    lons.long_name = 'Longitude'
    lons.valid_range = (np.min(lons[:]), np.max(lons[:]))

    vars = ((src01, wgt01), (src02, wgt02), (src03, wgt03), (src04, wgt04))
    for i, sw in enumerate(vars):
        s, w = sw
        sname = 'src{:02d}'.format(i+1)
        wname = 'wgt{:02d}'.format(i+1)
        s[:] = src.variables[sname][:]
        s.units = 1
        s.long_name = '{} Grid Index {} (Flattened)'.format(atmos_grid_name,i+1)
        s.valid_range = np.array(
            (np.min(src.variables[sname]), np.max(src.variables[sname])))
        
        w[:] = src.variables[wname][:]
        w.units = 1
        w.long_name = 'Salish Sea Grid Weights for {}'.format(sname)
        w.valid_range = np.array(
            (np.min(src.variables[wname]), np.max(src.variables[wname])))
    
    nc_tools.init_dataset_attrs(
        weights,
        'West Coast of Vancouver Island NEMO {} Atmospheric Forcing Interpolation Weights'.format(atmos_grid_name),
        [],
        netcdf4_weight,
    )

    weights.history = history
    weights.source = "https://bitbucket.org/salishsea/analysis-michael/src/tip/weights/improveweights.py"
