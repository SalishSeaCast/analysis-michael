import netCDF4 as nc
import numpy as np
import shutil, datetime, sys

w0='weights-gem2.5-ops.nc'         # Original weights file (NC4 with attributes)
w1='met_gem_weight.nc'             # New weights from get_nemo_weights (NC3, no attributes)
w2='weights-gem2.5-ops_201702.nc'  # Final weights file

# Copy the documented NC4 weights file
shutil.copyfile(w0,w2)

with nc.Dataset(w1, 'r') as f1, nc.Dataset(w2, 'r+') as f2:
    # Replace the weights with the new weights
    for v in ['src01','src02','src03','src04','wgt01','wgt02','wgt03','wgt04']:
        f2.variables[v][:] = f1.variables[v][:]
    
    # Update notes
    f2.history += ("\n [{}] Adjust weights for 201702 coordinates."
                   .format(datetime.datetime.today().strftime('%Y-%m-%d')))
