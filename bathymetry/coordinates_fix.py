from IPython import embed
import netCDF4 as nc
import numpy as np
import datetime

def coordinates_fix(infile, outfile):
    """
    Flip e1/e2 in coordinates.nc
    """
    fin  = nc.Dataset(infile, 'r')
    fout = nc.Dataset(outfile, 'w')

    # Copy global attributes
    for attr in fin.ncattrs():
        fout.setncattr(attr, fin.getncattr(attr))

    # Copy dimensions
    dim = fin.dimensions['x']
    fout.createDimension(dim.name, dim.size)
    dim = fin.dimensions['y']
    fout.createDimension(dim.name, dim.size)
    dim = fin.dimensions['time']
    fout.createDimension(dim.name, None)

    # Create identical output variables and copy attributes
    for k, v in fin.variables.items():
        fout.createVariable(v.name, v.datatype, v.dimensions, zlib=True, complevel=4, shuffle=False)
        for attr in v.ncattrs():
            fout.variables[v.name].setncattr(attr, v.getncattr(attr))

    # Copy variable data
    for k, v in fin.variables.items():
        if v.name.startswith('e'):
            # Swap e1,e2 pairs
            if v.name[1] == '1': dst ='e2'+v.name[2]
            if v.name[1] == '2': dst ='e1'+v.name[2]
            print("Writing {} to {} ...".format(v.name,dst))
            fout.variables[dst][:] = fin.variables[v.name][:]
        else:
            # Simply copy the data
            print("Copying var {} ...".format(v.name))
            fout.variables[v.name][:] = fin.variables[v.name][:]

    # Add note to history
    note ='[{}] Swapped e1/e2 variables in the original '
    note+='coordinates_seagrid_SalishSea.nc and added compression at level 4.'
    fout.setncattr('history', note.format(datetime.datetime.today().strftime('%Y-%m-%d')))

    fin.close()
    fout.close()

if __name__ == "__main__":
    path='/home/mdunphy/MEOPAR/NEMO-forcing/grid/'
    coordinates_fix(path+'coordinates_seagrid_SalishSea.nc',
                    path+'coordinates_seagrid_SalishSea2.nc')
