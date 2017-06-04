# This script takes a restart file from the results archive and produces an "expanded" version, where
# data is extrapolated in an ad-hoc manner into the land cells. The purpose is to use these expanded files
# with the AGRIF nesting tools, because they don't extrapolate well, and proceeding without "expanded" inputs
# usually means the AGRIF domain will crash due to bad values in the deep cells

from multiprocessing import Pool
from IPython import embed
import netCDF4 as nc
import numpy.ma as ma
import numpy as np
import os,sys,time

def vfill(data,ndeep):
    # input: 3d array
    k0,j0,i0 = np.where(data==0.0)
    n0 = i0.size
    for i in range(data.shape[2]):
        for j in range(data.shape[1]):
            deepen = ndeep
            for k in range(1,data.shape[0]):
                if data[k,j,i] == 0:
                    data[k,j,i] = data[k-1,j,i]
                    deepen -= 1
                if deepen < 1:
                    continue
    k0,j0,i0 = np.where(data==0.0)
    n1 = i0.size
    print("vfill: ndeep = {}, filled {:6d} missing points ({} to {})".format(ndeep,n0-n1,n0,n1))
    return data

def hfill(pn,data,k):
    ny,nx=data.shape[0],data.shape[1]
    data2 = np.copy(data)
    # Find missing values
    missing = data == 0.0
    if missing.any():
        # Convert to indices
        j0,i0 = np.where(missing)
        n0 = i0.size

        # Fill points
        for ci,cj in zip(i0,j0):
            if ci>0 and ci < nx-1 and cj>0 and cj < ny-1:   # not at the edges
                tmp = data[cj-1:cj+2,ci-1:ci+2]
            else:                                           # edges
                x1,x2=max(0,ci-1),min(nx-1,ci+1)
                y1,y2=max(0,cj-1),min(ny-1,cj+1)
                tmp = data[y1:y2+1,x1:x2+1]

            nnz = np.count_nonzero(tmp)
            if nnz>1:
                nonzeros = tmp != 0.0
                #nz = nonzeros.size - nnz
                #if nnz > nz:
                data2[cj,ci] = np.mean(tmp[nonzeros])
    
    j0,i0 = np.where(data2==0.0)
    n1 = i0.size
    print("hfill: pass = {}, level = {:2d}, filled {:6d} missing points ({} to {})".format(pn,k,n0-n1,n0,n1))    
    return (data2,k)

def rstexpand(infile, outfile, npass, ndeep, nrepeat):
    """
    Hack to fix the restart files...
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
    try:
        dim = fin.dimensions['z']
        fout.createDimension(dim.name, dim.size)
    except:
        dim = fin.dimensions['deptht']
        fout.createDimension(dim.name, dim.size)
    try:
        dim = fin.dimensions['t']
        fout.createDimension(dim.name, None)
    except:
        dim = fin.dimensions['time_counter']
        fout.createDimension(dim.name, None)

    # Copy variables
    for k, v in fin.variables.items():
        #if 'TRN' in v.name:
        #    continue
        # Create identical output variable, copy attributes
        fout.createVariable(v.name, v.datatype, v.dimensions, zlib=True, complevel=4, shuffle=False)
        for attr in v.ncattrs():
            fout.variables[v.name].setncattr(attr, v.getncattr(attr))

        ndim = len(fin.variables[v.name].shape)

        if ndim < 3 or v.name in ['rnf_b','rnf_hc_b','rnf_sc_b']:
            # Just copy the data
            print("Copying var {} ...".format(v.name))
            fout.variables[v.name][:] = fin.variables[v.name][:]
            
        elif ndim == 3: # t,y,x
            print("Fixing 2D var {} ...".format(v.name))
            data = fin.variables[v.name][:]
            for pn in range(npass):  # npasses
                t0=time.time()
                data2 = np.zeros(data.shape)
                data2[0,...] = hfill(pn+1,data[0,...],0)[0]
                data = data2
                print("processed in {:.2f} s".format(time.time()-t0))
            fout.variables[v.name][:]=data            

        elif ndim == 4: # t,z,y,x
            print("Fixing 3D var {} ...".format(v.name))
            data = fin.variables[v.name][:]
            if ma.isMaskedArray(data): data = data.data
            nz=data.shape[1]

            for repeat in range(nrepeat):
                # Do npass passes of horizontal filling
                for pn in range(npass):
                    t0=time.time()
                    data2 = np.zeros(data.shape)
                    with Pool(processes=8) as pool:  # start worker processes
                        def saveresult(x):
                            data2[0,x[1],...]=x[0]
                        for k in range(nz): # dispatch a process for each layer
                            pool.apply_async(hfill, args=(pn+1,data[0,k,...],k),callback=saveresult)
                        pool.close()
                        pool.join()
                    data = data2
                    print("hfill processed in {:.2f} s".format(time.time()-t0))
    
                # Vertical fill ndeep cells
                t0=time.time()
                data[0,...] = vfill(data[0,...],ndeep)
                print("vfill processed in {:.2f} s".format(time.time()-t0))
            # Store fixed data
            fout.variables[v.name][:]=data

    fin.close()
    fout.close()

if __name__ == "__main__":
    """
    Usage: python expand_restart.py /results/SalishSea/PICKONE/ddmmmyy/SalishSea_NNNNNNNN_restart.nc
    The expanded file will be written to $HOME/MEOPAR/NEMO-forcing-MD/PICKONE/ddmmmyy/SalishSea_NNNNNNNN_restart.nc

    Acquire a list of restarts by:
      find /results/SalishSea/hindcast -name "*restart*.nc" | sort -k2 -t _

    outdir = os.path.dirname(outfile)

    """
    npass = 2   # Number of horizontal fill passes per repeat
    ndeep = 2   # Number of vertical cells to fill per repeat
    nrepeat = 2 # number of times to run a horizontal+vertical fill

    # Input file
    infile = sys.argv[1]
    inprefix = '/results/SalishSea/'
    if not infile.startswith(inprefix) or not infile.endswith('.nc'):
        print('Please select a netcdf restart file under {} to expand (use full path)'.format(inprefix))
        exit()
    
    # Output file
    outprefix = os.path.expandvars('$HOME/MEOPAR/NEMO-forcing-MD/')
    outsuffix = os.path.relpath(infile,inprefix)
    outfile = os.path.join(outprefix, outsuffix)
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    rstexpand(infile, outfile, npass, ndeep, nrepeat)

