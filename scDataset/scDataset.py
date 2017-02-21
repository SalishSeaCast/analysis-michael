import netCDF4 as nc
import numpy as np
from IPython import embed
from collections import OrderedDict

class scDataset(object):
    def __init__(self, files):
        """
        !! EXPERIMENTAL !!

        Simple Concatenated Dataset

        scDataset is a partial implementation of MFDataset that automates the
        concatenation of netCDF variables split between multiple files. The aim
        is to provide a simple concatenated interface to variables where the
        first dimension is the unlimited dimension (usually "time_counter").

        Variables where the first dimension is not the unlimited time dimension
        (such as lon, lat, etc) are not concatenated. They are, however, made
        available in a "pass through" sense to the first dataset in the list.
        Thus, you can read those variables without opening a second Dataset class.

        The other MFDataset features are not implemented: attributes, etc are
        ignored for the concatenated variables. It may be possible to add the
        those features but the goal here is simply automated concatenation.

        Building this class was motivated by deficiencies in the other options
        for split-file concatenation:
         - xarray.open_mfdataset() loads the entire dataset into memory, which
           is both very slow and memory intensive.
         - netCDF4.MFDataset refuses to open NETCDF4 format files
        In the event that netCDF4.MFDataset is improved to work with NETCDF4
        files, this will become obsolete.

        !! EXPERIMENTAL !!

        Arguments: files is a list of netcdf files

        Example usage:

        # Get the (concatenated) output times
        ds = scDataset(files)
        t = ds.variables['time_counter'][:]
        ds.close()

        # Get temperature at all times and all depth at one location
        ds = scDataset(files)
        temper = ds.variables['votemper'][:,:,100,100]
        ds.close()

        # Load surface salinity at each time in a loop for plotting/animation
        ds = scDataset(files)
        for ti in range(ds.variables['vosaline'].shape[0]):
            print("Loading time "+str(ti))
            surfsal = ds.variables['vosaline'][ti,0,:,:]
            # make_a_plot(surfsal)
        ds.close()

        # Demo to show that normal Python indexing and slicing works
        ds = scDataset(files)
        t1 = ds.variables['votemper'][29:33:-1,-10:-1,100:130]
        print(t1.shape)
        ds.close()

        !! EXPERIMENTAL !!

        """
        # Open all of the datasets        
        self._datasets=[nc.Dataset(f,'r') for f in files]

        # Set a few class variables
        d0=self._datasets[0]
        self.description = d0.description
        self.file_format = d0.file_format
        self.filepath    = files

        # Find the time dimension name
        for dim in d0.dimensions.keys():
            if d0.dimensions[dim].isunlimited():
                timename = dim
                break

        # Get the indices
        fi = []   # file (dataset) index
        li = []   # local time index
        for di,d in enumerate(self._datasets):
            curlen = d.variables[timename].shape[0]
            fi += [di for x in range(curlen)]
            li += [x for x in range(curlen)]

        # The first dimension must be the unlimited dimension, else we
        # pass through to the first dataset
        self.variables = OrderedDict()
        for vname in d0.variables.keys():
            if d0.variables[vname].dimensions[0] == timename:
                # We concatenate this variable
                varlist = [ds.variables[vname] for ds in self._datasets]
                self.variables[vname] = scVariable(varlist, fi, li)
            else:
                # Passthrough this variable to the first file
                self.variables[vname] = d0.variables[vname]


    def close(self):
        """
        Closes the datasets. This may be useful when opening many datasets to
        keep the number of open files under control.
        """
        for ds in self._datasets:
            ds.close()
        self._datasets,self._vars,self.shape=[],[],[]

    def __del__(self):
        self.close()


class scVariable(object):
    """
    Builds a concatenated version of several netCDF Variable types
     - We aim to have correct indexing, and set a few class variables
       such as shape correctly. Attribute handling, etc is not implemented.
    """
    def __init__(self, varlist, fi, li):
        self._vars = varlist
        self._fi = fi
        self._li = li

        # Set a few class variables
        v0 = varlist[0]
        self.name       = v0.name
        self.dimensions = v0.dimensions
        self.dtype      = v0.dtype
        self.ndim       = v0.ndim
        self.size       = np.sum([v.size for v in varlist])
        self.shape      = (len(self._fi),) + v0.shape[1:]

    def __getitem__(self, items):
        """
        Implement Python indexing: int or slice accepted
        """
        # Make the input iterable
        if not isinstance(items, tuple):
            items=[items]

        # Check number of dimensions
        ndim = len(items)
        if self.ndim != ndim: print("Number of dimensions mismatch")

        # Find the time indices
        ti = items[0]      # global time indices to extract, may be int or slice
        fi = self._fi[ti]  # index of each file (dataset) to draw from
        li = self._li[ti]  # local time index for this dataset

        # For single time output (no concatenation), just draw from the right dataset
        if type(ti) is int:
            if ndim==1:
                out = self._vars[fi][li]
            if ndim==2:
                out = self._vars[fi][li, items[1]]
            if ndim==3:
                out = self._vars[fi][li, items[1], items[2]]
            if ndim==4:
                out = self._vars[fi][li, items[1], items[2], items[3]]
            return out

        # If we need to concatenate, then we need to determine the output
        # array size. This approach is an ugly hack but it works.
        sizo = [1]*ndim  # assume one in each dimension
        for ii,item in enumerate(items):
            if type(item) is not int:  # update output size at this dim if not an integer index
                tmp = [None]*self.shape[ii]   # build a dummy array
                sizo[ii] = len(tmp[item])     # index the dummy array, record length
        out = np.squeeze(np.zeros(sizo,self.dtype)) # remove singleton dimensions

        # Now we read each time index sequentially and fill the output array
        for ii in range(len(fi)):
            if ndim==1:
                out[ii]     = self._vars[ fi[ii] ][li[ii]]
            if ndim==2:
                out[ii,...] = self._vars[ fi[ii] ][li[ii], items[1]]
            if ndim==3:
                out[ii,...] = self._vars[ fi[ii] ][li[ii], items[1], items[2]]
            if ndim==4:
                out[ii,...] = self._vars[ fi[ii] ][li[ii], items[1], items[2], items[3]]
        return out
