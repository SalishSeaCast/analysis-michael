import netCDF4 as nc
import numpy as np
from collections import OrderedDict
from resource import getrlimit, RLIMIT_NOFILE

class scDataset(object):
    def __init__(self, files):
        """
        Simple Concatenated Dataset

        scDataset is a partial reimplementation of MFDataset that automates the
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

        Arguments: files is a list of netcdf files

        Example usage:

        # Get the (concatenated) output times
        with scDataset(files) as ds:
            t = ds.variables['time_counter'][:]

        # Get temperature at all times and all depth at one location
        with scDataset(files) as ds:
            temper = ds.variables['votemper'][:,:,100,100]

        # Load surface salinity at each time in a loop for plotting/animation
        with scDataset(files) as ds:
            for ti in range(ds.variables['vosaline'].shape[0]):
                print("Loading time "+str(ti))
                surfsal = ds.variables['vosaline'][ti,0,:,:]
                # make_a_plot(surfsal)

        # Demo to show that normal Python indexing and slicing works
        with scDataset(files) as ds:
            t1 = ds.variables['votemper'][29:33:-1,-10:-1,100:130]
            print(t1.shape)

        """
        # Initialize a dataset manager with the list of files
        self._dsmgr = _scDatasetManager(files)

        # Open the first dataset and set a few class variables
        d0 = self._dsmgr[0]
        self.description = d0.description
        self.file_format = d0.file_format
        self.filepath    = files

        # Find the time dimension name
        for dim in d0.dimensions:
            if d0.dimensions[dim].isunlimited():
                timedimname = dim
                break

        # Open each dataset, get time dimension size and set the indices fi and li
        fi = []  # file (dataset) index
        li = []  # local time index
        for di in range(len(files)):
            curlen = self._dsmgr[di].dimensions[timedimname].size
            fi += [di for x in range(curlen)]
            li += [x for x in range(curlen)]

        # First dimension must be unlimited, else use the first dataset
        self.variables = OrderedDict()
        vars0 = d0.variables
        for vname in vars0:
            if vars0[vname].dimensions[0] == timedimname:
                # We concatenate this variable
                self.variables[vname] = scVariable(vars0[vname], vname, self._dsmgr, fi, li)
            else:
                # Passthrough this variable to the first file
                self.variables[vname] = vars0[vname]

    def close(self):
        self._dsmgr.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __del__(self):
        self.close()


class _scDatasetManager(object):
    """
    Manages datasets by opening/closing them on demand
    """

    def __init__(self, files):
        self._files   = files
        self._MAXOPEN = getrlimit(RLIMIT_NOFILE)[0] // 5
        self._dslist  = [(-1, None)] * self._MAXOPEN

    def __getitem__(self, di):
        """
        Input is an integer index di, we return the corresponding nc.Dataset
        Also we cache open datasets to economize on opening/closing them
        """
        # Compute slot index (circular buffer; and index 0 is always kept open)
        si = 1 + (di - 1) % (self._MAXOPEN - 1) if di > 0 else 0

        # Check what is currently stored in slot si
        ci, ds = self._dslist[si]
        if ci != di:
            if ds is not None:
                # Repurpose slot si for the requested dataset
                ds.close()
            # Now open the requested dataset and store it in slot si
            ds = nc.Dataset(self._files[di], 'r')
            self._dslist[si] = (di, ds)
        return ds

    def close(self):
        for di, ds in self._dslist:
            if ds is not None:
                ds.close()
        self._dslist = []


class scVariable(object):
    """
    Builds a concatenated version of a netCDF Variable type
     - We aim to have correct indexing, and set a few class variables such as
       shape and dimensions correctly. Attribute handling, etc is not implemented.
    """
    def __init__(self, v0, vname, datasets, fi, li):
        self.ds = datasets
        self._fi = fi
        self._li = li

        # Set a few class variables
        self.name       = vname
        self.dimensions = v0.dimensions
        self.dtype      = v0.dtype
        self.ndim       = v0.ndim
        self.shape      = (len(self._fi), ) + v0.shape[1:]

    def __getitem__(self, initems):
        """
        Implement Python indexing: int, slice, ellipsis accepted
        """
        # Make the input iterable
        if not isinstance(initems, tuple):
            initems = initems,

        # Convert any ellipsis to slice
        items = [slice(None,None,None)]*self.ndim
        for i, item in enumerate(initems):
            if item is not Ellipsis:
                items[i] = item
            else:
                for j, item in enumerate(reversed(initems)):
                    if item is not Ellipsis:
                        items[self.ndim-j-1] = item
                    else:
                        break
                break

        # Find the time indices
        ti = items[0]      # global time indices to extract, may be int or slice
        fi = self._fi[ti]  # index of each file (dataset) to draw from
        li = self._li[ti]  # local time index for each dataset

        # For single time output (no concatenation), just draw from the right dataset
        if type(ti) is int or type(ti) is np.int64:
            if self.ndim == 1:
                out = self.ds[fi].variables[self.name][li]
            if self.ndim == 2:
                out = self.ds[fi].variables[self.name][li, items[1]]
            if self.ndim == 3:
                out = self.ds[fi].variables[self.name][li, items[1], items[2]]
            if self.ndim == 4:
                out = self.ds[fi].variables[self.name][li, items[1], items[2], items[3]]
            return out

        # If we need to concatenate, then we need to determine the output
        # array size. This approach is an ugly hack but it works.
        sizo = [1] * self.ndim  # assume one in each dimension
        rdim = []               # list of dimensions to remove
        for ii, item in enumerate(items):
            if type(item) is int or type(item) is np.int64:
                rdim += [ii]
            else:                             # update output size at this dim if not an integer index
                tmp = [None] * self.shape[ii] # build a dummy array
                sizo[ii] = len(tmp[item])     # index the dummy array, record length
        out = np.zeros(sizo, self.dtype)      # allocate output array with matching data type
        out = np.squeeze(out, axis=tuple(rdim))  # remove unwanted singleton dimensions

        # Now we read each time index sequentially and fill the output array
        for ii in range(len(fi)):
            if self.ndim == 1:
                out[ii] = self.ds[fi[ii]].variables[self.name][li[ii]]
            if self.ndim == 2:
                out[ii, ...] = self.ds[fi[ii]].variables[self.name][li[ii], items[1]]
            if self.ndim == 3:
                out[ii, ...] = self.ds[fi[ii]].variables[self.name][li[ii], items[1], items[2]]
            if self.ndim == 4:
                out[ii, ...] = self.ds[fi[ii]].variables[self.name][li[ii], items[1], items[2], items[3]]
        return out
