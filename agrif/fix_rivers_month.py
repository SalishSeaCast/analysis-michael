#!/usr/bin/env python
import netCDF4 as nc
import numpy as np
import sys

def fix_rodepth(infile):
    """
    Simple script to enforce rodepth = 1 at missing values
    """
    with nc.Dataset(infile, 'r+') as f:
        rodepth= f.variables['rodepth'][:]
        idx = rodepth == 0
        if np.any(idx):
            print("Switching zeros to 1")
            rodepth[idx] = 1
            f.variables['rodepth'][:] = rodepth
    
if __name__ == "__main__":
    fix_rodepth(sys.argv[1])
