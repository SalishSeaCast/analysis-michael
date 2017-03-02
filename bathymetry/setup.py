# To build, run
#  python setup.py build_ext --inplace

from distutils.core import setup, Extension
import numpy

# Define the extension module
slidingmean_c = Extension('slidingmean_c', sources=['slidingmean_c.c'],
                          include_dirs=[numpy.get_include()])

# run the setup
setup(ext_modules=[slidingmean_c])
