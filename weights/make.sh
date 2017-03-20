#!/bin/bash
#- From Fury...
# mpif90 -c grid.f90 -I/export/opt/netcdf/3.6.3/intel/medium/include -L/export/opt/netcdf/3.6.3/intel/medium/lib -lnetcdf 
# mpif90 -c map.F90 -I/export/opt/netcdf/3.6.3/intel/medium/include -L/export/opt/netcdf/3.6.3/intel/medium/lib -lnetcdf 
# mpif90 -c get_weight_nemo.F90 -I/export/opt/netcdf/3.6.3/intel/medium/include -L/export/opt/netcdf/3.6.3/intel/medium/lib -lnetcdf 
# mpif90 get_weight_nemo.o map.o grid.o -I/export/opt/netcdf/3.6.3/intel/medium/include -L/export/opt/netcdf/3.6.3/intel/medium/lib -lnetcdf 

#- On Ace-net
#LIBNETCDF=/usr/local/netcdf.intel-3.6.3
#mpif90 -c grid.f90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
#mpif90 -c map.F90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
#mpif90 -c get_weight_nemo.F90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
#mpif90 get_weight_nemo.o map.o grid.o -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf


#- On Mahone (ace-net)
#LIBNETCDF=/usr/local/netcdf.intel-3.6.3
#ifort -c grid.f90 -I ${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
#ifort -c map.F90 -I ${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
#ifort -c get_weight_nemo.F90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
#ifort get_weight_nemo.o map.o grid.o -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf

#- On Anchor (Dalhousie)
LIBNETCDF=/opt/netcdf-g95/                  
g95 -c grid.f90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
g95 -c map.F90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
g95 -c get_weight_nemo.F90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
g95 get_weight_nemo.o map.o grid.o -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf


#- On salish (UBC)
# LIBNETCDF=/usr
# mpif90 -c grid.f90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
# mpif90 -c map.F90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
# mpif90 -c get_weight_nemo.F90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
# mpif90 -o get_weight_nemo get_weight_nemo.o map.o grid.o -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf -lnetcdff


#gfortran -c grid.f90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
#gfortran -c map.F90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
#gfortran -c get_weight_nemo.F90 -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf
#gfortran get_weight_nemo.o map.o grid.o -I${LIBNETCDF}/include -L${LIBNETCDF}/lib -lnetcdf



