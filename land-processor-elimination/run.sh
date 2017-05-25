#!/bin/bash

TOOL=$HOME/MEOPAR/NEMO-3.6-code/NEMOGCM/TOOLS/MPP_PREP/mpp_optimiz_zoom_nc.exe

# Select mesh mask
#MESHMASK=$HOME/MEOPAR/NEMO-forcing/grid/mesh_mask_downbyone2.nc
MESHMASK=$HOME/MEOPAR/NEMO-forcing/grid/mesh_mask201702.nc

# The mpp_optimiz_zoom_nc tool expects the mesh mask file to have
# variable "Bathy_level", but the mesh mask produced by NEMO calls
# it "mbathy". So we rename the variable for use with the tool.
if [ -f $MESHMASK ]; then
  rm -f mesh_mask_renamed.nc
  ncrename -h -v mbathy,Bathy_level $MESHMASK mesh_mask_renamed.nc
else
  echo "Mesh mask not found"
  exit
fi

# Run the mpp_optimiz_zoom_nc tool
if [ -x $TOOL ]; then
  $TOOL
else
  echo "Tool not found"
  exit
fi

# Run script to convert processor.layout into lookup table and preferred decomposition list
python readpl.py

# Tidying
rm -f mesh_mask_renamed.nc processor.layout

