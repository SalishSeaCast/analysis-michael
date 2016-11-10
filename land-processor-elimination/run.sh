#!/bin/bash

MESHMASK=~/MEOPAR/NEMO-forcing/grid/mesh_mask_downbyone2.nc 
TOOL=~/MEOPAR/NEMO-3.6-r6770/NEMOGCM/TOOLS/MPP_PREP/mpp_optimiz_zoom_nc.exe

# The mpp_optimiz_zoom_nc tool expects the bathy file to have
# variable "Bathy_level", but the mesh mask produced by NEMO calls
# it "mbathy". So we rename the variable for use with the tool.
if [ -f $MESHMASK ]; then
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

# Run script to convert processor.layout into tables
python readpl.py

