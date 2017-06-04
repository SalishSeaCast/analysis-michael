#/bin/bash

set -e
#set -x

TOOLPATH=$HOME/MEOPAR/NestingTools/NEMOGCM/TOOLS/NESTING
WEIGHTSPATH=$HOME/MEOPAR/analysis-michael/weights
EMBEDPATH=$HOME/MEOPAR/analysis-michael/agrif

NEMOFORCING=$HOME/MEOPAR/NEMO-forcing
NEMOFORCINGMD=$HOME/MEOPAR/NEMO-forcing-MD

if [ $# -lt 2 ]; then
  echo "Usage: $0 <nesting_namelist> <action>"
  echo ""
  echo "where <nesting_namelist> is selected from \$HOME/MEOPAR/analysis-michael/agrif"
  echo "Example:"
  echo "    mkdir -p output && cd output"
  echo "    \$HOME/MEOPAR/analysis-michael/agrif/run.sh baynessound1 all"
  echo ""
  echo "Available actions: coords bathy weights restart restart_trc rivers"
  echo ""
  echo "Lastly ensure that the filename for <nesting_namelist> is not"
  echo "called 'namelist' as to avoid conflicting with the weights tool"
  exit 1;
fi

NMLIST=$EMBEDPATH/$1
ACTION=$2

case "$ACTION" in
#    all) # Make everything -- this is not so useful yet
#        for action in coords bathy weights rivers restart restart_trc; do
#         for action in coords bathy weights rivers; do
#            bash $0 $1 $action
#        done
#        ;;

    coords)
        # Build the coordinates file
        $TOOLPATH/agrif_create_coordinates.exe $NMLIST

        # Write AGRIF_FixdGrids.in
        python $EMBEDPATH/make_agfi.py $NMLIST

        # Store a copy of the main namelist here
        cp $NMLIST .
        ;;

    bathy)
        $TOOLPATH/agrif_create_bathy.exe $NMLIST

        # Enforce minimum depth on the child grid and fix nav_lon
        python $EMBEDPATH/fix_bathy.py 1_bathymetry_201702.nc
        ;;

    weights)
        # We need a surface forcing file called 'ops.nc'
        cp -p $NEMOFORCING/atmospheric/ops_y2015m06d13.nc ops.nc

        # The weights tool expects input to be called bathy_meter.nc
        ln -fs 1_bathymetry_201702.nc bathy_meter.nc

        # Run the tool
        cp -p $WEIGHTSPATH/namelist .
        $WEIGHTSPATH/get_weight_nemo

        # Convert to netcdf4 and improve the metadata
        python $EMBEDPATH/fix_weights.py
        
        # Tidying
        rm -f namelist met_gem_weight.nc bathy_meter.nc ops.nc

        # Set up 1_atmospheric
        mkdir -p 1_atmospheric
        mv 1_weights-gem2.5-ops.nc 1_atmospheric/1_weights-gem2.5-ops_201702.nc
        cd 1_atmospheric
        ln -fs ../../NEMO-forcing-MD/atmospheric/no_snow.nc
        for x in ../../NEMO-forcing-MD/atmospheric/ops*nc; do
          ln -fs $x
        done
        ;;

    rivers)
        # Symlink to where the river files are
        ln -fs $NEMOFORCINGMD/rivers rivers
        ln -fs $NEMOFORCING/rivers/bio_climatology bio_climatology

        # Create the data
        $TOOLPATH/agrif_create_data.exe $NMLIST

        # Tidying
        rm -f rivers
        rm -f bio_climatology

        # Set up 1_rivers
        mkdir -p 1_rivers
        mv 1_river_ConsTemp_month.nc 1_rivers_month_201702.nc 1_R201702DFraCElse* 1_rivers

        # Set up 1_bio_climatology
        mkdir -p 1_bio_climatology
        mv 1_riverBio* 1_bio_climatology
        ;;

    restart)
        $TOOLPATH/agrif_create_restart.exe $NMLIST
        ;;

    restart_trc)
        $TOOLPATH/agrif_create_restart_trc.exe $NMLIST
        ;;

    *)
        echo "Unknown Action"
        exit 1
        ;;
esac

