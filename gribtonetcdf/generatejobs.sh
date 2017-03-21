#!/bin/bash

# Generates a list of grib processing jobs

G2N=nowcast.workers.grib_to_netcdf
CFG=$HOME/MEOPAR/SalishSeaNowcast/config/wcvi.yaml

# Blank the list of tasks
>jobs.txt.in

cd /results/forcing/atmospheric/GEM2.5/GRIB
for x in ????????; do
   cd $HOME/MEOPAR/SalishSeaNowcast/nowcast/workers/
   d=${x:0:4}-${x:4:2}-${x:6:2}

   # Add job to list (serial approach)
   # echo "python -m $G2N --debug $CFG nowcast+ --run-date $d" >> jobs.txt.in

   # Add job to list (for use with gnu parallel)
   echo "source activate salishsea-nowcast && python -m $G2N --debug $CFG nowcast+ --run-date $d" >> jobs.txt.in
done

# Skip the first 85 days (Sept 11, 2014 through Dec 3, 2014, because the nowcast worker expects files from Dec 4, 2014 onwards)
tail -n +85 jobs.txt.in > jobs.txt
rm -f jobs.txt.in

