#!/bin/bash

# Copyright (C) 2009-2010 Andy Aschwanden and Ed Bueler
# PISM SeaRISE Greenland worked example

# before using this script:
#    1. run preprocess.sh to download SeaRISE "Present Day Greenland" master
#       dataset and adjust metadata
#    2. run spinup.sh to do spinup; creates g5km_0_rt.nc
#    3. run forecast.sh which takes g5km_0_rt.nc as input and computes SeaRISE
#       future runs into these files:
#           UAF!_G_D3_C?_E?_raw.nc
#           ts_UAF!_G_D3_C?_E?.nc
#    4. run this current script to produce
#           UAF!_G_D3_C?_E?.nc

# this script requires NCO:   http://nco.sourceforge.net/

set -e  # exit on error
set -x  # see commands as they are issued

E=1  # default model number
if [ $# -gt 0 ] ; then  # model number
  E="$1"
fi


# process files to remove fields, combine spatial and scalar series, and fix metadata
for NAME in "UAF${E}_G_D3_C1_E0" "UAF${E}_G_D3_C2_E0" "UAF${E}_G_D3_C3_E0" "UAF${E}_G_D3_C4_E0" "UAF1_G_D3_C1_S1" "UAF1_G_D3_C1_S2""UAF1_G_D3_C1_S3" "UAF1_G_D3_C1_M1" "UAF1_G_D3_C1_M2" "UAF1_G_D3_C1_M3"; do
  echo "(postprocess.sh)  working on deliverable $NAME.nc ..."
  echo "(postprocess.sh)    removing unreported fields and the vertical dimension ..."
  ncks -v artm,snowtemp,snowprecip,surftempoffset,sealevel,z,zb -x ${NAME}_raw_y*.nc \
         -o ${NAME}.nc                                     # creates draft of deliverable file
  echo "(postprocess.sh)    combining annual scalar time series with spatial series file ..."
  cp ts_y*_${NAME}.nc NEWTIME_ts_y*_${NAME}.nc
  ncrename -d t,tseries NEWTIME_ts_y*_${NAME}.nc
  ncrename -v t,tseries NEWTIME_ts_y*_${NAME}.nc
  ncecat -O NEWTIME_ts_y*_${NAME}.nc NEWTIME_ts_y*_${NAME}.nc # convert time to non-record dimension
  ncwa -O -a record NEWTIME_ts_y*_${NAME}.nc NEWTIME_ts_y*_${NAME}.nc # remove just-added record dimension
  ncks -A -v ivol,iareag,iareaf NEWTIME_ts_y*_${NAME}.nc -o ${NAME}.nc # actually combine
  rm NEWTIME_ts_y*_${NAME}.nc
  echo "(postprocess.sh)    fixing metadata and names ..."
  ncrename -v bwat,bwa ${NAME}.nc                          # fix "bwa" name
  ncatted -a pism_ssa_velocities_are_valid,global,d,, ${NAME}.nc  # why does this make it into -extra_file?
  ncatted -a institution,global,c,c,"University of Alaska, Fairbanks" ${NAME}.nc 
  echo "(postprocess.sh)    fixing mask to conform to spec, using postprocess_mask.py ..."
  ./postprocess_mask.py ${NAME}.nc
  echo "(postprocess.sh)    file $NAME.nc done "
done

# process the AR4 results to remove a bit more stuff
for NAME in "UAF${E}_G_D3_C2_E0" "UAF${E}_G_D3_C3_E0" "UAF${E}_G_D3_C3_E0"; do
  echo "(postprocess.sh)    removing extra variables from anomaly forcing ..."
  ncks -O -v delta_snowtemp,delta_snowprecip,snowprecip_inst -x ${NAME}.nc -o ${NAME}.nc
  echo "(postprocess.sh)    file $NAME.nc REALLY done"
done

