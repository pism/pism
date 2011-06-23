#!/bin/bash

# Copyright (C) 2009-2011 The PISM Authors
#
# before using this script:
#    1. run preprocess.sh to download SeaRISE "Present Day Greenland" master
#       dataset and adjust metadata
#    2. run spinup.sh to do spinup;
#    3. run experiments.sh to do
#       experiment runs into these files:
#           UAF1_G_D3_C?_E?_raw.nc
#           ts_UAF1_G_D3_C?_E?.nc
#    4. run this current script 
#         ./postprocess.sh
#       (or "./postprocess.sh ABC9" if runs are ABC9_...)
#       to produce
#           UAF1_G_D3_C?_E?.nc

# this script requires NCO:   http://nco.sourceforge.net/
# this script depends on postprocess_mask.py, and thus on python and netcdf4-python

set -e  # exit on error
set -x  # see commands as they are issued

MODEL=UAF1  # default initials and model number
if [ $# -gt 0 ] ; then
  MODEL="$1"
fi

INSTITUTION="University of Alaska, Fairbanks"  # will appear as global attribute in final file
if [ $# -gt 1 ] ; then
  INSTITUTION="$2"
fi

# process files to remove unneeded fields, combine spatial and scalar series, and fix metadata
for NAME in "${MODEL}_G_D3_C1_E0" \
            "${MODEL}_G_D3_C2_E0" "${MODEL}_G_D3_C3_E0" "${MODEL}_G_D3_C4_E0" \
            "${MODEL}_G_D3_C1_S1" "${MODEL}_G_D3_C1_S2" "${MODEL}_G_D3_C1_S3" \
            "${MODEL}_G_D3_C1_M1" "${MODEL}_G_D3_C1_M2" "${MODEL}_G_D3_C1_M3"; do

  echo "(postprocess.sh)  working on deliverable $NAME.nc ..."

  echo "(postprocess.sh)    removing unreported fields ..."
  # create draft of deliverable file and remove two early-diagnosis fields:
  ncks -v cbase,csurf,pism_overrides -x ${NAME}_raw_y*.nc -o ${NAME}.nc 
  echo "(postprocess.sh)    combining annual scalar time series with spatial series file ..."
  cp ts_y*_${NAME}.nc NEWTIME_ts_y*_${NAME}.nc
  ncrename -d t,tseries NEWTIME_ts_y*_${NAME}.nc
  ncrename -v t,tseries NEWTIME_ts_y*_${NAME}.nc
  ncecat -O NEWTIME_ts_y*_${NAME}.nc NEWTIME_ts_y*_${NAME}.nc # convert time to non-record dimension
  ncwa -O -a record NEWTIME_ts_y*_${NAME}.nc NEWTIME_ts_y*_${NAME}.nc # remove just-added record dimension
  #FIXME  do we want to preserve time bounds?
  ncks -A -v ivol,iareag,iareaf NEWTIME_ts_y*_${NAME}.nc -o ${NAME}.nc # actually combine
  rm NEWTIME_ts_y*_${NAME}.nc

  echo "(postprocess.sh)    fixing metadata and names ..."
  ncrename -v bwat,bwa ${NAME}.nc                          # fix "bwa" name
  #FIXME  note desired output gline_flx; what to do?
  ncatted -a institution,global,c,c,"${INSTITUTION}" ${NAME}.nc 

  echo "(postprocess.sh)    fixing mask to conform to spec, using postprocess_mask.py ..."
  ./postprocess_mask.py ${NAME}.nc

  echo "(postprocess.sh)    file $NAME.nc done "
done

