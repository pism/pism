#!/bin/bash

# Copyright (C) 2010 Ed Bueler

# Uses a SeaRISE-Greenland data set to illustrate the use of regional
# climate model (RCM) output to find PDD parameters.  The goal is for the
# surface mass balance from PISM's PDD model to closely fit the corresponding
# RACMO/GR regional climate model output, from Ettema et al.

# This is the top-level script.  See also README.

# This script uses NCO (http://nco.sourceforge.net/).

set -e  # exit on error

DATANAME=Greenland_5km_v1.1.nc
PISMDATA=pism_$DATANAME

./preprocess.sh # generates pism_Greenland_5km_v1.1.nc and base_config.nc

./boot.sh  # creates start.nc

# now put thickness map in tempthk.nc, w/o degenerate t axis:
ncecat -O -v thk start.nc tempthk.nc
ncwa -O -a t -v thk start.nc tempthk.nc

for THRESHOLD in 268 270 273
do
  for DDFSNOW in 0.001 0.003 0.009
  do
    for REFREEZE in 0.2 0.6 1.0
    do
      NAMEROOT=${THRESHOLD}_${DDFSNOW}_${REFREEZE}
      CONFIG=config_${THRESHOLD}_${DDFSNOW}_${REFREEZE}.nc
      echo "case ${NAMEROOT}:"
      echo "  creating -config_override file $CONFIG ..."
      rm -rf $CONFIG
      ncks -O base_config.nc $CONFIG      
      ncatted -O -a pdd_positive_threshold_temp,pism_overrides,m,d,$THRESHOLD $CONFIG
      ncatted -O -a pdd_factor_snow,pism_overrides,m,d,$DDFSNOW $CONFIG
      ncatted -O -a pdd_refreeze,pism_overrides,m,d,$REFREEZE $CONFIG

      CLIMATE0=out_$NAMEROOT.nc
      ./runcase.sh $CONFIG start.nc $CLIMATE0

      CLIMATE=clim_$NAMEROOT.nc
      echo "  removing some fields from $CLIMATE0 and adding in thk in prep for masking;"
      echo "  generating $CLIMATE ..."
      rm -rf $CLIMATE
      ncks -O $CLIMATE0 $CLIMATE
      ncks -O -x -v shelfbasetemp,shelfbasemassflux $CLIMATE $CLIMATE
      ncks -A -v thk tempthk.nc $CLIMATE  # put the thk variable in $CLIMATE
      rm -rf $CLIMATE0  # don't need this file any more

      echo "  masking $CLIMATE using climmask.py to remove ice-free areas from consideration"
      ./climmask.py -v acab,smelt,srunoff,saccum $CLIMATE
      
      echo "  computing objective function by comparing 'acab' in $CLIMATE"
      echo "    to 'smb' in $PISMDATA"
      #FIXME:  ./objective.py -v acab -s smb $CLIMATE $PISMDATA
      echo
    done
  done
done

rm -f tempthk.nc


