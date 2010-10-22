#!/bin/bash

# Copyright (C) 2010 Ed Bueler

# Uses a SeaRISE-Greenland data set to illustrate the use of regional
# climate model (RCM) output to find PDD parameters.  The goal is for the
# surface mass balance from PISM's PDD model to closely fit the corresponding
# RACMO/GR regional climate model output, from Ettema et al.

# This is the top-level script.  See also README.

# Suggested way to run:  $ ./dotune.sh >& out.dotune &

# This script uses NCO (http://nco.sourceforge.net/).

set -e  # exit on error

DATANAME=Greenland_5km_v1.1.nc
PISMDATA=pism_$DATANAME

DIFFSFILE=diffs.txt

./preprocess.sh # generates pism_Greenland_5km_v1.1.nc and base_config.nc

./boot.sh  # creates start.nc, which contains 'thk' used in masking in objective.py

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

      CLIMATE=clim_$NAMEROOT.nc
      ./runcase.sh $CONFIG start.nc $CLIMATE
      rm -rf $CONFIG  # don't need this file any more BECAUSE pism_overrides are
                      #   carried forward into $CLIMATE
      
      echo
      echo "  computing objective function by comparing 'acab' in $CLIMATE"
      echo "    to 'smb' in $PISMDATA and putting objective value in $DIFFSFILE"
      ./objective.py -v acab,smb -H start.nc $CLIMATE $PISMDATA $DIFFSFILE
      echo
    done
  done
done

rm -f tempthk.nc


