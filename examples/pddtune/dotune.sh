#!/bin/bash

# Copyright (C) 2010 Ed Bueler

# see README for role of this script

FIXME:  NOT EVEN RUN ONCE

SCRIPTNAME="#(dotune.sh)"

set -e  # exit on error

./preprocess.sh

./boot.sh

for THRESHOLD in 268 270 273.15 275
do
  for DDFSNOW in 0.001 0.003 0.009
  do
    for REFREEZE in 0.2 0.6 1.0
    do
      CONFIG=case_${THRESHOLD}_${SNOW}_${REFREEZE}.nc
      echo "case $CONFIG"
      ncks -O base_config.nc $CONFIG
      
      ncatted pdd_positive_threshold_temp,$THRESHOLD $CONFIG
      ncatted pdd_factor_snow,$SNOW $CONFIG
      ncatted pdd_refreeze,$REFREEZE $CONFIG
      ./runcase.sh $CONFIG start.nc FOO_out.nc
      
      ./climmask.py FOO_out.nc
      
      ./foo.py
    done
  done
done



