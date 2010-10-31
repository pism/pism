#!/bin/bash

# Copyright (C) 2010 Ed Bueler

# Uses a SeaRISE-Greenland data set to illustrate the use of regional
# climate model (RCM) output to find PDD parameters.  The goal is for the
# surface mass balance from PISM's PDD model to closely fit the corresponding
# RACMO/GR regional climate model output, from Ettema et al.

# This is the top-level script, but see the README.
# Suggested way to run:  $ ./dotune.sh >& out.dotune &

set -e  # exit on error

# do these first:
#  ./preprocess.sh          # generates pism_Greenland_5km_v1.1.nc
#  ./boot.sh >> log.boot    #     "     start.nc

# harder: 5^3 * 11 = 1375 basic cases; linesearch happens on each one so total perhaps 8000
THRESHOLDRANGE="268 269 270 271 273"           # default 273.15
DDFSNOWRANGE="0.001 0.002 0.003 0.005 0.009"   # default 0.003
REFREEZERANGE="0.4 0.5 0.6 0.7 0.8"            # default 0.6
STDDEVLAPSERANGE="-0.5 -0.3 -0.2 -0.1 -0.03 0.0 0.03 0.1 0.2 0.3 0.5" # default 0.0

# easier: 3^3 * 5 = 135 basic cases; linesearch happens on each one so total perhaps 1000
THRESHOLDRANGE="268 270 273"                   # default 273.15
DDFSNOWRANGE="0.001 0.002 0.003 0.005 0.009"   # default 0.003
REFREEZERANGE="0.4 0.6 0.8"                    # default 0.6
STDDEVLAPSERANGE="-0.3 -0.1 0.0 0.1 0.3"       # default 0.0

for THRESHOLD in $THRESHOLDRANGE
do
  for DDFSNOW in $DDFSNOWRANGE
  do
    for REFREEZE in $REFREEZERANGE
    do
      for STDDEVLAPSE in $STDDEVLAPSERANGE
      do
        ./linesearch.py --thresh=$THRESHOLD \
                        --snow=$DDFSNOW \
                        --refreeze=$REFREEZE \
                        --lapse=$STDDEVLAPSE \
                        --tol=0.0001 \
                        --diffsfile=diffs.txt --startfile=start.nc --deletenc
      done
    done
  done
done

