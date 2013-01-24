#!/bin/bash

# Copyright (C) 2013 The PISM Authors

# This is just a helper script to make running EISMINT II experiments easier.
# It adds suggested diagnostics which help compare to the published experiments.

SCRIPTNAME="#(runexp.sh)"

set -e  # exit on error

PISM_MPIDO="mpiexec -n "

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "runexp.sh 8" then NN = 8
  NN="$1"
fi

EXP=A  # default experiment
if [ $# -gt 1 ] ; then  # if user says "runexp.sh 8 F" then NN = 8 and EXP = F
  EXP="$2"
fi

MM=61  # default resolution
if [ $# -gt 2 ] ; then  # if user says "runexp.sh 8 F 121" then NN = 8 and EXP = F and MM = 121
  MM="$3"
fi

# choose -i or start from grid
if [ $# -gt 3 ] ; then  # if user says "runexp.sh 8 F 121 foo.nc" then NN = 8 and
                        #    EXP = F and MM is ignored and '-i foo.nc' is used
  GRIDORINPUT="-i $4"
else
  GRIDORINPUT="-Mx $MM -My $MM -Mz $MM"
fi

if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO="" 
fi

LENGTH=200000  # eismint default
SKIP=10        # adjust downward for low res, upward for high res
echo "$SCRIPTNAME  run into steady state with constant climate forcing for $RUNTIME a"
cmd="$PISM_MPIDO $NN pisms -eisII $EXP $GRIDORINPUT -ys 0 -y $LENGTH \
 -skip -skip_max $SKIP -o eisII$EXP$MM.nc -extra_file ex_eisII$EXP$MM.nc \
 -extra_vars thk,temppabase,csurf,cbar,cflx,diffusivity,bmelt,taud_mag \
 -extra_times 1000:1000:$LENGTH -ts_file ts_eisII$EXP$MM.nc \
 -ts_times 0:100:$LENGTH -ts_vars ivol,iarea,iareatemp,iareacold"
$PISM_DO $cmd

