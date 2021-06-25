#!/bin/bash

# Copyright (C) 2013--2016, 2021 The PISM Authors

# This is just a helper script to make running EISMINT II experiments easier.
# It adds suggested diagnostics which help compare to the published experiments.

SCRIPTNAME="#(runexp.sh)"

set -e  # exit on error

PISM_MPIDO="mpiexec -n "

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "runexp.sh 8" then NN=8
  NN="$1"
fi

EXP=A  # default experiment
if [ $# -gt 1 ] ; then  # if user says "runexp.sh 8 F" then NN=8, EXP=F
  EXP="$2"
fi

MHOR=61  # default resolution
if [ $# -gt 2 ] ; then  # if user says "runexp.sh 8 F 121" then NN=8, EXP=F, Mx=My=121
  MHOR="$3"
fi

MVER=$MHOR  # default vertical resolution
if [ $# -gt 3 ] ; then  # if user says "runexp.sh 8 F 121 201" then NN=8, EXP=F, Mx=My=121, Mz=201
  MVER="$4"
fi

DUR=2e5  # default duration of 200000 a
if [ $# -gt 4 ] ; then  # if user says "runexp.sh 8 F 121 201 1e4" then NN=8, EXP=F, Mx=My=121, Mz=201, -y 1e4
  DUR="$5"
fi

# choose -i or start from grid
if [ $# -gt 5 ] ; then  # if user says "runexp.sh 8 F X X foo.nc" then NN=8, EXP=F, next two are ignored, -y 1e4, and '-i foo.nc' is used
  GRIDORINPUT="-i $6"
else
  if [ "$EXP" = "F" ] ; then
    LZ="-Lz 6000"
  else
    LZ="-Lz 5000"
  fi
  GRIDORINPUT="-Mx $MHOR -My $MHOR -Mz $MVER $LZ"
fi

if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO="" 
fi

SKIP=5         # adjust upward for high res
ROOT=eisII$EXP$MHOR
echo "$SCRIPTNAME  run into steady state with constant climate forcing for $RUNTIME a"
cmd="$PISM_MPIDO $NN pismr -eisII $EXP $GRIDORINPUT -ys 0 -y $DUR \
 -skip -skip_max $SKIP -o $ROOT.nc -extra_file ex_$ROOT.nc \
 -extra_vars thk,temppabase,velsurf_mag,velbar_mag,flux_mag,diffusivity,bmelt,taud_mag \
 -extra_times 1000:1000:$DUR -ts_file ts_$ROOT.nc \
 -ts_times 0:100:$DUR"
$PISM_DO $cmd

