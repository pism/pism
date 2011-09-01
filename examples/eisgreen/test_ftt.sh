#!/bin/bash

NN=8  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "psearise.sh 8" then NN = 8
  NN="$1"
fi

# set MPIDO if using different MPI execution command, for example:
#  $ export PISM_MPIDO="aprun -n "
if [ -n "${PISM_MPIDO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO  (already set)"
else
  PISM_MPIDO="mpiexec -n "
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO"
fi

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var DO is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO="" 
fi

# prefix to pism (not to executables)
if [ -n "${PISM_PREFIX:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX  (already set)"
else
  PISM_PREFIX=""    # just a guess
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX"
fi

# set PISM_EXEC if using different executables, for example:
#  $ export PISM_EXEC="pismr -cold"
if [ -n "${PISM_EXEC:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC  (already set)"
else
  PISM_EXEC="pismr"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi


PISM="${PISM_PREFIX}${PISM_EXEC}"

cmd="$PISM_MPIDO $NN $PISM -ys -500.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
    -surface pdd -pdd_fausto \
    -o no_force.nc -ts_file ts_no_force.nc -ts_times -500:10:0"
#$PISM_DO $cmd

echo

cmd="$PISM_MPIDO $NN $PISM -ys -500.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
  -surface pdd,forcing -pdd_fausto -force_to_thk green20km_y1.nc \
  -o with_force.nc -ts_file ts_with_force.nc -ts_times -500:10:0"
#$PISM_DO $cmd

echo

cmd="$PISM_MPIDO $NN $PISM -ys -500.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
    -surface pdd,forcing -pdd_fausto -force_to_thk green20km_y1.nc -force_to_thk_alpha 0.005 \
    -o weak_force.nc -ts_file ts_weak_force.nc -ts_times -500:10:0"
#$PISM_DO $cmd


echo
echo "Test restartability"

cmd="$PISM_MPIDO $NN $PISM -ys -50.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -surface constant,forcing \
  -force_to_thk green20km_y1.nc -o foo.nc"
$PISM_DO $cmd

echo

cmd="$PISM_MPIDO $NN $PISM -ys -50.0 -ye 25 -skip 5 -i green_ssl2_110ka.nc -surface constant,forcing \
  -force_to_thk green20km_y1.nc -o joe.nc"
$PISM_DO $cmd

echo

cmd="$PISM_MPIDO $NN $PISM -ys -25.0 -ye 0 -skip 5 -i joe.nc -surface constant,forcing \
  -force_to_thk green20km_y1.nc -o bar.nc"
$PISM_DO $cmd

echo

cmd="nccmp.py -v thk foo.nc bar.nc"
$PISM_DO $cmd
