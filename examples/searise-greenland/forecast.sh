#!/bin/bash

# Copyright (C) 2009-2010 Andy Aschwanden and Ed Bueler

# PISM SeaRISE Greenland "forecast" script which merely does the SeaRISE
# "control run".  Before using this script, run preprocess.sh to download and
# adjust metadata, then run spinup.sh to do spinup.  The final state used
# here is assumed to be called "g5km_0.nc" (or set env var PISM_SPUNUP).  See
# the PISM User's Manual.
#
# Recommended way to run with N processors is "./forecast.sh N >& out.fore &"

if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(forecast.sh)"
fi

echo
echo "# =================================================================================="
echo "# PISM SeaRISE Greenland: forecast"
echo "# =================================================================================="
echo

set -e  # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "psearise-forecast.sh 8" then NN = 8
  NN="$1"
fi

echo   "$SCRIPTNAME              NN = $NN"

if [ -n "${PISM_SPUNUP:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_SPUNUP = $PISM_SPUNUP  (already set)"
else
  PISM_SPUNUP=g5km_0.nc
  echo "$SCRIPTNAME     PISM_SPUNUP = $PISM_SPUNUP"
fi

# set PISM_MPIDO if using different MPI execution command, for example:
#  $ export PISM_MPIDO="aprun -n "
if [ -n "${PISM_MPIDO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO  (already set)"
else
  PISM_MPIDO="mpiexec -n "
  echo "$SCRIPTNAME      PISM_MPIDO = $PISM_MPIDO"
fi

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO="" 
fi

# prefix to pismr (e.g. PISM_PREFIX=/home/username/pism-dev/bin/)
if [ -n "${PISM_PREFIX:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX  (already set)"
else
  PISM_PREFIX=  # just a guess
  echo "$SCRIPTNAME     PISM_PREFIX = $PISM_PREFIX"
fi

# set PISM_EXEC if using different executables or options, for example:
#  $ export PISM_EXEC="pismr -cold"
if [ -n "${PISM_EXEC:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC  (already set)"
else
  PISM_EXEC="pismr"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi

# cat prefix and exec together
PISM="${PISM_PREFIX}${PISM_EXEC} -ocean_kill -eta"

COUPLER="-atmosphere searise_greenland -surface pdd -pdd_fausto"

# use "control run" parameters from Bueler et al. submitted
PARAMS="-e 3 -pseudo_plastic_q 0.25 -plastic_pwfrac 0.98"

FULLPHYS="-ssa_sliding -thk_eff ${PARAMS}"

echo "$SCRIPTNAME      executable = '$PISM'"
echo "$SCRIPTNAME    full physics = '$FULLPHYS'"
echo "$SCRIPTNAME         coupler = '$COUPLER'"


pismopts="$PISM_MPIDO $NN $PISM -skip 50 $FULLPHYS -bed_def lc"

expackage="-extra_vars usurf,topg,thk,bmelt,bwat,dHdt,mask,uvelsurf,vvelsurf,wvelsurf,uvelbase,vvelbase,wvelbase,tempsurf,tempbase"

ENDTIME=500
EXTSTIMES=0:5:${ENDTIME}

OUTNAME=control_y${ENDTIME}.nc
EXNAME=ex_y${ENDTIME}.nc
TSNAME=ts_y${ENDTIME}.nc
echo
echo "$SCRIPTNAME  5km grid: control run from 0 to $ENDTIME years w save every 5 years:"
echo
cmd="$pismopts -i $PISM_SPUNUP $COUPLER -ye $ENDTIME \
  -extra_file $EXNAME -extra_times $EXTSTIMES $expackage \
  -ts_file $TSNAME -ts_times $EXTSTIMES -o $OUTNAME"
$PISM_DO $cmd
echo
echo "$SCRIPTNAME  steady-climate control run done"


