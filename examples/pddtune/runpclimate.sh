#!/bin/bash

# Copyright (C) 2010 Ed Bueler

# uses the PISM SeaRISE Greenland example to illustrate the use of regional
# climate model (RCM) output to find PDD parameters which produce closer fit of
# surface mass balance from PISM's PDD model to the RCM model output

# depends on all tools in examples/searise-greenland

# recommended way to run with N processors is " ./FOO.sh N >& out.FOO & "


SCRIPTNAME="#(runpclimate.sh)"


set -e  # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "psearise.sh 8" then NN = 8
  NN="$1"
fi

echo "$SCRIPTNAME              NN = $NN"

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


PISM_DATAVERSION=1.1
PISM_DATANAME=pism_Greenland_5km_v$PISM_DATAVERSION.nc

for INPUT in $PISM_DATANAME; do
  if [ -e "$INPUT" ] ; then  # check if file exist
    echo "$SCRIPTNAME           input   $INPUT (found)"
  else
    echo "$SCRIPTNAME           input   $INPUT (MISSING!!)"
    echo
    echo "     !!!!   RUN  preprocess.sh  TO GENERATE  $INPUT   !!!!"
    echo
  fi
done


# grids
TENKMGRID="-Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 41 -Mbz 11 -z_spacing equal -zb_spacing equal"
TWENTYKMGRID="-Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 41 -Mbz 11 -z_spacing equal -zb_spacing equal"

GRID=$TWENTYKMGRID
CS=20

# cat prefix and exec together
PISM="${PISM_PREFIX}pismr -ocean_kill -e 3"

# coupler settings
COUPLER="-atmosphere searise_greenland -surface pdd -pdd_fausto"

# bootstrap and do smoothing run for 1 year
SMOOTHRUNLENGTH=1
PRE0NAME=g${CS}km_pre${SMOOTHRUNLENGTH}.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run (for ${SMOOTHRUNLENGTH}a)"
echo
cmd="$PISM_MPIDO $NN $PISM -boot_file $PISM_DATANAME $GRID \
  $COUPLER -y ${SMOOTHRUNLENGTH} -o $PRE0NAME"
$PISM_DO $cmd


# quick look at climate in a particular period
CLIMSTARTTIME=1958
#CLIMENDTIME=2007
CLIMENDTIME=1977
CLIMATE0=g${CS}km_clim_${CLIMSTARTTIME}-${CLIMENDTIME}.nc
echo
echo "$SCRIPTNAME  running pclimate with monthly timesteps to generate $CLIMATE0"
echo "$SCRIPTNAME    to get ice surface climate in period $CLIMSTARTTIME-$CLIMENDTIME"
echo
cmd="$PISM_MPIDO $NN ${PISM_PREFIX}pclimate -i $PRE0NAME $COUPLER \
  -ys $CLIMSTARTTIME -ye $CLIMENDTIME -dt 0.0833333333 -o $CLIMATE0"
$PISM_DO $cmd


CLIMATE=will_break_it.nc
echo
echo "$SCRIPTNAME  removing some fields and adding in thk (in prep for masking), using NCO;"
echo "$SCRIPTNAME    generating $CLIMATE"
echo
cmd="ncks -O $CLIMATE0 $CLIMATE"
$PISM_DO $cmd
cmd="ncks -O -x -v shelfbasetemp,shelfbasemassflux $CLIMATE $CLIMATE"
$PISM_DO $cmd
cmd="ncecat -O -v thk $PRE0NAME tempthk.nc"
$PISM_DO $cmd
cmd="ncwa -O -a t -v thk $PRE0NAME tempthk.nc"
$PISM_DO $cmd
cmd="ncks -A -v thk tempthk.nc $CLIMATE"
$PISM_DO $cmd
cmd="rm -f tempthk.nc"
$PISM_DO $cmd

echo
echo "$SCRIPTNAME  masking using python script to remove ice-free areas"
echo
cmd="./climmask.py -v acab,artm,smelt,srunoff $CLIMATE"
$PISM_DO $cmd
echo

