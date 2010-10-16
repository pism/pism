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
#TWENTYKMGRID="-Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 41 -Mbz 11 -z_spacing equal -zb_spacing equal"

GRID=$TENKMGRID
CS=10

# cat prefix and exec together
PISMR="${PISM_PREFIX}pismr -ocean_kill -e 3"

# coupler settings: Fausto 2m air temp parameterization, but default PDD
#   (w/o Greve 2005 settings of PDD parameters)
COUPLER="-atmosphere searise_greenland -surface pdd"


# config file
CDLCONFIG=base_config.cdl
PCONFIG=base_config.nc
echo
echo "$SCRIPTNAME  generating $PCONFIG file from $CDLCONFIG"
cmd="rm -rf $PCONFIG"
$PISM_DO $cmd
cmd="ncgen -o $PCONFIG $CDLCONFIG"
$PISM_DO $cmd


# bootstrap and do very short smoothing run
SMOOTHRUNLENGTH=0.1
PRE0NAME=g${CS}km_pre${SMOOTHRUNLENGTH}.nc
echo
echo "$SCRIPTNAME  bootstrapping plus very short smoothing run (for ${SMOOTHRUNLENGTH}a)"
cmd="$PISM_MPIDO $NN $PISMR -boot_file $PISM_DATANAME $GRID \
  $COUPLER -y ${SMOOTHRUNLENGTH} -o $PRE0NAME"
$PISM_DO $cmd


# climate in one year, a year in the middle of the period where
#   Ettema et al precip and smb and melt and runoff are available
CLIMSTARTTIME=1990
CLIMENDTIME=1991
#DT=0.0833333333 # monthly = (1/12) of year
DT=1.0
CLIMATE0=g${CS}km_clim_one_year.nc
echo
echo "$SCRIPTNAME  running pclimate with monthly timesteps to generate $CLIMATE0"
echo "$SCRIPTNAME    to get ice surface climate in one year period"
cmd="$PISM_MPIDO $NN ${PISM_PREFIX}pclimate -i $PRE0NAME $COUPLER \
  -config_override $PCONFIG -ys $CLIMSTARTTIME -ye $CLIMENDTIME -dt $DT -o $CLIMATE0"
$PISM_DO $cmd

CLIMATE=mask$CLIMATE0
echo
echo "$SCRIPTNAME  removing some fields and adding in thk (in prep for masking), using NCO;"
echo "$SCRIPTNAME    generating $CLIMATE"
cmd="rm -rf $CLIMATE"
$PISM_DO $cmd
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
echo "$SCRIPTNAME  masking $CLIMATE using climmask.py to remove ice-free areas"
cmd="./climmask.py -v acab,smelt,srunoff $CLIMATE"
$PISM_DO $cmd
echo

echo "$SCRIPTNAME  file $CLIMATE done"

