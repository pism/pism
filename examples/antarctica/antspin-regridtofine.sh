#!/bin/bash

# Copyright (C) 2009-2015, 2017  PISM authors
##################################################################################
# Complete spinup of Antarctic ice sheet model by regriding to a finer resolution (15 km).
##################################################################################

SCRIPTNAME="#(antspin-regridtofine.sh)"

set -e  # exit on error

# naming files, directories, executables
RESDIR=
BOOTDIR=
PISM_EXEC=pismr
PISM_MPIDO="mpiexec -n "

# input data:
export PISM_INDATANAME=${BOOTDIR}pism_Antarctica_5km.nc

source set-physics.sh

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "antspinCC.sh 8" then NN = 8
  NN="$1"
fi
echo "$SCRIPTNAME              NN = $NN"
set -e  # exit on error

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO=""
fi
DO=$PISM_DO

echo "$SCRIPTNAME             PISM = $PISM_EXEC"
echo "$SCRIPTNAME         FULLPHYS = $FULLPHYS"
echo "$SCRIPTNAME          PIKPHYS = $PIKPHYS"
echo "$SCRIPTNAME PIKPHYS_COUPLING = $PIKPHYS_COUPLING"


# #######################################
## do a regridding to fine grid and reset year to 0 and run for 2000 model years:
# #######################################
GRID=$FIFTEENKMGRID
SKIP=$SKIPFIFTEENKM
GRIDNAME=15km

COARSENAME=run_30km.nc
RESDIR=

stage=run_regrid_${GRIDNAME}
INNAME=$COARSENAME
RESNAME=${RESDIR}${stage}.nc
TSNAME=${RESDIR}ts_${stage}.nc
RUNTIME=2000

EXTRANAME=${RESDIR}extra_${stage}.nc
exvars="thk,usurf,velbase_mag,velbar_mag,mask,diffusivity,tauc,bmelt,tillwat,temppabase,hardav,Href,gl_mask"
expackage="-extra_times 0:10:$RUNTIME -extra_vars $exvars"

echo
echo "$SCRIPTNAME  continue but regrid to $GRIDNAME and run for 2000 a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -skip -skip_max $SKIP \
    -i $PISM_INDATANAME -bootstrap $GRID \
    -regrid_file $INNAME -regrid_vars litho_temp,thk,enthalpy,tillwat,basal_melt_rate_grounded \
    $SIA_ENHANCEMENT $PIKPHYS_COUPLING $PIKPHYS $FULLPHYS \
    -ys 0 -y $RUNTIME \
    -ts_file $TSNAME -ts_times 0:1:$RUNTIME \
    -extra_file $EXTRANAME $expackage \
    -o $RESNAME -o_size big"

$DO $cmd

# one can regrid to 10km, 6.7km, 5km and so on, if the number of model years
# (RUNTIME) is appropriately shortened, and sufficient memory is available

echo
echo "$SCRIPTNAME  fine-grid part of spinup done"

