#!/bin/bash

# Copyright (C) 2009-2015, 2017, 2023, 2024  PISM authors
##################################################################################
# Complete spinup of Antarctic ice sheet model by regriding to a finer resolution (~15 km).
##################################################################################
#
# Set the environment variable RUN_LENGTH to adjust the run length for testing.

SCRIPTNAME="#(antspin-regridtofine.sh)"

log() {
echo "${SCRIPTNAME} ${*}"
}

set -e  # exit on error

USE_PROJ=${USE_PROJ:-false}

# naming files, directories, executables
RESDIR=
BOOTDIR=
PISM_EXEC=pism
PISM_MPIDO="mpiexec -n "

# input data:
export PISM_INDATANAME=${BOOTDIR}pism_Antarctica_5km.nc

source set-physics.sh

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "antspin-regridtofine.sh 8" then NN = 8
  NN="$1"
fi
log "             NN = $NN"
set -e  # exit on error

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  log "        PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO=""
fi
DO=$PISM_DO

log "            PISM = $PISM_EXEC"
log "        FULLPHYS = $FULLPHYS"
log "         PIKPHYS = $PIKPHYS"
log "PIKPHYS_COUPLING = $PIKPHYS_COUPLING"


# #######################################
## interpolate to a finer grid, reset year to 0 and run for 2000 model years:
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
RUN_LENGTH=${RUN_LENGTH:-2000}

EXTRANAME=${RESDIR}extra_${stage}.nc
exvars="thk,usurf,velbase_mag,velbar_mag,mask,diffusivity,tauc,bmelt,tillwat,temppabase,hardav,ice_area_specific_volume,cell_grounded_fraction"
expackage="-extra_times 0:10:$RUN_LENGTH -extra_vars $exvars"

# Note: switching to a finer grid may produce very high SIA diffusivities
echo
log "continue but regrid to $GRIDNAME and run for 2000 a"
cmd="$PISM_MPIDO $NN $PISM_EXEC
  -grid.recompute_longitude_and_latitude ${USE_PROJ}
  -skip
  -skip_max $SKIP
  -i $PISM_INDATANAME
  -bootstrap $GRID
  -regrid_file $INNAME
  -regrid_vars litho_temp,thk,enthalpy,tillwat,basal_melt_rate_grounded
   $SIA_ENHANCEMENT
   $PIKPHYS_COUPLING
   $PIKPHYS
   $FULLPHYS
  -stress_balance.sia.max_diffusivity 1e5
  -ys 0
  -y $RUN_LENGTH
  -ts_file $TSNAME
  -ts_times 1
  -extra_file $EXTRANAME
   $expackage
  -o $RESNAME
  -o_size big"

$DO $cmd

# one can regrid to 10km, 6.7km, 5km and so on, if the number of model years (RUN_LENGTH)
# is appropriately shortened, and sufficient memory is available

echo
log "fine grid part of spinup done"
