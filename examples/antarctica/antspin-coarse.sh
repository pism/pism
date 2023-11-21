#!/bin/bash

# Copyright (C) 2009-2018, 2023  PISM authors
##################################################################################
# Coarse grid spinup of Antarctic ice sheet model using data from Anne Le Brocq
# (README.md).  Uses PIK physics and enthalpy model
# and modified configuration parameters with constant climate.  Uses constant
# and precip and a parameterization for artm as in Martin et al (2011).
# WARNING: at finer resolutions (e.g. under 15 km), output is large!
##################################################################################
#
# Set environment variables NOMASS_RUN_LENGTH and CONSTANT_CLIMATE_RUN_LENGTH to adjust run
# lengths for testing.

SCRIPTNAME="#(antspin-coarse.sh)"

log() {
echo "${SCRIPTNAME} ${*}"
}

log "run preprocess.sh before this..."

set -e  # exit on error

log "  Coarse grid, constant-climate spinup script using ALBMAP v1 data"
log "     and -stress_balance ssa+sia and -pik"
log "  Run as './${SCRIPTNAME} NN' for NN procs and ~30km grid"

# naming files, directories, executables
RESDIR=
BOOTDIR=
PISM_EXEC=${PISM_EXEC:=pismr}
PISM_MPIDO="mpiexec -n "

# input data:
export PISM_INDATANAME=${BOOTDIR}pism_Antarctica_5km.nc

source set-physics.sh

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "antspin-coarse.sh 8" then NN = 8
  NN="$1"
fi
log "             NN = $NN"

# check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  log "        PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO=""
fi
DO=$PISM_DO

# these coarse grid defaults are for development/regression runs, not
# "production" or science
GRID=$THIRTYKMGRID
SKIP=$SKIPTHIRTYKM
GRIDNAME=30km

log "            PISM = $PISM_EXEC"
log "        FULLPHYS = $FULLPHYS"
log "         PIKPHYS = $PIKPHYS"
log "PIKPHYS_COUPLING = $PIKPHYS_COUPLING"


# #######################################
# bootstrap and SHORT smoothing run to 100 years
# #######################################
stage=earlyone
INNAME=$PISM_INDATANAME
RESNAMEONE=${RESDIR}${stage}_${GRIDNAME}.nc
RUN_LENGTH=1
echo
log "bootstrapping on $GRIDNAME grid plus SIA run for $RUN_LENGTH a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -skip -skip_max $SKIP -i ${INNAME} -bootstrap $GRID \
	$SIA_ENHANCEMENT $PIKPHYS_COUPLING -front_retreat_file ${INNAME} \
	-y $RUN_LENGTH -o $RESNAMEONE"
$DO $cmd
#exit # <-- uncomment to stop here

stage=smoothing
RESNAME=${RESDIR}${stage}_${GRIDNAME}.nc
RUN_LENGTH=100
echo
log "short SIA run for $RUN_LENGTH a"
cmd="$PISM_MPIDO $NN $PISM_EXEC -skip -skip_max $SKIP -i $RESNAMEONE \
	$SIA_ENHANCEMENT $PIKPHYS_COUPLING -front_retreat_file ${INNAME} \
	-y $RUN_LENGTH -o $RESNAME"
$DO $cmd

# #######################################
# SIA run with fixed geometry on a coarse grid for 200ka
# #######################################
stage=nomass
INNAME=$RESNAME
RESNAME=${RESDIR}${stage}_${GRIDNAME}.nc
TSNAME=${RESDIR}ts_${stage}_${GRIDNAME}.nc
NOMASS_RUN_LENGTH=${NOMASS_RUN_LENGTH:-200000}
EXTRANAME=${RESDIR}extra_${stage}_${GRIDNAME}.nc
expackage="-extra_times 0:1000:${NOMASS_RUN_LENGTH} -extra_vars bmelt,tillwat,velsurf_mag,temppabase,diffusivity,hardav"
echo
log "SIA only run with fixed geometry for ${NOMASS_RUN_LENGTH} years"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i $INNAME $PIKPHYS_COUPLING  \
    $SIA_ENHANCEMENT -no_mass \
    -ys 0 -y ${NOMASS_RUN_LENGTH} \
    -extra_file $EXTRANAME $expackage \
    -o $RESNAME"
$DO $cmd
#exit # <-- uncomment to stop here


# #######################################
# SSA+SIA run "into steady state" (100000 years) with constant climate forcing
# #######################################
stage=run
INNAME=$RESNAME
RESNAME=${RESDIR}${stage}_${GRIDNAME}.nc
TSNAME=${RESDIR}ts_${stage}_${GRIDNAME}.nc
CONSTANT_CLIMATE_RUN_LENGTH=${CONSTANT_CLIMATE_RUN_LENGTH:-100000}
EXTRANAME=${RESDIR}extra_${stage}_${GRIDNAME}.nc
exvars="thk,usurf,velbase_mag,velbar_mag,mask,diffusivity,tauc,bmelt,tillwat,temppabase,hardav,cell_grounded_fraction,ice_area_specific_volume,amount_fluxes,basal_mass_flux_grounded,basal_mass_flux_floating"
expackage="-extra_times 0:1000:${CONSTANT_CLIMATE_RUN_LENGTH} -extra_vars $exvars"

echo
log "SSA+SIA run \"into steady state\" with constant climate forcing for ${CONSTANT_CLIMATE_RUN_LENGTH} years"
cmd="$PISM_MPIDO $NN $PISM_EXEC -bootstrap ${GRID} -skip -skip_max $SKIP -i $INNAME \
    $SIA_ENHANCEMENT $PIKPHYS_COUPLING $PIKPHYS $FULLPHYS \
    -ys 0 -y ${CONSTANT_CLIMATE_RUN_LENGTH} \
    -ts_file $TSNAME -ts_times 0:1:${CONSTANT_CLIMATE_RUN_LENGTH} \
    -extra_file $EXTRANAME $expackage \
    -o $RESNAME -o_size big"
$DO $cmd

echo
log "coarse-grid part of spinup done"
