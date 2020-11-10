#!/bin/bash

# Copyright (C) 2018 PISM authors
# Based on bootstrapping and no-mass configuration from Albrecht et al. (2020)

###############################################################################
# This script tests ISMIP6 oceanic models in 2D from Jourdain et al. (2020):
# - ocean.ismip6: local parameterisation
# - ocean.ismip6nl: non-local parameterisation (needs drainage basin mask)
###############################################################################

SCRIPTNAME="#(run_paleo.sh)"

set -e  # exit on error

echo "$SCRIPTNAME   Run 30km grid paleo-climate spinup script using -pik options"
echo "$SCRIPTNAME   Run as './run_paleo.sh NN' for NN procs and 15km grid"

# get user and platform-specific variables like working_dir, pismcodedir,
# pism_exec and mpi command
source set_environment.sh

runname=`echo $PWD | awk -F/ '{print $NF}'`
#echo $runname
#codever=pism0.7_pik
thisdir=`echo $PWD`
outdir=$working_dir
PISM_EXEC=$pism_exec
#echo $outdir


###############################################################################
# INPUT DATA
export origfile=$input_data_dir/bedmap2/bedmap2_martos_MAR_monthly_pism30km.nc
export atmfile=$origfile
export oceanfile=$input_data_dir/schmidtko/schmidtko_pism30km_means_ismip6.nc

# PISM CONFIGURATION FILE
cp ./pism_config_ant_paleo.cdl $input_data_dir
export paramfile=$input_data_dir/pism_config_ant_paleo
ncgen3 -o ${paramfile}.nc ${paramfile}.cdl


# parallelization #############################################################
NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "run_paleo.sh 8" then NN = 8
  NN="$1"
fi
echo "$SCRIPTNAME              NN = $NN"
set -e  # exit on error


###### use only MPI if job is submitted #######################################
#if [ -n "${PISM_ON_CLUSTER:+1}" ]; then  # check if env var is set
  echo "This run was submitted, use MPI"
  PISM_MPIDO=$pism_mpi_do


#check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO=""
fi
DO=$PISM_DO

echo "PISM_MPIDO = $PISM_MPIDO"
PISM_DO="$PISM_MPIDO $NN $PISM_EXEC"


###############################################################################
source set_physics.sh

GRID=$THIRTYKMGRID
SKIP=$SKIPTHIRTYKM
GRIDNAME=30km

echo "TEST STOPS HERE"
###### output settings ########################################################
bootlength=100
nomasslength=200

# #############################################################################
# bootstrap and SHORT smoothing run to 20 years
# #############################################################################
stage=boot
INNAME=$origfile
RESNAMEBOOT=${outdir}/${stage}_${GRIDNAME}.nc
EXTRANAME=${outdir}/extra_${stage}_${GRIDNAME}.nc
expnomass="-extra_times 0:10:$bootlength -extra_vars bmelt,shelfbtemp,shelfbmassflux"

RUNTIME=$bootlength
echo
echo "$SCRIPTNAME  bootstrapping on $GRIDNAME grid plus SIA run for $RUNTIME years"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i ${INNAME} -bootstrap $GRID \
        $stress_opts $pre_calv_opts $pre_bed_opts \
      	$pre_ocean_opts $pre_atm_opts \
      	$tech_opts $config_opts \
        -extra_file $EXTRANAME $expnomass \
        -y $RUNTIME -o $RESNAMEBOOT" #-skip -skip_max $SKIP
#echo $DO $cmd
$DO $cmd >> ${outdir}/log/out_${stage}_${GRIDNAME}.out
exit # <-- uncomment to stop here


# #############################################################################
# run with -no_mass (no surface change) on coarse grid for 200yr
# #############################################################################
stage=nomass
INNAME=$RESNAMEBOOT
RESNAMENOMASS=${outdir}/${stage}_${GRIDNAME}.nc
TSNAME=${outdir}/ts_${stage}_${GRIDNAME}.nc
RUNTIME=$nomasslength
EXTRANAME=${outdir}/extra_${stage}_${GRIDNAME}.nc
expnomass="-extra_times 0:10:$nomasslength -extra_vars bmelt,shelfbtemp,shelfbmassflux"
echo
echo "$SCRIPTNAME  -no_mass (no surface change) SIA for $RUNTIME years"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i $INNAME -no_mass  \
     $stress_opts $pre_calv_opts $pre_bed_opts \
     $pre_ocean_opts $pre_atm_opts \
     $tech_opts $config_opts \
     -ys 0 -y $RUNTIME \
     -extra_file $EXTRANAME $expnomass \
     -o $RESNAMENOMASS"
#echo $DO $cmd
#$DO $cmd >> ${outdir}/log/out_${stage}_${GRIDNAME}.out
#exit # <-- uncomment to stop here
