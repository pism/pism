#!/bin/bash

# uses 2 files:
#    pism_Greenland_5km_v1.1.nc  from examples/std-greenland/
#    g2km_gridseq.nc             ditto; or similar
# these files are documented in the PISM User's Manual, chapter 1
#
# do
#   $ ./run-decoupled.sh 5 g2km_gridseq
# for 5 year run

set -e  # exit on error

# check if env var PISM_DO was set (i.e. set PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  echo "#   PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO=""
fi

DURATION=$1
INNAME=$2

MPIDO="mpiexec -n 6"

CLIMATE="-surface given -surface_given_file pism_Greenland_5km_v1.1.nc"
PHYS="-sia_e 3.0 -stress_balance ssa+sia -topg_to_phi 15.0,40.0,-300.0,700.0 -pseudo_plastic -pseudo_plastic_q 0.5 -till_effective_fraction_overburden 0.02 -tauc_slippery_grounding_lines"
CALVING="-calving ocean_kill -ocean_kill_file pism_Greenland_5km_v1.1.nc"

# run this to check for no shock: continue g2km_gridseq.nc run
NAME=cont.nc
cmd="$MPIDO pismr -i $INNAME -skip -skip_max 20 $CLIMATE $PHYS $CALVING -ts_file ts_$NAME -ts_times 0:yearly:$DURATION -y $DURATION -o $NAME"
#$PISM_DO $cmd
echo

# suitable for -hydrology routing,distributed runs which are decoupled:
EXVAR="mask,thk,topg,usurf,tillwat,bwat,hydrobmelt,bwatvel"
EXVARDIST="${EXVAR},bwp,bwprel,hydrovelbase_mag"

# -hydrology routing
NAME=routing-decoupled.nc
cmd="$MPIDO pismr -i $INNAME -no_mass -energy none -stress_balance none $CLIMATE -extra_file ex_$NAME -extra_times 0:monthly:$DURATION -extra_vars $EXVAR -ts_file ts_$NAME -ts_times 0:daily:$DURATION -hydrology routing -hydrology_bmelt_file $INNAME -ys 0 -y $DURATION -max_dt 0.03 -o $NAME"
$PISM_DO $cmd
echo

# -hydrology distributed
NAME=distributed-decoupled.nc
cmd="$MPIDO pismr -i $INNAME -no_mass -energy none -stress_balance none $CLIMATE -extra_file ex_$NAME -extra_times 0:monthly:$DURATION -extra_vars $EXVARDIST -ts_file ts_$NAME -ts_times 0:daily:$DURATION -hydrology distributed -hydrology_bmelt_file $INNAME -hydrology_velbase_mag_file $INNAME -ys 0 -y $DURATION -max_dt 0.03 -o $NAME"
$PISM_DO $cmd
echo

