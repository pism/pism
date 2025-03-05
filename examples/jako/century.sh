#!/bin/bash

# first run spinup.sh to produce a starting file (e.g. spunjako_0.nc)
# here is a basic Jakoshavn model run:
#   ./century NN Mx My spunjako_0.nc >> out.century &

# Mx=620, My=425 is 1 km grid
# Mx=310, My=213 is (approx) 2 km grid
# Mx=125, My=86 is (approx) 5 km grid

# if user says "century.sh 8 620 425 foo.nc" then NN=8, Mx=620, My=425, PREFILE=foo.nc
NN="$1"
Mx="$2"
My="$3"
PREFILE="$4"

if [ $# -lt 4 ] ; then  
  echo "century.sh ERROR: needs four arguments"
  exit
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
if [ -n "${PISM_BIN:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_BIN = $PISM_BIN  (already set)"
else
  PISM_BIN=""    # just a guess
  echo "$SCRIPTNAME     PISM_BIN = $PISM_BIN"
fi

# set PISM_EXEC if using different executables, for example:
#  $ export PISM_EXEC="pism -energy cold"
if [ -n "${PISM_EXEC:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC  (already set)"
else
  PISM_EXEC="pism -regional"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi


BOOT=jako.nc
CLIMATEFILE=g5km_climate.nc
BCFILE=g5km_bc.nc

CLIMATE="-surface given,forcing -surface_given_file $CLIMATEFILE -force_to_thickness_file $BOOT"

# regarding physics: match choices in spinup.sh
PHYS="-front_retreat_file $BOOT -pik -sia_e 1.0 -stress_balance ssa+sia -topg_to_phi -phi_min 15.0 -phi_max 40.0 -topg_min -300.0 -topg_max 700.0 -till_effective_fraction_overburden 0.02 -pseudo_plastic -pseudo_plastic_q 0.25 -tauc_slippery_grounding_lines"

SKIP=10

LENGTH=${RUN_LENGTH:-100}   # model years

Mz=401
Mbz=101
# inferior, but use if insufficient memory
#Mz=201
#Mbz=51

echo
cmd="$PISM_MPIDO $NN $PISM_EXEC -i $BOOT -bootstrap  \
  -Mx $Mx -My $My -Lz 4000 -Lbz 1000 -Mz $Mz -Mbz $Mbz -z_spacing equal \
  -no_model_strip 10 $PHYS \
  -ssa_dirichlet_bc -regrid_file $PREFILE -regrid_vars thk,basal_melt_rate_grounded,tillwat,enthalpy,litho_temp,vel_bc \
  $CLIMATE -y 0.01 -o jakofine_short.nc"
$PISM_DO $cmd

echo
cmd="$PISM_MPIDO $NN $PISM_EXEC -i jakofine_short.nc \
  -no_model_strip 10 $PHYS \
  -extra_file ex_jakofine.nc -extra_times 0:yearly:$LENGTH \
  -extra_vars mask,thk,velbase_mag,tillwat,tauc,dhdt,hardav,velsurf_mag,temppabase,diffusivity,bmelt,tempicethk_basal \
  -ts_file ts_jakofine.nc -ts_times 0:monthly:$LENGTH \
  -ssa_dirichlet_bc -regrid_file $BCFILE -regrid_vars vel_bc \
  $CLIMATE -ys 0 -ye $LENGTH -skip -skip_max $SKIP -o jakofine.nc"
$PISM_DO $cmd

