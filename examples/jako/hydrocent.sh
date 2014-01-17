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
  PISM_EXEC="pismo"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi


BOOT=jako.nc
CLIMATEFILE=g5km_climate.nc
BCFILE=g5km_bc.nc

CLIMATE="-surface given,forcing -surface_given_file $CLIMATEFILE -force_to_thk $BOOT"

PHYS="-calving ocean_kill -ocean_kill_file $BOOT -cfbc -kill_icebergs -sia_e 1.0 -ssa_sliding -topg_to_phi 5.0,30.0,-300.0,700.0 -hydrology_pressure_fraction 0.98 -pseudo_plastic -pseudo_plastic_q 0.25"

SKIP=10

LENGTH=100   # model years

Mz=201
Mbz=51
# strongly recommended if sufficient memory:
#Mz=401
#Mbz=101

echo
cmd="$PISM_MPIDO $NN $PISM_EXEC -boot_file $BOOT  \
  -Mx $Mx -My $My -Lz 4000 -Lbz 1000 -Mz $Mz -Mbz $Mbz -z_spacing equal \
  -no_model_strip 10 $PHYS \
  -ssa_dirichlet_bc -regrid_file $PREFILE -regrid_vars thk,Href,bmelt,bwat,enthalpy,litho_temp,vel_ssa_bc \
  $CLIMATE -y 0.01 -skip -skip_max $SKIP -o jakofine_short.nc"
#$PISM_DO $cmd

# insert a hydrology test case; this is hijacking the main purpose of century.sh ...

HYDROLENGTH=0.002
#HYDRO="-hydrology distributed -init_P_from_steady -hydrology_null_strip 10 -report_mass_accounting"
HYDRO="-hydrology distributed -hydrology_null_strip 10 -report_mass_accounting"

echo
cmd="$PISM_MPIDO $NN $PISM_EXEC -i jakofine_short.nc \
  -no_model_strip 10 $PHYS $HYDRO \
  -extra_file ex_jakofine_hydro.nc -extra_times 0:0.0005:$HYDROLENGTH \
  -extra_vars thk,cbase,bwat,bwatvel,bwp,effbwp,enwat,tauc,dhdt,hardav,csurf,temppabase,diffusivity,bmelt,tempicethk_basal \
  -ts_file ts_jakofine_hydro.nc -ts_times 0:0.0005:$HYDROLENGTH \
  -ssa_dirichlet_bc -regrid_file $BCFILE -regrid_vars vel_ssa_bc \
  $CLIMATE -ys 0 -ye $HYDROLENGTH -o jakofine_hydro.nc"
$PISM_DO $cmd

exit


echo
cmd="$PISM_MPIDO $NN $PISM_EXEC -i jakofine_short.nc \
  -no_model_strip 10 $PHYS \
  -extra_file ex_jakofine.nc -extra_times 0:yearly:$LENGTH \
  -extra_vars thk,cbase,bwat,tauc,dhdt,hardav,csurf,temppabase,diffusivity,bmelt,tempicethk_basal \
  -ts_file ts_jakofine.nc -ts_times 0:monthly:$LENGTH \
  -ssa_dirichlet_bc -regrid_file $BCFILE -regrid_vars vel_ssa_bc \
  $CLIMATE -ys 0 -ye $LENGTH -skip -skip_max $SKIP -o jakofine.nc"
$PISM_DO $cmd

