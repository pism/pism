#!/bin/bash

# do a basic Jakoshavn model run, presumably on a fine grid; run as
#   ./century NN Mx My >> out.century &

# the grid DEFAULTS below apply to:
#   $ ncdump -h jako.nc |head -n 4
#   netcdf jako {
#   dimensions:
#       y = 425 ;
#       x = 620 ;

NN=4     # default number of MPI processes
if [ $# -gt 0 ] ; then  # if user says "century.sh 8" then NN=8
  NN="$1"
fi

Mx=620   # grid default
if [ $# -gt 1 ] ; then  # if user says "century.sh 8 620" then NN=8 and Mx=620
  Mx="$2"
fi

My=425   # grid default
if [ $# -gt 2 ] ; then  # if user says "century.sh 8 620 425" then NN=8, Mx=620, My=425
  My="$3"
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
PREFILE=jako3km_0.nc
CLIMATEFILE=g5km_climate.nc
BCFILE=g5km_bc.nc

CLIMATE="-surface given,forcing -surface_given_file $CLIMATEFILE -force_to_thk $BOOT"

# assert: PHYS is same as in spinup.sh
PHYS="-ocean_kill $BOOT -cfbc -kill_icebergs -topg_to_phi 5.0,30.0,-300.0,700.0 -diffuse_bwat -thk_eff -ssa_sliding -plastic_pwfrac 0.98 -pseudo_plastic -pseudo_plastic_q 0.25"

SKIP=10

LENGTH=100   # model years
EXDT=1       # 1 year between saves, thus 100 frames

# recommended if sufficient memory (30GB total?)
Mz=401
Mbz=101
# lower memory choice for demo only (works if 12GB)
#Mz=101
#Mbz=26

echo
cmd="$PISM_MPIDO $NN $PISM -boot_file $BOOT  \
  -Mx $Mx -My $My -Lz 4000 -Lbz 1000 -Mz $Mz -Mbz $Mbz -z_spacing equal \
  -no_model_strip 10 $PHYS \
  -ssa_dirichlet_bc -regrid_file $PREFILE -regrid_vars thk,Href,bmelt,bwat,enthalpy,litho_temp,vel_ssa_bc \
  $CLIMATE -y 0.01 -skip -skip_max $SKIP -o jako1km_short.nc"
$PISM_DO $cmd

echo
cmd="$PISM_MPIDO $NN $PISM -i jako1km_short.nc \
  -no_model_strip 10 $PHYS \
  -extra_file ex_jako1km.nc -extra_times -$LENGTH:$EXDT:0 \
  -extra_vars thk,cbase,bwp,tauc,dhdt,hardav,csurf,temppabase,diffusivity,bmelt,tempicethk_basal \
  -ts_file ts_jako1km.nc -ts_times -$LENGTH:monthly:0 \
  -ssa_dirichlet_bc -regrid_file $BCFILE -regrid_vars vel_ssa_bc \
  $CLIMATE -ys 0 -ye $LENGTH -skip -skip_max $SKIP -o jako1km.nc"
$PISM_DO $cmd
