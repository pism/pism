#!/bin/bash

# do a basic Jakoshavn spinup by:
#   ./spinup NN Mx My >> out.spin &

# if user says "spinup.sh 8 620 425" then NN=8, Mx=620, My=425
NN="$1"
Mx="$2"
My="$3"

if [ $# -lt 3 ] ; then  
  echo "spinup.sh ERROR: needs three arguments"
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

PISM=$PISM_EXEC

BOOT=jako.nc
CLIMATEFILE=g5km_climate.nc
BCFILE=g5km_bc.nc

CLIMATE="-surface given,forcing -surface_given_file $CLIMATEFILE -force_to_thk $BOOT"

# regarding physics: '-hydrology_pressure_fraction 0.98 -pseudo_plastic -pseudo_plastic_q 0.25' matches
#   spinup.sh in examples/searise-greenland/spinup.sh.  but here we do not use -sia_e 3.0
#   and we use a particular calving model that is different
PHYS="-calving ocean_kill -ocean_kill_file $BOOT -cfbc -kill_icebergs -diffuse_bwat -thk_eff -sia_e 1.0 -ssa_sliding -topg_to_phi 5.0,30.0,-300.0,700.0 -hydrology_pressure_fraction 0.98 -pseudo_plastic -pseudo_plastic_q 0.25"

SKIP=5

LENGTH=2000   # model years
EXDT=20    # 20 year between saves, thus 100 frames

cmd="$PISM_MPIDO $NN $PISM -boot_file $BOOT -no_model_strip 10 \
  -Mx $Mx -My $My -Lz 4000 -Lbz 1000 -Mz 201 -Mbz 51 -z_spacing equal \
  -no_model_strip 10 $PHYS \
  -extra_file ex_spunjako_0.nc -extra_times -$LENGTH:$EXDT:0 \
  -extra_vars thk,cbase,bwp,tauc,dhdt,hardav,csurf,temppabase,diffusivity,bmelt,tempicethk_basal \
  -ts_file ts_spunjako_0.nc -ts_times -$LENGTH:yearly:0 \
  -ssa_dirichlet_bc -regrid_file $BCFILE -regrid_vars bmelt,bwat,enthalpy,litho_temp,vel_ssa_bc \
  $CLIMATE -ys -$LENGTH -ye 0 -skip -skip_max $SKIP -o spunjako_0.nc"
$PISM_DO $cmd

# NOTES:
# useful diagnostic:  -ssa_view_nuh
# good postprocess:
#ncap -O -s "dtau=taud_mag-tauc" spunjako_0.nc spunjako_0.nc 

