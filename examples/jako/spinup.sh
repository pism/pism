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
if [ -n "${PISM_BIN:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME     PISM_BIN = $PISM_BIN  (already set)"
else
  PISM_BIN=""    # just a guess
  echo "$SCRIPTNAME     PISM_BIN = $PISM_BIN"
fi

# set PISM_EXEC if using different executables, for example:
#  $ export PISM_EXEC="pismr -energy cold"
if [ -n "${PISM_EXEC:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC  (already set)"
else
  PISM_EXEC="pismr -regional"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi

PISM=$PISM_EXEC

BOOT=jako.nc
CLIMATEFILE=g5km_climate.nc
BCFILE=g5km_bc.nc

CLIMATE="-surface given,forcing -surface_given_file $CLIMATEFILE -force_to_thickness_file $BOOT"

# regarding physics: '-till_effective_fraction_overburden 0.02' plus
#   '-pseudo_plastic -pseudo_plastic_q 0.25' plus '-tauc_slippery_grounding_lines'
#   matches default choices in spinup.sh in examples/std-greenland/
# but here we do not use -sia_e 3.0 but instead -sia_e 1.0
PHYS="-front_retreat_file $BOOT -pik -sia_e 1.0 -stress_balance ssa+sia -topg_to_phi 15.0,40.0,-300.0,700.0 -till_effective_fraction_overburden 0.02 -pseudo_plastic -pseudo_plastic_q 0.25 -tauc_slippery_grounding_lines"

SKIP=5

LENGTH=${SPINUP_LENGTH:-2000}   # model years
EXDT=20    # 20 year between saves, thus 100 frames

cmd="$PISM_MPIDO $NN $PISM -i $BOOT -bootstrap -no_model_strip 10 \
  -Mx $Mx -My $My -Lz 4000 -Lbz 1000 -Mz 201 -Mbz 51 -z_spacing equal \
  -no_model_strip 10 $PHYS \
  -regional.zero_gradient \
  -spatial_file ex_spunjako_0.nc -spatial_times -$LENGTH:$EXDT:0 \
  -spatial_vars thk,velbase_mag,tillwat,tauc,dHdt,hardav,velsurf_mag,temppabase,diffusivity,bmelt,tempicethk_basal \
  -ts_file ts_spunjako_0.nc -ts_times -$LENGTH:yearly:0 \
  -ssa_dirichlet_bc -regrid_file $BCFILE -regrid_vars basal_melt_rate_grounded,tillwat,enthalpy,litho_temp,vel_bc \
  $CLIMATE -ys -$LENGTH -ye 0 -skip -skip_max $SKIP -o spunjako_0.nc"
$PISM_DO $cmd

# NOTES:
# useful diagnostic:  -ssa_view_nuh
# good postprocess:
#ncap -O -s "dtau=taud_mag-tauc" spunjako_0.nc spunjako_0.nc 

