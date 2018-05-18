#!/bin/bash

# Copyright (C) 2010-2015, 2018 Andy Aschwanden


if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(psg_flowline.sh)"
fi

echo
echo "# =================================================================================="
echo "# PISM Storglaciaren Flow Line Model"
echo "# =================================================================================="
echo

set -e # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "psg_flowline.sh 8" then NN = 8
  NN="$1"
fi

echo "$SCRIPTNAME              NN = $NN"

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
  PISM_EXEC="pismr"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi

echo

PCONFIG=psg_config.nc

# cat prefix and exec together
PISM="${PISM_BIN}${PISM_EXEC} -config_override $PCONFIG -o_order zyx"


DATANAME=storglaciaren_flowline.nc
PISM_DATANAME=pism_$DATANAME
INNAME=$PISM_DATANAME

# coupler settings
COUPLER="-surface given -surface_given_file $PISM_DATANAME"
COUPLER_FORCING="-surface given,forcing -surface_given_file $PISM_DATANAME"

# force-to-thickness
FTALPHA=0.05

# grid parameters
FINEGRID="-periodicity y -Mx 792 -My 3 -Mz 201 -Lz 300 -z_spacing equal"  # 5 m grid
FS=5
FINESKIP=5000
COARSEGRID="-periodicity y -Mx 114 -My 3 -Mz 101 -Lz 500 -z_spacing equal"  # 35 m grid
CS=35
COARSESKIP=1000

GRID=$COARSEGRID
SKIP=$COARSESKIP
GS=$CS
echo ""
if [ $# -gt 1 ] ; then
  if [ $2 -eq "2" ] ; then  # if user says "psg_flowline.sh N 1" then use 5m grid:
    echo "$SCRIPTNAME grid: ALL RUNS ON $FS m"
    echo "$SCRIPTNAME       WARNING: VERY LARGE COMPUTATIONAL TIME"
    GRID=$FINEGRID
    SKIP=$FINESKIP
    GS=$FS
  fi
else
    echo "$SCRIPTNAME grid: ALL RUNS ON $CS m"
fi
echo ""



EB="-sia_e 0.3"

#PARAMS="-pseudo_plastic -pseudo_plastic_q 0.75 -pseudo_plastic_uthreshold 10.0 -yield_stress mohr_coulomb -plastic_phi 40.0"
PARAMS="-pseudo_plastic -pseudo_plastic_q 0.25 -pseudo_plastic_uthreshold 100.0 -yield_stress mohr_coulomb -plastic_phi 40.0"

PETSCSTUFF="-ssafd_ksp_type gmres -ssafd_ksp_norm_type unpreconditioned -ssafd_ksp_pc_side right -ssafd_pc_type asm -ssafd_sub_pc_type lu"


FULLPHYS="-stress_balance ssa+sia $PARAMS $PETSCSTUFF"

SMOOTHRUNLENGTH=1
NOMASSRUNLENGTH=500

STEP=1

EXVARS="enthalpybase,temppabase,tempicethk,bmelt,tillwat,usurf,velsurf_mag,mask,hardav,thk" # add mask, so that check_stationarity.py ignores ice-free areas.

PREFIX=psg_flowline_

# bootstrap and do smoothing run to 1 year
OUTNAME=$PREFIX${GS}m_pre$SMOOTHRUNLENGTH.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run for ${SMOOTHRUNLENGTH}a"
cmd="$PISM_MPIDO $NN $PISM $EB -skip -skip_max $SKIP -i $INNAME -bootstrap $GRID \
  $COUPLER -y ${SMOOTHRUNLENGTH} -o $OUTNAME"
$PISM_DO $cmd

# run with -no_mass (no surface change) 
INNAME=$OUTNAME
OUTNAME=$PREFIX${GS}m_steady.nc
EXNAME=ex_${OUTNAME}
EXTIMES=0:25:${NOMASSRUNLENGTH}
echo
echo "$SCRIPTNAME  -no_mass (no surface change) sia run to achieve approximate enthalpy equilibrium, for ${NOMASSRUNLENGTH}a"
cmd="$PISM_MPIDO $NN $PISM $EB -i $INNAME $COUPLER \
  -no_mass -y ${NOMASSRUNLENGTH} \
  -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES -o $OUTNAME"
$PISM_DO $cmd


# use the force-to-thickness mechanism to infer the mass balance
STARTYEAR=0
RUNLENGTH=10
ENDTIME=$(($STARTYEAR + $RUNLENGTH))
INNAME=$OUTNAME
OUTNAME=ssa_ftt_${RUNLENGTH}a.nc
OUTNAMEFULL=$PREFIX${GS}m_$OUTNAME
TSNAME=ts_${OUTNAME}
TSTIMES=$STARTYEAR:$STEP:$ENDTIME

# create the force-to-thickness mask
ncap2 -O -s "ftt_mask = 0*mask+1" $INNAME $INNAME

echo
echo "$SCRIPTNAME  SSA run with force-to-thickness for $RUNLENGTH years on ${GS}m grid"
cmd="$PISM_MPIDO $NN $PISM -bootstrap $GRID -i $INNAME -regrid_file $INNAME $EB -skip -skip_max $SKIP $COUPLER_FORCING $FULLPHYS\
     -force_to_thickness_file $INNAME -force_to_thickness_alpha $FTALPHA \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -ys $STARTYEAR -y $RUNLENGTH -o_size big -o $OUTNAMEFULL"
$PISM_DO $cmd
echo
$PISM_DO flowline.py -c -o $OUTNAME $OUTNAMEFULL


COUPLER_ELEV="-surface elevation -ice_surface_temp -6,0,1395,1400 -climatic_mass_balance -3,2.5,1200,1450,1615 -climatic_mass_balance_limits -3,0"

STARTYEAR=0
RUNLENGTH=25
ENDTIME=$(($STARTYEAR + $RUNLENGTH))
INNAME=$OUTNAMEFULL
OUTNAME=ssa_${RUNLENGTH}a.nc
OUTNAMEFULL=$PREFIX${GS}m_$OUTNAME
TSNAME=ts_${OUTNAMEFULL}
EXNAME=ex_${OUTNAMEFULL}
TSTIMES=$STARTYEAR:$STEP:$ENDTIME
EXTIMES=$STARTYEAR:$STEP:$ENDTIME

echo
echo "$SCRIPTNAME  SSA run with elevation-dependent mass balance for $RUNLENGTH years on ${GS}m grid"
cmd="$PISM_MPIDO $NN $PISM $EB -skip -skip_max $SKIP -i $INNAME $COUPLER_ELEV $FULLPHYS \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTYEAR -y $RUNLENGTH -o_size big -o $OUTNAMEFULL"
$PISM_DO $cmd
echo
$PISM_DO flowline.py -c -o $OUTNAME $OUTNAMEFULL
