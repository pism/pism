#!/bin/bash

# Copyright (C) 2010-2011 Andy Aschwanden


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

set -e  # exit on error

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
  PISM_EXEC="pismr"
  echo "$SCRIPTNAME       PISM_EXEC = $PISM_EXEC"
fi

echo

PCONFIG=sg_config.nc

#EB="-ice_type custom -ice_custom_softness 2.2e-24"        # "standard" value for Storglaciaren
EB="-e 0.3" #ice type options, enhancement factor
#EB="-e 1"

# cat prefix and exec together
PISM="${PISM_PREFIX}${PISM_EXEC} -cts -config_override $PCONFIG"


DATANAME=sg_flowline.nc
PISM_DATANAME=pism_$DATANAME
INNAME=$PISM_DATANAME

# coupler settings
COUPLER_SIMPLE="-surface constant"
COUPLER_FORCING="-surface constant,forcing"


# some options to make it converge. The latter always requires to have petsc compiled with the parallel
# direct solver "mumps"
PETSCSTUFF="-pc_type asm -sub_pc_type lu -ksp_type lgmres -ksp_right_pc"
#PETSCSTUFF="-pc_type lu -pc_factor_mat_solver_package mumps"

# grid parameters
FINEGRID="-periodicity y -Mx 792 -My 3 -Mz 201 -Mbz 31 -Lz 300 -Lbz 300 -z_spacing equal -zb_spacing equal"  # 5 m grid
FS=5
FINESKIP=750
COARSEGRID="-periodicity y -Mx 114 -My 3 -Mz 51 -Mbz 31 -Lz 300 -Lbz 300 -z_spacing equal -zb_spacing equal"  # 35 m grid
CS=35
COARSESKIP=250

# force-to-thickness
FTALPHA=0.05

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




SMOOTHRUNLENGTH=1
NOMASSRUNLENGTH=250

EXVARS="enthalpybase,temppabase,tempicethk,bmelt,bwat,usurf,csurf,mask,hardav" # add mask, so that check_stationarity.py ignores ice-free areas.

PREFIX=psg

# bootstrap and do smoothing run to 1 year
PRE0NAME=$PREFIX${GS}m_pre$SMOOTHRUNLENGTH.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run (for ${SMOOTHRUNLENGTH}a)"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -boot_file $INNAME $GRID \
  $COUPLER_SIMPLE -y ${SMOOTHRUNLENGTH} -o $PRE0NAME"
$PISM_DO $cmd

# run with -no_mass (no surface change) for 
PRE1NAME=$PREFIX${GS}m_steady.nc
EX1NAME=ex_${PRE1NAME}
EXTIMES=0:25:${NOMASSRUNLENGTH}
echo
echo "$SCRIPTNAME  -no_mass (no surface change) SIA run to achieve approximate enthalpy equilibrium, for ${NOMASSSIARUNLENGTH}a"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIPTWENTYKM -i $PRE0NAME $COUPLER_SIMPLE \
  -no_mass -y ${NOMASSRUNLENGTH} \
  -extra_file $EX1NAME -extra_vars $EXVARS -extra_times $EXTIMES -o $PRE1NAME"
$PISM_DO $cmd

# smoothing for 1 year
PRE2NAME=$PREFIX${GS}m_post$SMOOTHRUNLENGTH.nc
echo
echo "$SCRIPTNAME  smoothing with SIA for ${SMOOTHRUNLENGTH}a"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -i $PRE1NAME $COUPLER_SIMPLE -y $SMOOTHRUNLENGTH -o_size big -o $PRE2NAME"
$PISM_DO $cmd


STARTYEAR=0
RUNLENGTH=250

TSTIMESTEP=1

ENDTIME=$(($STARTYEAR + $RUNLENGTH))
INNAME=$PRE2NAME
OUTNAME=$PREFIX${GS}m_sia_steadystate.nc
TSNAME=ts_${OUTNAME}
TSTIMES=$STARTYEAR:$TSTIMESTEP:$ENDTIME
echo
echo "$SCRIPTNAME  SIA steady-state run for $RUNLENGTH years on ${GS}m grid"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -i $INNAME $COUPLER_SIMPLE  \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -ys $STARTYEAR -y $RUNLENGTH -o_size big -o $OUTNAME"
$PISM_DO $cmd


RUNLENGTH=50
ENDTIME=$(($STARTYEAR + $RUNLENGTH))
INNAME=$PRE2NAME
OUTNAME=sia_ftt.nc
OUTNAMEFULL=$PREFIX${GS}m_$OUTNAME
TSNAME=ts_${OUTNAME}
TSTIMES=$STARTYEAR:$TSTIMESTEP:$ENDTIME
echo
echo "$SCRIPTNAME  SIA run with force-to-thickness for $RUNLENGTH years on ${GS}m grid"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -i $INNAME $COUPLER_FORCING  \
     -force_to_thk $INNAME -force_to_thk_alpha $FTALPHA \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -ys $STARTYEAR -y $RUNLENGTH -o_size big -o $OUTNAMEFULL"
$PISM_DO $cmd
$PISM_DO flowline.py -c -o $OUTNAME $OUTNAMEFULL

ENDTIME=$(($STARTYEAR + $RUNLENGTH))
INNAME=$PRE2NAME
OUTNAME=$PREFIX${GS}m_ssa.nc
TSNAME=ts_${OUTNAME}
TSTIMES=$STARTYEAR:$TSTIMESTEP:$ENDTIME
echo
echo "$SCRIPTNAME  SSA run for $RUNLENGTH years on ${GS}m grid"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -i $INNAME $COUPLER_SIMPLE $TILLPHI $FULLPHYS \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -ys $STARTYEAR -y $RUNLENGTH -o_size big -o $OUTNAME"
$PISM_DO $cmd


#TILLPHI="-topg_to_phi 5,45,1300,1600,15"


phi=20
e=0.3
q=0.75
uth=10

TILLPHI="-plastic_phi $phi"
EB="-e $e"
PARAMS="-pseudo_plastic_q $q -plastic_pwfrac 0.98 -pseudo_plastic_uthreshold $uth"
FULLPHYS="-ssa_sliding -thk_eff $PARAMS $PETSCSTUFF"

RUNLENGTH=50
ENDTIME=$(($STARTYEAR + $RUNLENGTH))

OUTNAME=ssa_ftt_phi_${phi}_e_${e}_uth_${uth}_q_${q}.nc
OUTNAMEFULL=$PREFIX${GS}m_$OUTNAME
TSNAME=ts_${OUTNAME}
TSTIMES=$STARTYEAR:$TSTIMESTEP:$ENDTIME
echo
echo "$SCRIPTNAME  SSA run with force-to-thickness for $RUNLENGTH years on ${GS}m grid"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -i $INNAME  $COUPLER_FORCING $TILLPHI $FULLPHYS \
                     -force_to_thk $INNAME -force_to_thk_alpha $FTALPHA \
                     -ts_file $TSNAME -ts_times $TSTIMES \
                     -ys $STARTYEAR -y $RUNLENGTH -o_size big -o $OUTNAMEFULL"
$PISM_DO $cmd
$PISM_DO flowline.py -c -o $OUTNAME $OUTNAMEFULL
