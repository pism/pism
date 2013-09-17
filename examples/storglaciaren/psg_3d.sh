#!/bin/bash

# Copyright (C) 2011, 2012 Andy Aschwanden and Ed Bueler

set -e # exit on error

echo "# PISM Storglaciaren 3d Model"

if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(psg_3d.sh)"
fi

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

PCONFIG=psg_config.nc

# cat prefix and exec together
PISM="${PISM_PREFIX}${PISM_EXEC} -config_override $PCONFIG"


DATANAME=storglaciaren_3d.nc
PISM_DATANAME=pism_$DATANAME
INNAME=$PISM_DATANAME

EB="-e_sia 0.3"
PARAMS="-plastic_phi"
FULLPHYS="-ssa_sliding -thk_eff $PARAMS"
COUPLER="-surface given" # FIXME  should be using PSElevation as in flowline example
COUPLER_FORCING="-surface given,forcing -surface_given_file $PISM_DATANAME"

# 40 m grid
GRID="-Mx 94 -My 51 -Mz 151 -Mbz 1 -Lz 300 -z_spacing equal"
GS=40
SKIP=200

# 20 m grid
#GRID="-Mx 186 -My 101 -Mz 301 -Mbz 1 -Lz 300 -z_spacing equal"
#GS=20
#SKIP=500

# 10 m grid
#GRID="-Mx 371 -My 201 -Mz 301 -Mbz 1 -Lz 300 -z_spacing equal"
#GS=10
#SKIP=500

# bootstrap and do smoothing run
OUTNAME=psg_3d_${GS}m_pre1.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run for 1 a"
cmd="$PISM_MPIDO $NN $PISM -skip $SKIP -boot_file $INNAME $GRID \
  $COUPLER -y 1 -o $OUTNAME"
$PISM_DO $cmd

# FIXME:  reasonable to run SIA somewhat longer to equilibrium to establish stable thermodynamical state before any attempt to invert surface velocities; this is the start of it
# FIXME:  also reasonable to use hybrid model with ad hoc specification of basal yield stress
INNAME=$OUTNAME
#incSTEADYNAME=inc_psg_3d_${GS}m_Tsteady.nc  # FIXME: missing a field because of bug
STEADYNAME=psg_3d_${GS}m_steady.nc
echo
echo "$SCRIPTNAME  running toward thermodynamical steady state"
cmd="$PISM_MPIDO $NN $PISM -i $INNAME \
  $COUPLER -y 500 -no_mass -o $STEADYNAME"
$PISM_DO $cmd

echo
echo "$SCRIPTNAME  done"

# We use the force-to-thickness mechanism to infer the mass balance

STARTYEAR=0
RUNLENGTH=100
ENDTIME=$(($STARTYEAR + $RUNLENGTH))
INNAME=$OUTNAME
OUTNAME=ssa_ftt_${RUNLENGTH}a.nc
OUTNAMEFULL=$PREFIX${GS}m_$OUTNAME
TSNAME=ts_${OUTNAME}
TSTIMES=$STARTYEAR:$STEP:$ENDTIME
echo
echo "$SCRIPTNAME  SSA run with force-to-thickness for $RUNLENGTH years on ${GS}m grid"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -i $INNAME $COUPLER_FORCING $FULLPHYS\
     -force_to_thk $PISM_DATANAME \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -ys $STARTYEAR -y $RUNLENGTH -o_size big -o $OUTNAMEFULL"
$PISM_DO $cmd
echo
