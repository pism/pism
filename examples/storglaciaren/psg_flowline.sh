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

PCONFIG=psg_config.nc

# cat prefix and exec together
PISM="${PISM_PREFIX}${PISM_EXEC} -cts -config_override $PCONFIG"


DATANAME=storglaciaren_flowline.nc
PISM_DATANAME=pism_$DATANAME
INNAME=$PISM_DATANAME

# coupler settings
COUPLER="-surface elevation"

# grid parameters
FINEGRID="-periodicity y -Mx 792 -My 3 -Mz 201 -Mbz 31 -Lz 500 -Lbz 300 -z_spacing equal -zb_spacing equal"  # 5 m grid
FS=5
FINESKIP=1500
COARSEGRID="-periodicity y -Mx 114 -My 3 -Mz 51 -Mbz 31 -Lz 1000 -Lbz 300 -z_spacing equal -zb_spacing equal"  # 35 m grid
CS=35
COARSESKIP=500

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


phi=40
e=0.3
q=0.75
uth=10

TILLPHI="-plastic_phi $phi"
EB="-e $e"
PARAMS="$TILLPHI -pseudo_plastic_q $q -plastic_pwfrac 0.98 -pseudo_plastic_uthreshold $uth"
# Use the FEM solver for the SSA, as the FD solver show convergence issues
FULLPHYS="-ssa_method fem -ssa_sliding -thk_eff $PARAMS"


SMOOTHRUNLENGTH=1
NOMASSRUNLENGTH=250

TSTIMESTEP=1

EXVARS="enthalpybase,temppabase,tempicethk,bmelt,bwat,usurf,csurf,mask,hardav" # add mask, so that check_stationarity.py ignores ice-free areas.

PREFIX=psg

# bootstrap and do smoothing run to 1 year
OUTNAME=$PREFIX${GS}m_pre$SMOOTHRUNLENGTH.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run for ${SMOOTHRUNLENGTH}a"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -boot_file $INNAME $GRID \
  $COUPLER -y ${SMOOTHRUNLENGTH} -o $OUTNAME"
$PISM_DO $cmd

# run with -no_mass (no surface change) 
INNAME=$OUTNAME
OUTNAME=$PREFIX${GS}m_steady.nc
EXNAME=ex_${OUTNAME}
EXTIMES=0:25:${NOMASSRUNLENGTH}
echo
echo "$SCRIPTNAME  -no_mass (no surface change) sia run to achieve approximate enthalpy equilibrium, for ${NOMASSRUNLENGTH}a"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -i $INNAME $COUPLER \
  -no_mass -y ${NOMASSRUNLENGTH} \
  -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES -o $OUTNAME"
$PISM_DO $cmd


STARTYEAR=0
RUNLENGTH=100
ENDTIME=$(($STARTYEAR + $RUNLENGTH))
INNAME=$OUTNAME
OUTNAME=ssa.nc
OUTNAMEFULL=$PREFIX${GS}m_$OUTNAME
TSNAME=ts_${OUTNAME}
TSTIMES=$STARTYEAR:$TSTIMESTEP:$ENDTIME
echo
echo "$SCRIPTNAME  SIA run with force-to-thickness for $RUNLENGTH years on ${GS}m grid"
cmd="$PISM_MPIDO $NN $PISM $EB -skip $SKIP -i $INNAME $COUPLER $FULLPHYS \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -ys $STARTYEAR -y $RUNLENGTH -o_size big -o $OUTNAMEFULL"
$PISM_DO $cmd
$PISM_DO flowline.py -c -o $OUTNAME $OUTNAMEFULL
