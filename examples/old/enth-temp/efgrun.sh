#!/bin/bash

# Copyright (C) 2009-2011 Andy Aschwanden and Ed Bueler

# Enthalpy Formulation for Greenland (EFG) demonstrates the applicability
# of the enthalpy formulation to the Greenland Ice Sheet
#
# recommended way to run with N processors is " ./efgrun.sh N >& out.efg & "
# which gives a viewable (with "less", for example) transcript in out.efg

if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(efgrun.sh)"
fi

set -e  # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "efgrun.sh 8" then NN = 8
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

# executable
PISM_EXEC="pismr"

echo



# preprocess.sh generates eis_green_smoothed.nc; run it first
if [ -n "${PISM_DATANAME:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME   PISM_DATANAME = $PISM_DATANAME  (already set)"
else
  PISM_DATANAME=../eisgreen/eis_green_smoothed.nc
fi

PISM_CONFIG=../eisgreen/eismint_config.nc

for INPUT in $PISM_DATANAME $PISM_CONFIG; do
  if [ -e "$INPUT" ] ; then  # check if file exist
    echo "$SCRIPTNAME           input   $INPUT (found)"
  else
    echo "$SCRIPTNAME           input   $INPUT (MISSING!!)"
    echo
    echo "$SCRIPTNAME           please run preprocess.sh in ../eisgreen/, exiting"
    echo
    exit
  fi
done

INNAME=$PISM_DATANAME

EXVARS="enthalpybase,temppabase,tempicethk_basal,tempicethk,bmelt,bwat,usurf,csurf,mask,hardav,diffusivity" # add mask, so that check_stationarity.py ignores ice-free areas.
PREFIX=efg

CRUNLENGTH=230000
FRUNLENGTH=20000
EXTIMESTEP=1000
TSTIMESTEP=10

GRIDBASE="-Lz 4000 -Lbz 2000 -Mbz 51"
CGRID="${GRIDBASE} -Mx 83 -My 141 -Mz 101"   # coarse grid = 20km
FGRID="${GRIDBASE} -Mx 165 -My 281 -Mz 201"  # fine grid = 10km

for METHOD in "enth" "temp" "vark" "varc" "varck"; do

    if [ ${METHOD} == "temp" ]; then
        COEM="-cold"
        REGRIDVARS="temp,thk,litho_temp,bwat,bmelt"
    else
        COEM=
        REGRIDVARS="enthalpy,thk,litho_temp,bwat,bmelt"
    fi

    if [ ${METHOD} == "vark" ] ; then
        OPTIONS="-vark"
    elif [ ${METHOD} == "varck" ] ; then
        OPTIONS="-vark -varc"
    elif [ ${METHOD} == "varc" ] ; then
        OPTIONS="-varc"
    else
        OPTIONS=""
    fi

    # cat prefix and exec together
    PISM="${PISM_PREFIX}${PISM_EXEC} $OPTIONS -cts -atmosphere eismint_greenland -surface pdd -config_override $PISM_CONFIG"

    echo "$SCRIPTNAME      executable = '$PISM'"
    echo ""

    echo "$SCRIPTNAME ++++ doing 20km $METHOD run ++++"
    SKIP=10
    ENDTIME=$CRUNLENGTH
    OUTNAME=${PREFIX}_20km_${METHOD}.nc
    EXNAME=ex_$OUTNAME
    EXTIMES=0:$EXTIMESTEP:$ENDTIME
    TSNAME=ts_$OUTNAME
    TSTIMES=0:$TSTIMESTEP:$ENDTIME
    echo "$SCRIPTNAME  $METHOD run for $CRUNLENGTH years on 20km grid"
    cmd="$PISM_MPIDO $NN $PISM $COEM -skip $SKIP -boot_file $INNAME $CGRID \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys 0 -ye $ENDTIME -o_size big -o $OUTNAME"
    $PISM_DO $cmd
    echo

    echo "$SCRIPTNAME ++++ doing 10km $METHOD run ++++"
    SKIP=50
    STARTTIME=$ENDTIME
    ENDTIME=$(($STARTTIME + $FRUNLENGTH))
    REGRIDNAME=$OUTNAME
    OUTNAME=${PREFIX}_10km_${METHOD}.nc
    EXNAME=ex_$OUTNAME
    EXTIMES=$STARTTIME:$EXTIMESTEP:$ENDTIME
    TSNAME=ts_$OUTNAME
    TSTIMES=$STARTTIME:$TSTIMESTEP:$ENDTIME
    echo "$SCRIPTNAME  $METHOD run for $FRUNLENGTH years on 10km grid, regridding from 20km result"
    cmd="$PISM_MPIDO $NN $PISM $COEM -skip $SKIP -boot_file $INNAME $FGRID \
     -regrid_file $REGRIDNAME -regrid_vars $REGRIDVARS \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o_size big -o $OUTNAME"
    $PISM_DO $cmd   
    echo

done
