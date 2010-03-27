#!/bin/bash

# Copyright (C) 2009-2010 Andy Aschwanden and Ed Bueler

# PISM SeaRISE Greenland
#
#
# before using this script, run preprocess.sh to download and adjust metadata
# on SeaRISE "Present Day Greenland" master dataset
#
# recommended way to run with N processors is " ./psearise-spinup.sh N >& out.psea & "
# which gives a viewable (with "less", for example) transcript in out.psea

if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(spinup.sh)"
fi

echo
echo "# =================================================================================="
echo "# PISM SeaRISE Greenland: spinup"
echo "# =================================================================================="
echo

set -e  # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "psearise.sh 8" then NN = 8
  NN="$1"
fi

echo "$SCRIPTNAME              NN = $NN"
set -e  # exit on error

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

# preprocess.sh generates pism_*.nc files; run it first
if [ -n "${PISM_DATANAME:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME   PISM_DATANAME = $PISM_DATANAME  (already set)"
else
  PISM_DATAVERSION=0.93
  PISM_DATANAME=pism_Greenland_5km_v$PISM_DATAVERSION.nc
fi
if [ -n "${PISM_TEMPSERIES:+1}" ] ; then
  echo "$SCRIPTNAME PISM_TEMPSERIES = $PISM_TEMPSERIES  (already set)"
else
  PISM_TEMPSERIES=pism_dT.nc
fi
if [ -n "${PISM_SLSERIES:+1}" ] ; then
  echo "$SCRIPTNAME   PISM_SLSERIES = $PISM_SLSERIES  (already set)"
else
  PISM_SLSERIES=pism_dSL.nc
fi

for INPUT in $PISM_DATANAME $PISM_TEMPSERIES $PISM_SLSERIES; do
  if [ -e "$INPUT" ] ; then  # check if file exist
    echo "$SCRIPTNAME           input   $INPUT (found)"
  else
    echo "$SCRIPTNAME           input   $INPUT (MISSING!!)"
    echo
    echo "$SCRIPTNAME           please run ./preprocess.sh, exiting"
    echo
    exit
  fi
done

INNAME=$PISM_DATANAME
OUTFINALNAME=g5km_0.nc

# run lengths and starting time for paleo
SMOOTHRUNLENGTH=100
NOMASSSIARUNLENGTH=50000
PALEOSTARTYEAR=-125000

# grids
TWENTYKMGRID="-Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 41 -Mbz 16"
TENKMGRID="-Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 31"
FIVEKMGRID="-Mx 301 -My 561 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 31"

# skips
SKIPTWENTYKM=2
SKIPTENKM=20
SKIPFIVEKM=50

# "standard" spinup 10km / 5km
COARSEGRID=$TWENTYKMGRID
FINEGRID=$TWENTYKMGRID
COARSESKIP=$SKIPTWENTYKM
FINESKIP=$SKIPTWENTYKM
CS=20 # km
FS=20 # km

COARSEENDTIME=-40000 # BP

echo ""
if [ $# -gt 1 ] ; then
  if [ $2 -eq "1" ] ; then  # if user says "spinup.sh N 1" then COARSE:
    echo "$SCRIPTNAME grid: ALL RUNS ON 10km"
    echo "$SCRIPTNAME       WARNING: LARGE COMPUTATIONAL TIME"
    COARSEGRID=$TENKMGRID
    FINEGRID=$TENKMGRID
    COARSESKIP=$SKIPTENKM
    FINESKIP=$SKIPTENKM
    CS=10 # km
    FS=10 # km
  fi
  if [ $2 -eq "2" ] ; then  # if user says "spinup.sh N 2" then SHORT PALEO and VERY COARSE:
    echo "$SCRIPTNAME grid: ON 10km TILL ${COARSEENDTIME}a THEN ON 5km"
    echo "$SCRIPTNAME       WARNING: VERY LARGE COMPUTATIONAL TIME"
    COARSEGRID=$TENKMGRID
    FINEGRID=$FIVEKMGRID
    COARSESKIP=$SKIPFIVEKM
    FINESKIP=$SKIPFIVEKM
    CS=10 # km
    FS=5  # km
  fi
else
    echo "$SCRIPTNAME grid: ALL RUNS ON 20km"
fi
echo ""

echo "$SCRIPTNAME     coarse grid = '$COARSEGRID' (= $CS km)"
echo "$SCRIPTNAME       fine grid = '$FINEGRID' (= $FS km)"

# cat prefix and exec together
PISM="${PISM_PREFIX}${PISM_EXEC} -ocean_kill"

# coupler settings for pre-spinup
COUPLER_SIMPLE="-atmosphere searise_greenland -surface pdd -pdd_fausto"
# coupler settings for spin-up (i.e. with forcing)
COUPLER_FORCING="-atmosphere searise_greenland,forcing -surface pdd -pdd_fausto -dTforcing $PISM_TEMPSERIES -dSLforcing $PISM_SLSERIES"

# default choices in parameter study; see Bueler & Brown (2009) re "tillphi"
TILLPHI="-topg_to_phi 5.0,20.0,-300.0,700.0,10.0"

# use "best run" parameters from Bueler et al. (GRL)
GRLPARAMS="-e 1 -pseudo_plastic_q 0.10 -plastic_pwfrac 0.99"

FULLPHYS="-ssa -super -plastic -thk_eff ${GRLPARAMS}"

echo "$SCRIPTNAME      executable = '$PISM'"
echo "$SCRIPTNAME         tillphi = '$TILLPHI'"
echo "$SCRIPTNAME    full physics = '$FULLPHYS'"
echo "$SCRIPTNAME  simple coupler = '$COUPLER_SIMPLE'"
echo "$SCRIPTNAME forcing coupler = '$COUPLER_FORCING'"



# bootstrap and do smoothing run to 100 years
PRE0NAME=g${CS}km_pre100.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run (for ${SMOOTHRUNLENGTH}a)"
cmd="$PISM_MPIDO $NN $PISM -skip $COARSESKIP -boot_from $INNAME $COARSEGRID \
  $COUPLER_SIMPLE -y ${SMOOTHRUNLENGTH} -o $PRE0NAME"
$PISM_DO $cmd


# run with -no_mass (no surface change) for 50ka
PRE1NAME=g${CS}km_steady.nc
EX1NAME=ex_${PRE1NAME}
EXTIMES=0:500:${NOMASSSIARUNLENGTH}
EXVARS="enthalpybase,temppabase,mask" # add mask, so that check_stationarity.py ignores ice-free areas.
echo
echo "$SCRIPTNAME  -no_mass (no surface change) SIA run to achieve approximate temperature equilibrium, for ${NOMASSSIARUNLENGTH}a"
cmd="$PISM_MPIDO $NN $PISM -skip $COARSESKIP -i $PRE0NAME $COUPLER_SIMPLE \
  -no_mass -y ${NOMASSSIARUNLENGTH} \
  -extra_file $EX1NAME -extra_vars $EXVARS -extra_times $EXTIMES -o $PRE1NAME"
$PISM_DO $cmd

# smoothing for 100 years
PRE2NAME=g${CS}km_SIA.nc
echo
echo "$SCRIPTNAME  smoothing with SIA for ${SMOOTHRUNLENGTH}a"
cmd="$PISM_MPIDO $NN $PISM -skip $COARSESKIP -i $PRE1NAME $COUPLER_SIMPLE -y $SMOOTHRUNLENGTH -o $PRE2NAME"
$PISM_DO $cmd



# pre-spinup done; ready to use paleoclimate forcing ...
# spinup is split into 5 substages

# start back at 125ka BPE, use GRIP and SPECMAP forcing, run till SPLIT1TIME
ENDTIME=$COARSEENDTIME
OUTNAME=g${CS}km_m40ka.nc
SNAPSNAME=snaps_g${CS}km_m40ka.nc
TSNAME=ts_g${CS}km_m40ka.nc
echo
echo "$SCRIPTNAME  paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      except bed deformation, from $PALEOSTARTYEAR a to ${ENDTIME}a"
cmd="$PISM_MPIDO $NN $PISM -skip $COARSESKIP -i $PRE2NAME $TILLPHI $FULLPHYS $COUPLER_FORCING \
     -ts_file $TSNAME -ts_times $PALEOSTARTYEAR:1:$ENDTIME \
     -save_file $SNAPSNAME -save_times -120000:5000:-45000 \
     -ys $PALEOSTARTYEAR -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd

# uncomment to stop here:
#exit

STARTTIME=$ENDTIME
ENDTIME=-30000 # BP
STARTNAME=$OUTNAME
OUTNAME=g${FS}km_m30ka.nc
SNAPSNAME=g${FS}km_m35ka.nc
SNAPSTIME=-35000
TSNAME=ts_g${FS}km_m30ka.nc
echo
echo "$SCRIPTNAME  regrid and do paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation, from ${STARTTIME}a BPE to ${ENDTIME}a BPE"
cmd="$PISM_MPIDO $NN $PISM -skip $FINESKIP -boot_from $INNAME $FINEGRID $TILLPHI $FULLPHYS -bed_def_lc $COUPLER_FORCING\
     -regrid_file $STARTNAME -regrid_vars litho_temp,thk,enthalpy,bwat  \
     -ts_file $TSNAME -ts_times $STARTTIME:1:$ENDTIME \
      -save_file $SNAPSNAME -save_times $SNAPSTIME \
     -ys $STARTTIME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd

# to stop here:
#exit

pism5kmopts="$PISM_MPIDO $NN $PISM -skip $FINESKIP $FULLPHYS -bed_def_lc $COUPLER_FORCING"

STARTTIME=$ENDTIME
ENDTIME=-20000 # BP
STARTNAME=$OUTNAME
OUTNAME=g${FS}km_m20ka.nc
SNAPSNAME=g${FS}km_m25ka.nc
SNAPSTIME=-25000
TSNAME=ts_g${FG}km_m20ka.nc
echo
echo "$SCRIPTNAME  paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation, from ${STARTTIME}a BPE to ${ENDTIME}a BPE"
cmd="$pism5kmopts -ts_file $TSNAME -ts_times $STARTTIME:1:$ENDTIME \
     -save_file $SNAPSNAME -save_times $SNAPSTIME -i $STARTNAME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd

STARTTIME=$ENDTIME
ENDTIME=-10000 # BP
STARTNAME=$OUTNAME
OUTNAME=g${FS}km_m10ka.nc
SNAPSNAME=g${FS}km_m15ka.nc
SNAPSTIME=-15000
TSNAME=ts_g${FG}km_m10ka.nc
echo
echo "$SCRIPTNAME  paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation, from ${STARTTIME}a BPE to ${ENDTIME}a BPE"
cmd="$pism5kmopts -ts_file $TSNAME -ts_times $STARTTIME:1:$ENDTIME \
     -save_file $SNAPSNAME -save_times $SNAPSTIME -i $STARTNAME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd

STARTTIME=$ENDTIME
ENDTIME=0 # BP
STARTNAME=$OUTNAME
OUTNAME=g${FS}km_0.nc
SNAPSNAME=g${FS}km_m5ka.nc
SNAPSTIME=-5000
TSNAME=ts_g${FS}km_0.nc
echo
echo "$SCRIPTNAME  paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation, from ${STARTTIME}a BPE to ${ENDTIME}a BPE (present)"
cmd="$pism5kmopts -ts_file $TSNAME -ts_times $STARTTIME:1:$ENDTIME \
     -save_file $SNAPSNAME -save_times $SNAPSTIME -i $STARTNAME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd


echo
echo "$SCRIPTNAME  spinup done"
