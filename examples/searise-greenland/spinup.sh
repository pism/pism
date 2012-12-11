#!/bin/bash

# Copyright (C) 2009-2012 The PISM Authors

# PISM SeaRISE Greenland spinup using modeled paleoclimate
#
# before using this script, run preprocess.sh to download and adjust metadata
# on SeaRISE "Present Day Greenland" master dataset
#
# recommended way to run with N processors is " ./spinup.sh N >& out.spinup & "
# which gives a viewable (with "less", for example) transcript in out.spinup

# seconds per year, from UDUNITS
SECPERA=3.15569259747e7

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
if [ $# -gt 0 ] ; then  # if user says "spinup.sh 8" then NN = 8
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

# preprocess.sh generates pism_*.nc files; run it first
if [ -n "${PISM_DATANAME:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME   PISM_DATANAME = $PISM_DATANAME  (already set)"
else
  PISM_DATAVERSION=1.1
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

PISM_CONFIG=searise_config.nc

for INPUT in $PISM_DATANAME $PISM_TEMPSERIES $PISM_SLSERIES $PISM_CONFIG; do
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

# run lengths and starting time for paleo
SMOOTHRUNLENGTH=100
NOMASSSIARUNLENGTH=50000
PALEOSTARTYEAR=-125000
FTTENDTIME=-100

# grids
VDIMS="-Lz 4000 -Lbz 2000"
COARSEVGRID="${VDIMS} -Mz 101 -Mbz 11 -z_spacing equal"
FINEVGRID="${VDIMS} -Mz 201 -Mbz 21 -z_spacing equal"
FINESTVGRID="${VDIMS} -Mz 401 -Mbz 41 -z_spacing equal"

TWENTYKMGRID="-Mx 76 -My 141 ${COARSEVGRID}"
TENKMGRID="-Mx 151 -My 281 ${FINEVGRID}"
FIVEKMGRID="-Mx 301 -My 561 ${FINESTVGRID}"

# skips
SKIPTWENTYKM=10
SKIPTENKM=50
SKIPFIVEKM=200

# defaults to coarse grid choices
COARSEGRID=$TWENTYKMGRID
FINEGRID=$TWENTYKMGRID
COARSESKIP=$SKIPTWENTYKM
FINESKIP=$SKIPTWENTYKM
CS=20 # km
FS=20 # km

COARSEENDTIME=-5000 # BP

echo ""
if [ $# -gt 1 ] ; then
  if [ $2 -eq "1" ] ; then  # if user says "spinup.sh N 1" then MEDIUM:
    echo "$SCRIPTNAME grid: ALL RUNS ON 10km"
    echo "$SCRIPTNAME       WARNING: LARGE COMPUTATIONAL TIME"
    COARSEGRID=$TENKMGRID
    FINEGRID=$TENKMGRID
    COARSESKIP=$SKIPTENKM
    FINESKIP=$SKIPTENKM
    CS=10 # km
    FS=10 # km
  fi
  if [ $2 -eq "2" ] ; then  # if user says "spinup.sh N 2" then FINE:
    echo "$SCRIPTNAME grid: ON 10km TILL ${COARSEENDTIME}a THEN ON 5km"
    echo "$SCRIPTNAME       WARNING: VERY LARGE COMPUTATIONAL TIME"
    COARSEGRID=$TENKMGRID
    FINEGRID=$FIVEKMGRID
    COARSESKIP=$SKIPTENKM
    FINESKIP=$SKIPFIVEKM
    CS=10 # km
    FS=5  # km
  fi
else
    echo "$SCRIPTNAME grid: ALL RUNS ON 20km"
fi
echo ""

echo "$SCRIPTNAME     coarse grid = '$COARSEGRID' (= $CS km), with -skip = $COARSESKIP)"
echo "$SCRIPTNAME       fine grid = '$FINEGRID' (= $FS km), with -skip = $FINESKIP)"

# cat prefix and exec together
PISM="${PISM_PREFIX}${PISM_EXEC} -config_override $PISM_CONFIG -climatic_mass_balance_cumulative"

# coupler settings for pre-spinup
COUPLER_SIMPLE="-atmosphere searise_greenland -surface pdd -ocean_kill $INNAME"
# coupler settings for spin-up (i.e. with forcing)
COUPLER_FORCING="-atmosphere searise_greenland,delta_T -surface pdd -paleo_precip $PISM_TEMPSERIES -atmosphere_delta_T_file $PISM_TEMPSERIES -ocean constant,delta_SL -ocean_delta_SL_file $PISM_SLSERIES -ocean_kill $INNAME"
# coupler settings for spin-up (i.e. with forcing) and force-to-thickness
COUPLER_FTT="-atmosphere searise_greenland,delta_T -surface pdd,forcing -paleo_precip $PISM_TEMPSERIES -atmosphere_delta_T_file $PISM_TEMPSERIES -ocean constant,delta_SL -ocean_delta_SL_file $PISM_SLSERIES -ocean_kill $INNAME"

# default choices in parameter study; see Bueler & Brown (2009) re "tillphi"
TILLPHI="-topg_to_phi 5.0,20.0,-300.0,700.0"

FULLPHYS="-ssa_sliding ${TILLPHI}"

echo "$SCRIPTNAME      executable = '$PISM'"
echo "$SCRIPTNAME    full physics = '$FULLPHYS'"
echo "$SCRIPTNAME  simple coupler = '$COUPLER_SIMPLE'"
echo "$SCRIPTNAME forcing coupler = '$COUPLER_FORCING'"


# bootstrap and do smoothing run to 100 years
PRE0NAME=g${CS}km_pre100.nc
echo
echo "$SCRIPTNAME  bootstrapping plus short smoothing run (for ${SMOOTHRUNLENGTH}a)"
cmd="$PISM_MPIDO $NN $PISM -skip -skip_max  $COARSESKIP -boot_file $INNAME $COARSEGRID \
  $COUPLER_SIMPLE -y ${SMOOTHRUNLENGTH} -o $PRE0NAME"
$PISM_DO $cmd


# quick look at climate in 500 a recent period; see also delta_T in pism_dT.nc
CLIMSTARTTIME=-500
PRE0CLIMATE=g${CS}km_climate${CLIMSTARTTIME}a.nc
PCLIM="${PISM_PREFIX}pclimate"
echo
echo "$SCRIPTNAME  running pclimate to show climate in modern period [${CLIMSTARTTIME} a,0 a], using current geometry and 10 year subintervals"
cmd="$PISM_MPIDO $NN $PCLIM -i $PRE0NAME $COUPLER_FORCING -times $CLIMSTARTTIME:10:0 -o $PRE0CLIMATE"
$PISM_DO $cmd


# run with -no_mass (no surface change) for 50ka
PRE1NAME=g${CS}km_steady.nc
EX1NAME=ex_${PRE1NAME}
EXTIMES=0:500:${NOMASSSIARUNLENGTH}
EXVARS="diffusivity,temppabase,tempicethk_basal,bmelt,bwp,csurf,hardav,mask" # check_stationarity.py can be applied to ex_${PRE1NAME}
echo
echo "$SCRIPTNAME  -no_mass (no surface change) SIA run to achieve approximate temperature equilibrium, for ${NOMASSSIARUNLENGTH}a"
cmd="$PISM_MPIDO $NN $PISM -skip -skip_max  $COARSESKIP -i $PRE0NAME $COUPLER_SIMPLE \
  -no_mass -ys 0 -y ${NOMASSSIARUNLENGTH} \
  -extra_file $EX1NAME -extra_vars $EXVARS -extra_times $EXTIMES -o $PRE1NAME"
$PISM_DO $cmd


# pre-spinup done; ready to use paleoclimate forcing for real spinup ...

EXSTEP=500
EXFSTEP=10
TSSTEP=yearly

STARTTIME=$PALEOSTARTYEAR

ENDTIME=$COARSEENDTIME

ET=$(($ENDTIME/-1))
OUTNAME=g${CS}km_m${ET}a.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME
EXVARS="diffusivity,temppabase,tempicethk_basal,bmelt,bwp,csurf,hardav,mask,dHdt,cbase,tauc,thk,topg,usurf,climatic_mass_balance_cumulative"
echo
echo "$SCRIPTNAME  paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation, from $PALEOSTARTYEAR a to ${ENDTIME}a"
cmd="$PISM_MPIDO $NN $PISM -skip -skip_max  $COARSESKIP -i $PRE1NAME $FULLPHYS -bed_def lc $COUPLER_FORCING \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd

# exit    # uncomment to stop here

STARTTIME=$ENDTIME
ENDTIME=0 # BP
STARTNAME=$OUTNAME

# ######################################
# "regular" run
# ######################################

OUTNAME=g${FS}km_0.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME
echo
echo "$SCRIPTNAME  regular run"
echo "$SCRIPTNAME  regrid to fine grid and do paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation,"
echo "$SCRIPTNAME      from ${STARTTIME}a BPE to ${ENDTIME}a BPE"
cmd="$PISM_MPIDO $NN $PISM -skip -skip_max  $FINESKIP -boot_file $INNAME $FINEGRID $FULLPHYS \
     -bed_def lc $COUPLER_FORCING \
     -regrid_file $STARTNAME -regrid_vars litho_temp,thk,enthalpy,bwat,bmelt -regrid_bed_special  \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd


# ######################################
# "force-to-thickness" run
# ######################################

STARTTIME=$COARSEENDTIME
ENDTIME=$FTTENDTIME
ET=$(($FTTENDTIME/-1))
OUTNAME=g${FS}km_m${ET}a_ftt.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$FTTENDTIME
EXNAME=ex_$OUTNAME
EXTIMES=$STARTTIME:$EXSTEP:$FTTENDTIME
echo
echo "$SCRIPTNAME  force-to-thickness run"
echo "$SCRIPTNAME  regrid to fine grid and do paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation and modified surface mass balance,"
echo "$SCRIPTNAME      from ${STARTTIME}a BPE to ${ENDTIME}a BPE"
cmd="$PISM_MPIDO $NN $PISM -skip -skip_max  $FINESKIP -boot_file $INNAME $FINEGRID $FULLPHYS \
     -bed_def lc $COUPLER_FTT \
     -force_to_thk $INNAME -force_to_thk_alpha 0.005 \
     -regrid_file $STARTNAME -regrid_vars litho_temp,thk,enthalpy,bwat -regrid_bed_special  \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $FTTENDTIME -o $OUTNAME"
$PISM_DO $cmd

STARTTIME=$FTTENDTIME
ENDTIME=0
STARTNAME=$OUTNAME
OUTNAME=g${FS}km_0_ftt.nc
TSNAME=ts_$OUTNAME
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXTIMES=$STARTTIME:$EXFSTEP:$ENDTIME
echo
echo "$SCRIPTNAME  force-to-thickness run finishes with frequent diagnostic information"
echo "$SCRIPTNAME  do paleo-climate forcing run with full physics,"
echo "$SCRIPTNAME      including bed deformation and modified surface mass balance,"
echo "$SCRIPTNAME      from ${STARTTIME}a BPE to ${ENDTIME}a BPE"
cmd="$PISM_MPIDO $NN $PISM -skip -skip_max  $FINESKIP -i $STARTNAME $FULLPHYS \
     -bed_def lc $COUPLER_FTT \
     -force_to_thk $INNAME -force_to_thk_alpha 0.005 \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -extra_file $EXNAME -extra_vars $EXVARS -extra_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o $OUTNAME"
$PISM_DO $cmd

echo
echo "$SCRIPTNAME  some postprocessing"
echo
# calculate yearly-averages of climatic_mass_balance and dHdt using ncap2 sleight of hand.
cmd="ncap2 -O -s '*sz_idt=time.size(); climatic_mass_balance[\$time,\$x,\$y]= 0.f; dHdt[\$time,\$x,\$y]= 0.f; for(*idt=1 ; idt<sz_idt ; idt++) {climatic_mass_balance(idt,:,:)=(climatic_mass_balance_cumulative(idt,:,:)-climatic_mass_balance_cumulative(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA; dHdt(idt,:,:)=(thk(idt,:,:)-thk(idt-1,:,:))/(time(idt)-time(idt-1))*$SECPERA;}' $EXNAME $EXNAME"
$PISM_DO $cmd

echo
# adjust meta data for new fields
cmd="ncatted -a units,climatic_mass_balance,o,c,'m year-1' -a units,dHdt,o,c,'m year-1' \
      -a long_name,climatic_mass_balance,o,c,'surface mass balance' \
      -a long_name,dHdt,o,c,'rate of change of ice thickness' \
      -a grid_mapping,climatic_mass_balance,o,c,'mapping' \
      -a grid_mapping,dHdt,o,c,'mapping' \
      -a cell_methods,climatic_mass_balance,o,c,'time: mean (interval: $EXFSTEP years)' \
      -a cell_methods,dHdt,o,c,'time: mean (interval: $EXFSTEP years)' $EXNAME"
$PISM_DO $cmd
echo
# now extract last climatic_mass_balance record
cmd="ncks -A -v climatic_mass_balance -d time,$ENDTIME. $EXNAME $OUTNAME"
$PISM_DO $cmd

echo
echo "$SCRIPTNAME  spinup done"

