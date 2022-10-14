#!/bin/bash

# Copyright (C) 2013, 2014, 2015 Andy Aschwanden
#
# *****************************************************************************
# Relax Greenland Topography
# *****************************************************************************
#
#

# recommended way to run with N processors is " ./run-relax.sh N >& out.relax & "
# which gives a viewable (with "less", for example) transcript in out.relax

SCRIPTNAME="#(run-relax.sh)"

echo
echo "# =================================================================================="
echo "# PISM run for relaxation experiment "
echo "# =================================================================================="
echo

set -e  # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "psearise.sh 8" then NN = 8
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

# set output format:
#  $ export PISM_OFORMAT="netcdf4_parallel "
if [ -n "${PISM_OFORMAT:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME                      PISM_OFORMAT = $PISM_OFORMAT  (already set)"
else
  PISM_OFORMAT="netcdf3"
  echo "$SCRIPTNAME                      PISM_OFORMAT = $PISM_OFORMAT"
fi
OFORMAT=$PISM_OFORMAT

echo

# preprocess.sh generates pism_*.nc files; run it first
if [ -n "${PISM_DATANAME:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME   PISM_DATANAME = $PISM_DATANAME  (already set)"
else
  PISM_DATANAME=pism_Greenland_5km_v1.1.nc
fi

if [ -n "${PISM_TARGETNAME:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME   PISM_TARGETNAME = $PISM_TARGETNAME  (already set)"
else
  PISM_TARGETNAME=target_$PISM_DATANAME
fi

for INPUT in $PISM_DATANAME $PISM_TARGETNAME; do
  if [ -e "$INPUT" ] ; then  # check if file exist
    echo "$SCRIPTNAME           input   $INPUT (found)"
  else
    echo "$SCRIPTNAME           input   $INPUT (MISSING!!)"
    echo
    echo "     !!!!   RUN  preprocess.sh  TO GENERATE  $INPUT   !!!!"
    echo
  fi
done

# grids
GRID20KM="-Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 21 -z_spacing equal"
GRID10KM="-Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 201 -Mbz 21 -z_spacing equal"
GRID5KM="-Mx 301 -My 561 -Lz 4000 -Lbz 2000 -Mz 201 -Mbz 21 -z_spacing equal"
GRID3KM="-Mx 501 -My 935 -Lz 4000 -Lbz 2000 -Mz 401 -Mbz 21 -z_spacing equal"
GRID2FKM="-Mx 601 -My 1121 -Lz 4000 -Lbz 2000 -Mz 401 -Mbz 21 -z_spacing equal"
GRID2KM="-Mx 751 -My 1401 -Lz 4000 -Lbz 2000 -Mz 401 -Mbz 21 -z_spacing equal"
GRID1KM="-Mx 1501 -My 2801 -Lz 4000 -Lbz 2000 -Mz 401 -Mbz 21 -z_spacing equal "


# grid spacings
GS20KM=20
GS10KM=10
GS5KM=5
GS3KM=3
GS2FKM=2.5
GS2KM=2
GS1KM=1

# skips
SKIP20KM=10
SKIP10KM=50
SKIP5KM=200
SKIP3KM=500
SKIP2FKM=750
SKIP2KM=1000
SKIP1KM=2000

# defaults to coarse grid choices
GRID=$GRID20KM
SKIP=$SKIP20KM
GS=$GS20KM

echo ""
if [ $# -gt 1 ] ; then
  if [ $2 -eq "1" ] ; then  # if user says "run-relax.sh N 1" then use:
    echo "$SCRIPTNAME grid: RUNS ON 20km and 10km"
    echo "$SCRIPTNAME       WARNING: MEDIUM COMPUTATIONAL TIME"
    GRID=$GRID10KM
    SKIP=$SKIP10KM
    GS=$GS10KM
  fi
  if [ $2 -eq "2" ] ; then  # if user says "run-relax.sh N 2" then use:
    echo "$SCRIPTNAME grid: RUNS ON 20km, 10km, 5km and 5km"
    echo "$SCRIPTNAME       WARNING: LARGE COMPUTATIONAL TIME"
    GRID=$GRID5KM
    SKIP=$SKIP5KM
    GS=$GS5KM
  fi
  if [ $2 -eq "3" ] ; then  # if user says "run-relax.sh N 3":
    echo "$SCRIPTNAME grid: RUNS ON 20km, 10km, 5km and 2.5km"
    echo "$SCRIPTNAME       WARNING: VERY LARGE COMPUTATIONAL TIME"
    GRID=$GRID2FKM
    SKIP=$SKIP2FKM
    GS=$GS2FKM
  fi
  if [ $2 -eq "4" ] ; then  # if user says "run-relax.sh N 4":
    echo "$SCRIPTNAME grid: RUNS ON 10km, 5km, 2.5km and 2km"
    echo "$SCRIPTNAME       WARNING: EXTREMELY LARGE COMPUTATIONAL TIME"
    GRID=$GRID2KM
    SKIP=$SKIP2KM
    GS=$GS2KM
  fi
  if [ $2 -eq "5" ] ; then  # if user says "run-relax.sh N 5":
    echo "$SCRIPTNAME grid: RUNS ON 10km, 5km, 2.5km, 2km and 1km"
    echo "$SCRIPTNAME       WARNING: EXTREMELY LARGE COMPUTATIONAL TIME"
    GRID=$GRID1KM
    SKIP=$SKIP1KM
    GS=$GS1KM
  fi
else
    echo "$SCRIPTNAME grid: ALL RUNS ON 20km"
fi
echo ""

OCEAN="-calving ocean_kill"
# FIXME (CK): OCEAN (above) is never used
COUPLER="-surface given -surface_given_file $PISM_TARGETNAME"
PISM="${PISM_BIN}${PISM_EXEC} -energy none -bed_def lc"
# output file size
OSIZE="big"

BOOTNAME=$PISM_DATANAME

STARTTIME=0
ENDTIME=100

OUTNAME=g${GS}km_removeice.nc
TSNAME=ts_$OUTNAME
TSSTEP=yearly
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXSTEP=100
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME
EXVARS="diffusivity,thk,mask,lat,lon,taud_mag,topg,usurf"


cmd="$PISM_MPIDO $NN $PISM -skip -skip_max $SKIP -i $BOOTNAME -bootstrap 
     $GRID $COUPLER \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -spatial_file $EXNAME -spatial_vars $EXVARS -spatial_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o_size $OSIZE -o $OUTNAME -o_format $OFORMAT"

$PISM_DO $cmd


COUPLER="-no_mass"

RUNLENGTH=50000
STARTTIME=$ENDTIME
ENDTIME=$(($STARTTIME + $RUNLENGTH))

INNAME=$OUTNAME
OUTNAME=g${GS}km_relax.nc
TSNAME=ts_$OUTNAME
TSSTEP=yearly
TSTIMES=$STARTTIME:$TSSTEP:$ENDTIME
EXNAME=ex_$OUTNAME
EXSTEP=100
EXTIMES=$STARTTIME:$EXSTEP:$ENDTIME
EXVARS="diffusivity,thk,mask,lat,lon,taud_mag,topg,usurf"

echo
cmd="$PISM_MPIDO $NN $PISM -skip -skip_max $SKIP -i $INNAME 
     $COUPLER \
     -ts_file $TSNAME -ts_times $TSTIMES \
     -spatial_file $EXNAME -spatial_vars $EXVARS -spatial_times $EXTIMES \
     -ys $STARTTIME -ye $ENDTIME -o_size $OSIZE -o $OUTNAME -o_format $OFORMAT"

$PISM_DO $cmd
