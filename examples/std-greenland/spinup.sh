#!/bin/bash

# Copyright (C) 2009-2015, 2017, 2018 The PISM Authors

# PISM Greenland spinup using either constant present-day climate or modeled
# paleoclimate.  See README.md.

# Before using this script, run preprocess.sh to download and adjust metadata
# in the SeaRISE "Present Day Greenland" master dataset.

set -e  # exit on error

GRIDLIST="{40, 20, 10, 5, 3, 2}"
CLIMLIST="{const, paleo}"
DYNALIST="{sia, hybrid}"

# preprocess.sh generates pism_*.nc files; run it first
PISM_DATAVERSION=1.1
PISM_DATANAME=pism_Greenland_5km_v$PISM_DATAVERSION.nc
PISM_TEMPSERIES=pism_dT.nc
PISM_SLSERIES=pism_dSL.nc

if [ $# -lt 5 ] ; then
  cat - <<EOF
spinup.sh ERROR: needs 5 or 6 or 7 positional arguments ... ENDING NOW

usage:

    spinup.sh PROCS CLIMATE DURATION GRID DYNAMICS [OUTFILE] [BOOTFILE]

  where:
    PROCS     = 1,2,3,... is number of MPI processes
    CLIMATE   in $CLIMLIST
    DURATION  = model run time in years; does '-ys -DURATION -ye 0'
    GRID      in $GRIDLIST (km)
    DYNAMICS  in $DYNALIST; sia is non-sliding; default = sia
    OUTFILE   optional name of output file; default = unnamed.nc
    BOOTFILE  optional name of input file; default = $PISM_DATANAME

consider setting optional environment variables (see script for meaning):
    EXSTEP       spacing in years between -extra_files outputs; defaults to 100
    EXVARS       desired -extra_vars; defaults to 'diffusivity,temppabase,
                   tempicethk_basal,bmelt,tillwat,velsurf_mag,mask,thk,topg,usurf'
                   plus ',hardav,velbase_mag,tauc' if DYNAMICS=hybrid
    NODIAGS      if set, DON'T use -ts_file or -extra_file
    USEPIK       if set, add -pik -subgl
    PARAM_PPQ    sets (hybrid-only) option -pseudo_plastic_q \$PARAM_PPQ
                   [default=0.25]
    PARAM_SIAE   sets option -sia_e \$PARAM_SIAE   [default=3.0]
    PARAM_TEFO   sets (hybrid-only) option -till_effective_fraction_overburden
                   \$PARAM_TEFO   [default=0.02]
    PARAM_TTPHI  sets (hybrid-only) option -topg_to_phi \$PARAM_TTPHI
                   [default=15.0,40.0,-300.0,700.0]
    PARAM_NOSGL  if set, DON'T use -tauc_slippery_grounding_lines
    PISM_DO      set to 'echo' if no run desired; defaults to empty
    PISM_MPIDO   defaults to 'mpiexec -n'
    PISM_BIN  set to path to pismr executable if desired; defaults to empty
    PISM_EXEC    defaults to 'pismr'
    REGRIDFILE   set to file name to regrid from; defaults to empty (no regrid)
    REGRIDVARS   desired -regrid_vars; applies *if* REGRIDFILE set;
                   defaults to 'basal_melt_rate_grounded,enthalpy,litho_temp,thk,tillwat'

example usage 1:

    $ ./spinup.sh 4 const 1000 20 sia

  Does spinup with 4 processors, constant-climate, 1000 year run, 20 km
  grid, and non-sliding SIA stress balance.  Bootstraps from and outputs to
  default files.

example usage 2:

    $ PISM_DO=echo ./spinup.sh 128 paleo 100.0 5 hybrid out.nc boot.nc &> foo.sh

  Creates a script foo.sh for spinup with 128 processors, simulated paleo-climate,
  5 km grid, sliding with SIA+SSA hybrid, output to {out.nc,ts_out.nc,ex_out.nc},
  and bootstrapping from boot.nc.
EOF
  exit
fi

if [ -n "${SCRIPTNAME:+1}" ] ; then
  echo "[SCRIPTNAME=$SCRIPTNAME (already set)]"
  echo ""
else
  SCRIPTNAME="#(spinup.sh)"
fi

if [ $# -gt 7 ] ; then
  echo "$SCRIPTNAME WARNING: ignoring arguments after argument 7 ..."
fi

NN="$1" # first arg is number of processes
DURATION=$3
RUNSTARTEND="-ys -$DURATION -ye 0"

# set coupler from argument 2
if [ "$2" = "const" ]; then
  climname="constant-climate"
  INLIST=""
  COUPLER="-surface given -surface_given_file $PISM_DATANAME"
elif [ "$2" = "paleo" ]; then
  climname="paleo-climate"
  INLIST="$PISM_TEMPSERIES $PISM_SLSERIES"
  COUPLER=" -bed_def lc -atmosphere searise_greenland,delta_T,paleo_precip -surface pdd -atmosphere_paleo_precip_file $PISM_TEMPSERIES -atmosphere_delta_T_file $PISM_TEMPSERIES -sea_level constant,delta_sl -ocean_delta_sl_file $PISM_SLSERIES"
else
  echo "invalid second argument; must be in $CLIMLIST"
  exit
fi

# decide on grid and skip from argument 4
COARSESKIP=10
FINESKIP=20
FINESTSKIP=50
VDIMS="-Lz 4000 -Lbz 2000 -skip -skip_max "
COARSEVGRID="-Mz 101 -Mbz 11 -z_spacing equal ${VDIMS} ${COARSESKIP}"
FINEVGRID="-Mz 201 -Mbz 21 -z_spacing equal ${VDIMS} ${FINESKIP}"
FINESTVGRID="-Mz 401 -Mbz 41 -z_spacing equal ${VDIMS} ${FINESTSKIP}"
if [ "$4" -eq "40" ]; then
  dx=40
  myMx=38
  myMy=71
  vgrid=$COARSEVGRID
elif [ "$4" -eq "20" ]; then
  dx=20
  myMx=76
  myMy=141
  vgrid=$COARSEVGRID
elif [ "$4" -eq "10" ]; then
  dx=10
  myMx=151
  myMy=281
  vgrid=$FINEVGRID
elif [ "$4" -eq "5" ]; then
  # "native" resolution in data file, with 561 x 301 grid
  dx=5
  myMx=301
  myMy=561
  vgrid=$FINEVGRID
elif [ "$4" -eq "3" ]; then
  dx=3
  myMx=501
  myMy=934
  vgrid=$FINEVGRID
elif [ "$4" -eq "2" ]; then
  dx=2
  myMx=750
  myMy=1400
  vgrid=$FINESTVGRID
else
  echo "invalid fourth argument: must be in $GRIDLIST"
  exit
fi

grid="-Mx $myMx -My $myMy $vgrid -grid.recompute_longitude_and_latitude false -grid.registration corner"

# set stress balance from argument 5
if [ -n "${PARAM_SIAE:+1}" ] ; then  # check if env var is already set
  PHYS="-calving ocean_kill -ocean_kill_file ${PISM_DATANAME} -sia_e ${PARAM_SIAE}"
else
  PHYS="-calving ocean_kill -ocean_kill_file ${PISM_DATANAME} -sia_e 3.0"
fi
if [ -n "${USEPIK:+1}" ] ; then  # check if env var is already set
  PHYS="${PHYS} -pik -subgl"
fi

# done forming $PHYS if "$5" = "sia"
if [ "$5" = "hybrid" ]; then
  if [ -z "${PARAM_TTPHI}" ] ; then  # check if env var is NOT set
    PARAM_TTPHI="15.0,40.0,-300.0,700.0"
  fi
  if [ -z "${PARAM_PPQ}" ] ; then  # check if env var is NOT set
    PARAM_PPQ="0.25"
  fi
  if [ -z "${PARAM_TEFO}" ] ; then  # check if env var is NOT set
    PARAM_TEFO="0.02"
  fi
  if [ -z "${PARAM_NOSGL}" ] ; then  # check if env var is NOT set
    SGL="-tauc_slippery_grounding_lines"
  else
    SGL=""
  fi
  PHYS="${PHYS} -stress_balance ssa+sia -topg_to_phi ${PARAM_TTPHI} -pseudo_plastic -pseudo_plastic_q ${PARAM_PPQ} -till_effective_fraction_overburden ${PARAM_TEFO} ${SGL}"
else
  if [ "$5" = "sia" ]; then
    echo "$SCRIPTNAME  sia-only case: ignoring PARAM_TTPHI, PARAM_PPQ, PARAM_TEFO ..."
  else
    echo "invalid fifth argument; must be in $DYNALIST"
    exit
  fi
fi

# set output filename from argument 6
if [ -z "$6" ]; then
  OUTNAME=unnamed.nc
else
  OUTNAME=$6
fi

# set bootstrapping input filename from argument 6
if [ -z "$7" ]; then
  INNAME=$PISM_DATANAME
else
  INNAME=$7
fi
INLIST="${INLIST} $INNAME $REGRIDFILE"

# now we have read options ... we know enough to report to user ...
echo
echo "# ======================================================================="
echo "# PISM std Greenland spinup:"
echo "#    $NN processors, $DURATION a run, $dx km grid, $climname, $5 dynamics"
echo "# ======================================================================="

# actually check for input files
for INPUT in $INLIST; do
  if [ -e "$INPUT" ] ; then  # check if file exist
    echo "$SCRIPTNAME           input   $INPUT (found)"
  else
    echo "$SCRIPTNAME           input   $INPUT (MISSING!!)"
    echo
    echo "$SCRIPTNAME  ***WARNING***  you may need to run ./preprocess.sh to generate standard input files!"
    echo
  fi
done

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

# set EXSTEP to default if not set
if [ -n "${EXSTEP:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME          EXSTEP = $EXSTEP  (already set)"
else
  EXSTEP="100"
  echo "$SCRIPTNAME          EXSTEP = $EXSTEP"
fi

# set EXVARS list to defaults if not set
if [ -n "${EXVARS:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME          EXVARS = $EXVARS  (already set)"
else
  EXVARS="diffusivity,temppabase,tempicethk_basal,bmelt,tillwat,velsurf_mag,mask,thk,topg,usurf"
  if [ "$5" = "hybrid" ]; then
    EXVARS="${EXVARS},hardav,velbase_mag,tauc"
  fi
  echo "$SCRIPTNAME          EXVARS = $EXVARS"
fi

# if REGRIDFILE set then form regridcommand
if [ -n "${REGRIDFILE:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME      REGRIDFILE = $REGRIDFILE"
  if [ -n "${REGRIDVARS:+1}" ] ; then  # check if env var is already set
    echo "$SCRIPTNAME      REGRIDVARS = $REGRIDVARS  (already set)"
  else
    REGRIDVARS='litho_temp,thk,enthalpy,tillwat,basal_melt_rate_grounded'
    # note: other vars which are "state":  Href, dbdt, shelfbtemp, shelfbmassflux
    echo "$SCRIPTNAME      REGRIDVARS = $REGRIDVARS"
  fi
  regridcommand="-regrid_file $REGRIDFILE -regrid_vars $REGRIDVARS"
else
  regridcommand=""
fi

# show remaining setup options:
PISM="${PISM_BIN}${PISM_EXEC}"
echo "$SCRIPTNAME      executable = '$PISM'"
echo "$SCRIPTNAME         coupler = '$COUPLER'"
echo "$SCRIPTNAME        dynamics = '$PHYS'"

# set up diagnostics
if [ -z "${NODIAGS}" ] ; then  # check if env var is NOT set
  TSNAME=ts_$OUTNAME
  TSTIMES=-$DURATION:yearly:0
  EXNAME=ex_$OUTNAME
  EXTIMES=-$DURATION:$EXSTEP:0
  # check_stationarity.py can be applied to $EXNAME
  DIAGNOSTICS="-ts_file $TSNAME -ts_times $TSTIMES -extra_file $EXNAME -extra_times $EXTIMES -extra_vars $EXVARS"
else
  DIAGNOSTICS=""
fi

# construct command
cmd="$PISM_MPIDO $NN $PISM -i $INNAME -bootstrap ${grid} $RUNSTARTEND $regridcommand $COUPLER $PHYS $DIAGNOSTICS -o $OUTNAME"
echo
$PISM_DO $cmd

