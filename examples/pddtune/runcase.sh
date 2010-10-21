#!/bin/bash

# Copyright (C) 2010 Ed Bueler

# see README for role of this script

SCRIPTNAME="#(runcase.sh)"

set -e  # exit on error

NN=2  # default number of processors

CONFIG=$1   # script's first argument is name of -config_override file
INFILE=$2   # script's second argument is PISM file with Greenland geometry
            #   and other needed info to run pclimate
OUTFILE=$3  # script's third argument is the name of the output of pclimate,
            #   which will be evaluated against smb from Ettema et al.

# change this if using other MPI, or to "" if you don't want MPI
MPIDO="mpiexec -n"

# coupler settings: Fausto 2m air temp parameterization, but default PDD
#   (w/o Greve/Fausto settings of PDD parameters)
COUPLER="-atmosphere searise_greenland -surface pdd"

CLIMSTARTTIME=1990
CLIMENDTIME=1991
#DT=0.0833333333 # monthly = (1/12) of year
DT=1.0
echo
echo "$SCRIPTNAME  running pclimate with -config_override file $CONFIG to generate"
echo "$SCRIPTNAME    $OUTFILE containing ice surface mass balance"

$MPIDO $NN pclimate -i $INFILE $GRID $COUPLER \
  -config_override $CONFIG -ys $CLIMSTARTTIME -ye $CLIMENDTIME -dt $DT -o $OUTFILE

echo

