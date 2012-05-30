#!/bin/bash

# Copyright (C) 2010 Ed Bueler

# see README for role of this script

SCRIPTNAME="#(boot.sh)"

set -e  # exit on error

NN=8  # default number of processors

PISMDATA=pism_Greenland_5km_v1.1.nc
OUTFILE=start.nc

## grids: choose one; there are only modest computational requirements for this
##        no-ice-dynamics anyway job 
FIVEKMGRID="-Mx 301 -My 561 -Lz 4000 -Lbz 2000 -Mz 41 -Mbz 11 -z_spacing equal -zb_spacing equal"
TENKMGRID="-Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 41 -Mbz 11 -z_spacing equal -zb_spacing equal"
TWENTYKMGRID="-Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 41 -Mbz 11 -z_spacing equal -zb_spacing equal"

## actual choice of grid
GRID=$TENKMGRID
GRIDNAME=10

# cat prefix and exec together
MPIDO="mpiexec -n"
PISMR="pismr -ocean_kill -e 3"

# coupler settings: Fausto 2m air temp parameterization, but default PDD
#   (w/o Greve/Fausto settings of PDD parameters)
COUPLER="-atmosphere searise_greenland -surface pdd"

# bootstrap and do very short smoothing run
RUN=1.0  # do one year; essentially no ice dynamics but slightly the smooths surface
echo
echo "$SCRIPTNAME  bootstrapping on $GRIDNAME km grid, for very short run (${RUN}a)"
$MPIDO $NN $PISMR -boot_file $PISMDATA $GRID $COUPLER -y $RUN -o $OUTFILE

echo

