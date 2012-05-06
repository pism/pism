#!/bin/bash

# Copyright (C) 2011-2012 the PISM authors

# downloads 5km data, including a precomputed PISM state for the whole
# ice sheet and the SeaRISE 5km Greenland data set, from which we will only
# grab the "smb" = Surface Mass Balance field, which was produced by
# the regional atmosphere model RACMO

# depends on wget and NCO (ncrename, ncap, ncwa, ncks)

set -e  # exit on error

echo "fetching pre-computed whole ice-sheet result on 5km grid"
URL=http://www.pism-docs.org/download
WHOLE=g5km_0_ftt.nc
wget -nc ${URL}/$WHOLE
echo "... done"
echo

BCFILE=g5km_bc.nc
echo "creating PISM-readable boundary conditions file $BCFILE"
echo "   from whole ice sheet model result ..."
ncks -O -v u_ssa,v_ssa,bmelt,bwat,enthalpy,litho_temp $WHOLE $BCFILE
# rename u_ssa and v_ssa so that they are specified as b.c.
ncrename -O -v u_ssa,u_ssa_bc -v v_ssa,v_ssa_bc $BCFILE
echo "... done with creating bc file $BCFILE"
echo

