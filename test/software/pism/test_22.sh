#!/bin/bash

source ../functions.sh

# Test name:
test="Test #22: SeaRISE-Greenland regression: 40km, paleo, SIA+SSA, topography, tillphi, enthalpy."
# The list of files to delete when done.
files="unnamed.nc g40km_m124980.nc g40km_SIA.nc pism_dT.nc"
dir=`pwd`

run_test ()
{
    cleanup

    set +e
    
    # prepare
    cp ../../../examples/searise-greenland/test/pism_dT.nc .
    cp ../../../examples/searise-greenland/test/g40km_SIA.nc.gz .
    cp ../../../examples/searise-greenland/test/g40km_m124980.nc.gz .
    gunzip g40km_SIA.nc.gz
    gunzip g40km_m124980.nc.gz

    # run pismr with SeaRISE-Greenland options:
    run -n 2 pismr -ocean_kill -e 3 -skip 5 -i g40km_SIA.nc \
      -topg_to_phi 5.0,20.0,-300.0,700.0,10.0 -ssa_sliding -thk_eff \
      -pseudo_plastic_q 0.25 -plastic_pwfrac 0.98 \
      -atmosphere searise_greenland,forcing -dTforcing pism_dT.nc \
      -surface pdd -pdd_fausto -ocean constant -ys -125000 -y 20

    # compare key variables in result:
    run nccmp.py -v bmelt,bwat,tauc,tauc,thk,csurf,cbase,mask unnamed.nc g40km_m124980.nc
    if [ $? != 0 ];
    then
	fail "recomputed and stored results of SeaRISE-Greenland runs are different."
	return 1
    fi

    pass
    return 0
}

run_test

