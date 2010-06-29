#!/bin/bash

source ../functions.sh

# Test name:
test="Test #22: SeaRISE-Greenland regression: paleo, SIA+SSA, topog., tillphi, enthalpy."
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

    # compare; goal is 7 digit agreement; comment at right gives scale
    run nccmp.py -v temp_pa -t 1e-7 unnamed.nc g40km_m124980.nc # 1 degree C
    if [ $? != 0 ]; then
	fail "recomputed and stored temp_pa differ by more than 1e-7 relative."
	return 1;  fi
    run nccmp.py -v liqfrac -t 1e-9 unnamed.nc g40km_m124980.nc # 0.01 [pure]
    if [ $? != 0 ]; then
	fail "recomputed and stored liqfrac differ by more than 1e-7 relative."
	return 2;  fi
    run nccmp.py -v thk -t 1e-4 unnamed.nc g40km_m124980.nc # 1000 m
    if [ $? != 0 ]; then
	fail "recomputed and stored thk differ by more than 1e-7 relative."
	return 3;  fi
    run nccmp.py -v csurf -t 1e-5 unnamed.nc g40km_m124980.nc # 100 m/a
    if [ $? != 0 ]; then
	fail "recomputed and stored csurf differ by more than 1e-7 relative."
	return 4;  fi
    run nccmp.py -v cbase -t 1e-5 unnamed.nc g40km_m124980.nc # 100 m/a
    if [ $? != 0 ]; then
	fail "recomputed and stored cbase differ by more than 1e-7 relative."
	return 5;  fi

    pass
    return 0
}

run_test

