#!/bin/bash

source ../functions.sh

# Test name:
test="Test #23: PST regression P1: SIA+SSA hybrid, flat bed, cold ice."
# The list of files to delete when done.
files="simp_exper.nc ser_pst_unnamedpst.nc P0A_small.nc P1_small50.nc"
dir=`pwd`

run_test ()
{
    cleanup

    set +e
    
    # prepare
    cp ../../../examples/pst/test/P0A_small.nc.gz .
    cp ../../../examples/pst/test/P1_small50.nc.gz .
    gunzip P0A_small.nc.gz
    gunzip P1_small50.nc.gz

    # run pisms for PST exper:
    run -n 2 pisms -pst -P1 -i P0A_small.nc -ys 0 -y 50 -skip 5

    # compare; goal is 7 digit agreement; comment at right gives scale
    run nccmp.py -v bwat -t 1e-7 simp_exper.nc P1_small50.nc # 1 m
    if [ $? != 0 ]; then
	fail "recomputed and stored bwat differ by more than 1e-7 relative."
	return 1; fi
    run nccmp.py -v tauc -t 1e-0 simp_exper.nc P1_small50.nc # 10 bar = 1e6 Pa
    if [ $? != 0 ]; then
	fail "recomputed and stored tauc differ by more than 1e-6 relative." # LOOSER
	return 2; fi
    run nccmp.py -v thk -t 1e-3 simp_exper.nc P1_small50.nc # 1000 m
    if [ $? != 0 ]; then
	fail "recomputed and stored thk differ by more than 1e-6 relative." # LOOSER
	return 3; fi
    run nccmp.py -v csurf -t 1e-3 simp_exper.nc P1_small50.nc # 100 m/a
    if [ $? != 0 ]; then
	fail "recomputed and stored csurf differ by more than 1e-5 relative." # LOOSER
	return 4; fi
    run nccmp.py -v cbase -t 1e-3 simp_exper.nc P1_small50.nc # 100 m/a
    if [ $? != 0 ]; then
	fail "recomputed and stored cbase differ by more than 1e-5 relative." # LOOSER
	return 5; fi

    pass
    return 0
}

run_test

