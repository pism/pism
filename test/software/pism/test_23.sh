#!/bin/bash

source ../functions.sh

# Test name:
test="Test #23: PST regression P1 on 25km grid: SIA+SSA hybrid, flat bed, cold ice."
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

    # compare key variables in result:
    run nccmp.py -v bmelt,bwat,tauc,tauc,thk,csurf,cbase simp_exper.nc P1_small50.nc
    if [ $? != 0 ];
    then
	fail "recomputed and stored results of short P1 runs are different."
	return 1
    fi

    pass
    return 0
}

run_test

