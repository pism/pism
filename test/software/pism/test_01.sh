#!/bin/bash

source ../functions.sh

test="Test #1: pismr exact restartability."
files="verify.nc foo.nc joe.nc bar.nc"
dir=`pwd`

test_01 ()
{
    cleanup

    # generate an interesting file
    run pismv -test G -y 10 -o verify.nc

    # run for ten years, fixed time step
    run pismr -i verify.nc -max_dt 1 -y 10 -o foo.nc

    # chain two five year runs, fixed time step
    run pismr -i verify.nc -max_dt 1 -y 5 -o joe.nc
    run pismr -i joe.nc -max_dt 1 -y 5 -o bar.nc

    # Compare output files at year 10:
    run nccmp.py -t 1e-6 foo.nc bar.nc
    if [ $? != 0 ];
    then
	fail "Output files are different."
	return 1
    fi

    pass
    return 0
}

test_01
