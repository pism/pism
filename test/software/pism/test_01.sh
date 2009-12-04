#!/bin/bash

source ../functions.sh

test="Test #1: pismr exact restartability."
files="verify.nc foo.nc joe.nc bar.nc"
dir=`pwd`

OPTS="-max_dt 1"

test_01 ()
{
    cleanup

    set -e

    # generate an interesting file
    run -n 1 pismv -test G -y 10 -o verify.nc

    # run for ten years, fixed time step
    run -n 1 pismr -i verify.nc $OPTS -y 10 -o foo.nc

    # chain two five year runs, fixed time step
    run -n 1 pismr -i verify.nc $OPTS -y 5 -o joe.nc
    run -n 1 pismr -i joe.nc    $OPTS -y 5 -o bar.nc

    set +e

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
