#!/bin/bash

source ../functions.sh

test="Test #3: no information loss on -y 0 runs (ignoring diagnostic variables)."
files="foo.nc bar.nc"
dir=`pwd`

test_03 ()
{
    cleanup

    set -e

    # Create a file to start from:
    run -n 2 pismv -test G -Mx 61 -My 61 -Mz 61 -y 100 -verbose 1 -o foo.nc

    # Run for 0 years:
    run -n 2 pismr -i foo.nc -y 0 -o bar.nc

    set +e

    # Compare, excluding irrelevant diagnostic variables:
    run nccmp.py -x -v usurf,dHdt,cbar,cflx,cbase,csurf,wvelsurf,taud,tauc,dhdt foo.nc bar.nc
    if [ $? != 0 ];
    then
	fail "foo.nc and bar.nc are different."
	return 1
    fi

    pass
    return 0
}

test_03
