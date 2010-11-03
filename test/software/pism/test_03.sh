#!/bin/bash

source ../functions.sh

test="Test #3: no information loss on -y 0 runs (ignoring diagnostic variables)."
files="foo.nc bar.nc baz.nc"
dir=`pwd`

OPTS="-surface constant -o_size small"

test_03 ()
{
    cleanup

    set -e

    # Create a file to start from:
    run -n 2 pisms -Mx 61 -My 61 -Mz 61 -y 10 -no_cold -verbose 1 -o foo.nc

    # Run for a year using pismr (so that all the parameters, including
    # rheology, are the same):
    run -n 2 pismr -i foo.nc -y 1 $OPTS -o bar.nc

    # Run for 0 years:
    run -n 2 pismr -i bar.nc -y 0 $OPTS -o baz.nc

    set +e

    # Compare, excluding irrelevant diagnostic variables:
    run nccmp.py -x -v t bar.nc baz.nc
    if [ $? != 0 ];
    then
	fail "foo.nc and bar.nc are different."
	return 1
    fi

    pass
    return 0
}

test_03
