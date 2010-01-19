#!/bin/bash

source ../functions.sh

test="Test #12: penth exact restartability."
files="pre.nc blah.nc foo.nc joe.nc bar.nc"
dir=`pwd`

OPTS="-max_dt 1"

run_test ()
{
    cleanup

    set -e

    # generate an interesting file
    run -n 1 pismv -test G -y 10 -o pre.nc

    # bootstrap it to get variable 'enthalpy' into blah.
    run -n 1 penth -boot_from pre.nc -Mx 31 -My 31 -Mz 31 -Lz 4000 -y 1 -o blah.nc

    # run for 10 years, fixed time step
    run -n 1 penth -i blah.nc $OPTS -y 10 -o foo.nc

    # chain two 5-year runs, fixed time step
    run -n 1 penth -i blah.nc $OPTS -y 5  -o joe.nc
    run -n 1 penth -i joe.nc  $OPTS -y 5  -o bar.nc

    set +e

    # Compare output files at year 20:
    run nccmp.py foo.nc bar.nc
    if [ $? != 0 ];
    then
	fail "Output files are different."
	return 1
    fi

    pass
    return 0
}

run_test
