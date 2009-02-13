#!/bin/bash

source ../functions.sh

test="Test #4: regridding: coarse -> fine -> coarse."
files="foo.nc bar.nc baz.nc"
dir=`pwd`

test_04 ()
{
    cleanup

    run pismv -test G -Mx 11 -My 11 -Mz 11 -y 0 -o foo.nc
    run pismr -boot_from foo.nc -Mx 21 -My 21 -Mz 21 -y 0 -o bar.nc
    run pismr -boot_from bar.nc -Mx 11 -My 11 -Mz 11 -y 0 -o baz.nc

    run nccmp.py -t 1e-16 -v topg foo.nc baz.nc
    if [ ! $? ];
    then
	fail "foo.nc and baz.nc are different."
	return 1
    fi

    pass
    return 0
}

test_04