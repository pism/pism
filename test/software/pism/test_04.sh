#!/bin/bash

source ../functions.sh

test="Test #4: regridding during bootstrapping: coarse -> fine -> coarse."
files="foo.nc bar.nc baz.nc"
dir=`pwd`

test_04 ()
{
    cleanup

    # Create a file to bootstrap from:
    run pismv -test G -Mx 11 -My 11 -Mz 11 -y 0 -o foo.nc

    # Coarse -> fine:
    run pismr -boot_from foo.nc -Mx 21 -My 21 -Mz 21 -Lz 4000 -y 0 -o bar.nc

    # Fine -> coarse:
    run pismr -boot_from bar.nc -Mx 11 -My 11 -Mz 11 -Lz 4000 -y 0 -o baz.nc

    # Compare:
    run nccmp.py -v topg foo.nc baz.nc
    if [ $? != 0 ];
    then
	fail "foo.nc and baz.nc are different."
	return 1
    fi

    pass
    return 0
}

test_04
