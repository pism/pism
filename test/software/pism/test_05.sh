#!/bin/bash

source ../functions.sh

test="Test #5: bootstrapping from files with different variable orders."
files="foo.nc bar.nc baz.nc"
dir=`pwd`

test_05 ()
{
    cleanup

    # Create a file to bootstrap from (with a non-trivial bed topography):
    run pisms -eisII I -Mx 61 -My 61 -Mz 201 -y 0 -o foo.nc

    # Bootstrap from this file and run for 0 years:
    run mpiexec -n 2 pismr -boot_from foo.nc -y 0 -o bar.nc

    # Change the variable order in foo.nc to z,y,x:
    run ncpdq -O -a z,y,x foo.nc foo.nc

    # Bootstrap from this file and run for 0 years:
    run mpiexec -n 2 pismr -boot_from foo.nc -y 0 -o baz.nc

    # Compare bar.nc and baz.nc:
    run nccmp.py bar.nc baz.nc
    if [ ! $? ];
    then
	fail "files bar.nc and baz.nc are different"
	return 1
    fi

    pass
    return 0
}

test_05
