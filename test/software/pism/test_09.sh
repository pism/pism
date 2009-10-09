#!/bin/bash

source ../functions.sh

test="Test #9: 3D regridding from files with different variable orders."
files="foo.nc bar.nc baz.nc"
dir=`pwd`

run_test ()
{
    cleanup

    set -e

    # Create a file to bootstrap from (with a non-trivial bed topography):
    run -n 1 pisms -eisII I -Mx 61 -My 61 -Mz 201 -Mbz 21 -Lbz 1000 -y 0 -o foo.nc

    # Bootstrap from this file and run for 0 years:
    run -n 2 pismr -boot_from foo.nc -Mx 61 -My 61 -Mz 201 -Lz 5000 -regrid_from foo.nc -regrid_vars bBehHLT -y 0 -o bar.nc

    # Change the variable order in foo.nc to z,y,x:
    run ncpdq -O -a z,y,x foo.nc foo.nc

    # Bootstrap from this file and run for 0 years:
    run -n 2 pismr -boot_from foo.nc -Mx 61 -My 61 -Mz 201 -Lz 5000 -regrid_from foo.nc -regrid_vars bBehHLT -y 0 -o baz.nc

    set +e

    # Compare bar.nc and baz.nc:
    run nccmp.py bar.nc baz.nc
    if [ $? != 0 ];
    then
	fail "files bar.nc and baz.nc are different"
	return 1
    fi

    pass
    return 0
}

run_test
