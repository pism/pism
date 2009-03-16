#!/bin/bash

source ../functions.sh

# Test name:
test="Test #8: regridding: coarse -> fine -> coarse (in the vertical direction)."
# The list of files to delete when done.
files="coarse1.nc coarse2.nc fine1.nc fine2.nc"
dir=`pwd`

run_test ()
{
    cleanup

    set -e

    # Create a file to regrid from:
    run -n 1 pismv -test G -Mx 11 -My 11 -Mz 11 -y 0 -o coarse1.nc
    # Create another file with a finer grid:
    run -n 1 pismv -test G -Mx 11 -My 11 -Mz 21 -y 0 -o fine1.nc

    # Coarse -> fine:
    run -n 1 pismr -i fine1.nc -regrid_from coarse1.nc -regrid_vars T -y 0 -o fine2.nc
    # Fine -> coarse:
    run -n 1 pismr -i coarse1.nc -regrid_from fine2.nc -regrid_vars T -y 0 -o coarse2.nc

    set +e

    # Compare:
    run nccmp.py -t 1e-16 -v temp coarse1.nc coarse2.nc
    if [ $? != 0 ];
    then
	fail "files coarse1.nc and coarse2.nc are different"
    fi

    pass
    return 0
}

run_test
