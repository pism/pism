#!/bin/bash

source ../functions.sh

# Test name:
test="Test # 8: regridding: coarse -> fine -> coarse (in the vertical direction)."
# The list of files to delete when done.
files="coarse1.nc coarse2.nc fine1.nc fine2.nc"
dir=`pwd`

OPTS="-surface constant -y 0"
run_test ()
{
    cleanup

    set -e

    # Create a file to regrid from:
    run -n 1 pismv -test G -Mx 11 -My 11 -Mz 11 -y 0 -o coarse1.nc
    # Create another file with a finer grid:
    run -n 1 pismv -test G -Mx 11 -My 11 -Mz 21 -y 0 -o fine1.nc

    # Coarse -> fine:
    run -n 1 pismr -i fine1.nc -regrid_file coarse1.nc -regrid_vars enthalpy $OPTS -o fine2.nc
    # Fine -> coarse:
    run -n 1 pismr -i coarse1.nc -regrid_file fine2.nc -regrid_vars enthalpy $OPTS -o coarse2.nc

    set +e

    # Compare:
    run nccmp.py -t 1e-16 -v enthalpy coarse1.nc coarse2.nc
    if [ $? != 0 ];
    then
	fail "files coarse1.nc and coarse2.nc are different"
	return 0
    fi

    pass
    return 0
}

run_test
