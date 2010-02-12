#!/bin/bash

source ../functions.sh

# Test name:
test="Test #7: comparing regridding via -boot_from and -regrid_from."
# The list of files to delete when done:
files="foo.nc bar.nc baz.nc"
dir=`pwd`

OPTS="-y 0 -surface constant"
run_test ()
{
    cleanup

    set -e

    # Create the file to regrid and bootstrap from:
    run -n 1 pismv -test G -Lx 4000 -Ly 4000 -Lz 4000 -Mx 41 -My 41 -Mz 41 -y 0 -o foo.nc

    # Bootstrap from this file:
    run -n 1 pismr -boot_from foo.nc -Lx 2000 -Ly 2000 -Lz 4000 -Mx 41 -My 41 -Mz 41 $OPTS -o bar.nc

    # Overwrite topg using -regrig_from and save the result to baz.nc:
    run -n 1 pismr -i bar.nc -regrid_from foo.nc -regrid_vars topg $OPTS -o baz.nc

    set +e

    # Compare:
    run nccmp.py -v topg bar.nc baz.nc
    if [ $? != 0 ];
    then
	fail "files bar.nc and baz.nc are different"
	return 1
    fi

    pass
    return 0
}

run_test
