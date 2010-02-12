#!/bin/bash

source ../functions.sh

# Test name:
test="Test #11: automatic vertical grid extension."
# The list of files to delete when done.
files="foo.nc bar.nc baz.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # run with Lz set too low:
    run -n 2 pisms -eisII A -y 1000 -Mmax 0.925 -Lz 900 -z_spacing equal -o foo.nc 

    # run with Lz set just right:
    run -n 2 pisms -eisII A -y 1000 -Mz 33 -Mmax 0.925 -Lz 960 -z_spacing equal -o bar.nc 

    # regrid from the extended grid onto the one in bar.nc:
    run -n 2 pismr -i bar.nc -surface constant -regrid_from foo.nc -regrid_vars temp -y 0 -o baz.nc

    # compare results
    run nccmp.py -v temp bar.nc baz.nc

    if [ $? != 0 ];
    then
	fail "temperature fields in bar.nc and baz.nc are different"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
