#!/bin/bash

source ../functions.sh

# Test name:
test="Test #31: testing whether runtime viewers break or not."
# The list of files to delete when done.
files="simp_exper.nc"
dir=`pwd`

run_test ()
{
    cleanup

    # do stuff
    run pisms -eisII A -y 1000 -view_sounding temp,litho_temp -view_map velsurf,thk -Mbz 11 -Lbz 1000 -o_size small

    if [ $? != 0 ];
    then
	fail "pisms with runtime viewers on crashed"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
