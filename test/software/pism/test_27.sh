#!/bin/bash

source ../functions.sh

# Test name:
test="Test #27: comparing restarting with \"-i\" to \"-boot_file ... -regrid_file ...\""
# The list of files to delete when done.
files="foo.nc bar.nc"
dir=`pwd`

bedrock="-Mbz 11 -Lbz 1000"
opts="-y 0 -surface constant -Mx 61 -My 61 -Mz 31 -Lz 4000 $bedrock -o_size small"
run_test ()
{
    cleanup

    set -e

    # create foo.nc (at about 6500 years we get some basal melting...)
    run -n 2 pisms -no_cold -y 6500 $bedrock -o_size small -o foo.nc

    # bootstrap from it, re-gridding all the variables we can
    run pismr -boot_file foo.nc -regrid_file foo.nc $opts -o bar.nc -no_temp

    set +e

    # compare results (foo.nc and bar.nc contain model_state variables only)
    run nccmp.py foo.nc bar.nc
    if [ $? != 0 ];
    then
	fail "foo.nc and bar.nc are different."
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
