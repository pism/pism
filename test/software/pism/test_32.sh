#!/bin/bash

source ../functions.sh

# Test name:
test="Test #32: EISMINT-Ross experiment regression"
# The list of files to delete when done.
files="ross.nc riggs.nc rossComputed.nc ross.txt chi.txt diff.txt"
dir=`pwd`

run_test ()
{
    cleanup

    cp ../../../examples/eisross/test/ross.nc.gz .
    cp ../../../examples/eisross/test/riggs.nc.gz .
    gunzip ross.nc.gz
    gunzip riggs.nc.gz

    $MPIDO -n 2 pross -boot_file ross.nc -Mx 147 -My 147 -Mz 3 -Lz 1.5e3 -ssaBC ross.nc \
        -riggs riggs.nc -o rossComputed.nc > ross.txt

    cat ross.txt | grep "Chi^2" > chi.txt

    # compare results
    diff chi.txt - > diff.txt <<EOF
Chi^2 statistic for computed results compared to RIGGS is   5384.665
EOF

    if [ $? != 0 ];
    then
	fail "calculated and stored Chi^2 don't match"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
