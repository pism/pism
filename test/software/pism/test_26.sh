#!/bin/bash

source ../functions.sh

# Test name:
test="Test #26: Schoof (2003) bed roughness SIA parameterization regression."
# The list of files to delete when done.
files="brout.txt diff.txt"
dir=`pwd`

run_test ()
{
    cleanup

    bedrough_test > brout.txt
    # compare results
    diff brout.txt - > diff.txt <<END-OF-OUTPUT
PISMBedSmoother TEST
  smoothing domain:  Nx = 2, Ny = 2
  original bed    :  min elev =  -500.000000 m,  max elev =   500.000000 m
  smoothed bed    :  min elev =  -372.992474 m,  max elev =   372.992474 m
  Schoof's theta  :  min      =  0.714730065,    max      =  0.988484365
END-OF-OUTPUT

    if [ $? != 0 ];
    then
	fail "bedrough_test output does not match the one stored"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
