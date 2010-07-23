#!/bin/bash

source ../functions.sh

# Test name:
test="Test #24: GPBLD flow law regression."
# The list of files to delete when done.
files="flowtable.txt diff.txt"
dir=`pwd`

run_test ()
{
    cleanup

    flowlaw_test -ice_type gpbld > flowtable.txt
    # compare results
    diff flowtable.txt - > diff.txt <<END-OF-OUTPUT
flow law:   "gpbld"
flowtable:  [pressure =   2.00e+07 throughout]
  (stress)   (enthalpy)    (temp)     =   (flow)
  1.00e+04   560000.000    271.570000 = 4.062203e-15
  1.00e+04   550000.000    271.570000 = 1.559473e-15
  1.00e+04   540000.000    268.790443 = 2.421466e-16
  1.00e+04   530000.000    263.812842 = 7.593050e-17
  1.00e+04   520000.000    258.835241 = 3.332262e-17
  5.00e+04   560000.000    271.570000 = 1.015551e-13
  5.00e+04   550000.000    271.570000 = 3.898683e-14
  5.00e+04   540000.000    268.790443 = 6.053665e-15
  5.00e+04   530000.000    263.812842 = 1.898262e-15
  5.00e+04   520000.000    258.835241 = 8.330654e-16
  1.00e+05   560000.000    271.570000 = 4.062203e-13
  1.00e+05   550000.000    271.570000 = 1.559473e-13
  1.00e+05   540000.000    268.790443 = 2.421466e-14
  1.00e+05   530000.000    263.812842 = 7.593050e-15
  1.00e+05   520000.000    258.835241 = 3.332262e-15
  1.50e+05   560000.000    271.570000 = 9.139957e-13
  1.50e+05   550000.000    271.570000 = 3.508814e-13
  1.50e+05   540000.000    268.790443 = 5.448299e-14
  1.50e+05   530000.000    263.812842 = 1.708436e-14
  1.50e+05   520000.000    258.835241 = 7.497589e-15
END-OF-OUTPUT

    if [ $? != 0 ];
    then
	fail "flowtable output does not match the one stored"
	# the return statement *is* needed here, because 'fail' does not
	# terminate the test execution
	return 1
    fi

    pass
    return 0
}

run_test
