#!/bin/bash

source ../functions.sh

# Test name:
test="Test #30: Cold part of Paterson-Budd flow law regression."
# The list of files to delete when done.
files="flowtable.txt diff.txt"
dir=`pwd`

run_test ()
{
    cleanup

    flowlaw_test -ice_type arr > flowtable.txt
    # compare results
    diff flowtable.txt - > diff.txt <<END-OF-OUTPUT
flow law:   "arr"
flowtable:  [pressure =   2.00e+07 throughout]
  (stress)   (enthalpy)    (temp)     =   (flow)
  1.00e+04   560000.000    271.570000 = 1.040083e-16
  1.00e+04   550000.000    271.570000 = 1.040083e-16
  1.00e+04   540000.000    268.790443 = 7.901848e-17
  1.00e+04   530000.000    263.812842 = 4.761380e-17
  1.00e+04   520000.000    258.835241 = 2.813686e-17
  5.00e+04   560000.000    271.570000 = 2.600208e-15
  5.00e+04   550000.000    271.570000 = 2.600208e-15
  5.00e+04   540000.000    268.790443 = 1.975462e-15
  5.00e+04   530000.000    263.812842 = 1.190345e-15
  5.00e+04   520000.000    258.835241 = 7.034215e-16
  1.00e+05   560000.000    271.570000 = 1.040083e-14
  1.00e+05   550000.000    271.570000 = 1.040083e-14
  1.00e+05   540000.000    268.790443 = 7.901848e-15
  1.00e+05   530000.000    263.812842 = 4.761380e-15
  1.00e+05   520000.000    258.835241 = 2.813686e-15
  1.50e+05   560000.000    271.570000 = 2.340187e-14
  1.50e+05   550000.000    271.570000 = 2.340187e-14
  1.50e+05   540000.000    268.790443 = 1.777916e-14
  1.50e+05   530000.000    263.812842 = 1.071310e-14
  1.50e+05   520000.000    258.835241 = 6.330794e-15
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
