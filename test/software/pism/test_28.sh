#!/bin/bash

source ../functions.sh

# Test name:
test="Test #28: Hooke flow law regression."
# The list of files to delete when done.
files="flowtable.txt diff.txt"
dir=`pwd`

run_test ()
{
    cleanup

    flowlaw_test -ice_type hooke > flowtable.txt
    # compare results
    diff flowtable.txt - > diff.txt <<END-OF-OUTPUT
flow law:   "hooke"
flowtable:  [pressure =   2.00e+07 throughout]
  (stress)   (enthalpy)    (temp)     =   (flow)
  1.00e+04   560000.000    271.570000 = 5.323971e-15
  1.00e+04   550000.000    271.570000 = 5.323971e-15
  1.00e+04   540000.000    268.790443 = 3.028986e-16
  1.00e+04   530000.000    263.812842 = 1.430122e-16
  1.00e+04   520000.000    258.835241 = 7.091070e-17
  5.00e+04   560000.000    271.570000 = 1.330993e-13
  5.00e+04   550000.000    271.570000 = 1.330993e-13
  5.00e+04   540000.000    268.790443 = 7.572464e-15
  5.00e+04   530000.000    263.812842 = 3.575304e-15
  5.00e+04   520000.000    258.835241 = 1.772767e-15
  1.00e+05   560000.000    271.570000 = 5.323971e-13
  1.00e+05   550000.000    271.570000 = 5.323971e-13
  1.00e+05   540000.000    268.790443 = 3.028986e-14
  1.00e+05   530000.000    263.812842 = 1.430122e-14
  1.00e+05   520000.000    258.835241 = 7.091070e-15
  1.50e+05   560000.000    271.570000 = 1.197893e-12
  1.50e+05   550000.000    271.570000 = 1.197893e-12
  1.50e+05   540000.000    268.790443 = 6.815218e-14
  1.50e+05   530000.000    263.812842 = 3.217773e-14
  1.50e+05   520000.000    258.835241 = 1.595491e-14
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
