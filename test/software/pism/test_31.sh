#!/bin/bash

source ../functions.sh

# Test name:
test="Test #31: Warm part of Paterson-Budd flow law regression."
# The list of files to delete when done.
files="flowtable.txt diff.txt"
dir=`pwd`

run_test ()
{
    cleanup

    flowlaw_test -ice_type arrwarm > flowtable.txt
    # compare results
    diff flowtable.txt - > diff.txt <<END-OF-OUTPUT
flow law:   "arrwarm"
flowtable:  [pressure =   2.00e+07 throughout]
  (stress)   (enthalpy)    (temp)     =   (flow)
  1.00e+04   560000.000    271.570000 = 3.181966e-16
  1.00e+04   550000.000    271.570000 = 3.181966e-16
  1.00e+04   540000.000    268.790443 = 1.683549e-16
  1.00e+04   530000.000    263.812842 = 5.206775e-17
  1.00e+04   520000.000    258.835241 = 1.539252e-17
  5.00e+04   560000.000    271.570000 = 7.954916e-15
  5.00e+04   550000.000    271.570000 = 7.954916e-15
  5.00e+04   540000.000    268.790443 = 4.208873e-15
  5.00e+04   530000.000    263.812842 = 1.301694e-15
  5.00e+04   520000.000    258.835241 = 3.848129e-16
  1.00e+05   560000.000    271.570000 = 3.181966e-14
  1.00e+05   550000.000    271.570000 = 3.181966e-14
  1.00e+05   540000.000    268.790443 = 1.683549e-14
  1.00e+05   530000.000    263.812842 = 5.206775e-15
  1.00e+05   520000.000    258.835241 = 1.539252e-15
  1.50e+05   560000.000    271.570000 = 7.159424e-14
  1.50e+05   550000.000    271.570000 = 7.159424e-14
  1.50e+05   540000.000    268.790443 = 3.787986e-14
  1.50e+05   530000.000    263.812842 = 1.171524e-14
  1.50e+05   520000.000    258.835241 = 3.463316e-15
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
