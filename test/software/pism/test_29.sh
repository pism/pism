#!/bin/bash

source ../functions.sh

# Test name:
test="Test #29: hybrid (GK-PB) flow law regression."
# The list of files to delete when done.
files="flowtable.txt diff.txt"
dir=`pwd`

run_test ()
{
    cleanup

    flowlaw_test -ice_type hybrid > flowtable.txt
    # compare results
    diff flowtable.txt - > diff.txt <<END-OF-OUTPUT
flow law:   "hybrid"
flowtable:  [pressure =   2.00e+07 throughout]
  (stress)   (enthalpy)    (temp)     =   (flow)
  1.00e+04   560000.000    271.570000 = 2.446307e-14
  1.00e+04   550000.000    271.570000 = 2.446307e-14
  1.00e+04   540000.000    268.790443 = 1.081282e-14
  1.00e+04   530000.000    263.812842 = 2.432228e-15
  1.00e+04   520000.000    258.835241 = 5.854443e-16
  5.00e+04   560000.000    271.570000 = 9.155356e-14
  5.00e+04   550000.000    271.570000 = 9.155356e-14
  5.00e+04   540000.000    268.790443 = 3.914932e-14
  5.00e+04   530000.000    263.812842 = 8.183158e-15
  5.00e+04   520000.000    258.835241 = 1.680155e-15
  1.00e+05   560000.000    271.570000 = 1.625522e-13
  1.00e+05   550000.000    271.570000 = 1.625522e-13
  1.00e+05   540000.000    268.790443 = 6.909137e-14
  1.00e+05   530000.000    263.812842 = 1.426425e-14
  1.00e+05   520000.000    258.835241 = 2.840796e-15
  1.50e+05   560000.000    271.570000 = 2.306246e-13
  1.50e+05   550000.000    271.570000 = 2.306246e-13
  1.50e+05   540000.000    268.790443 = 9.786462e-14
  1.50e+05   530000.000    263.812842 = 2.014066e-14
  1.50e+05   520000.000    258.835241 = 3.970141e-15
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
