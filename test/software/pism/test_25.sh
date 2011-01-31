#!/bin/bash

source ../functions.sh

# Test name:
test="Test #25: Paterson-Budd flow law regression."
# The list of files to delete when done.
files="flowtable.txt diff.txt"
dir=`pwd`

run_test ()
{
    cleanup

    flowlaw_test -ice_type pb > flowtable.txt
    # compare results
    diff flowtable.txt - > diff.txt <<END-OF-OUTPUT
flow law:   "pb"
pressure = 1.785e+07 Pa = (hydrostatic at depth 2000.00 m)
flowtable:
  (dev stress)   (abs temp) (liq frac) =   (flow)
      1.00e+04      241.740      0.000 = 4.657918e-18
      1.00e+04      266.740      0.000 = 1.451147e-16
      1.00e+04      271.740      0.000 = 4.542999e-16
      1.00e+04      271.740      0.005 = 4.542999e-16
      5.00e+04      241.740      0.000 = 1.164479e-16
      5.00e+04      266.740      0.000 = 3.627868e-15
      5.00e+04      271.740      0.000 = 1.135750e-14
      5.00e+04      271.740      0.005 = 1.135750e-14
      1.00e+05      241.740      0.000 = 4.657918e-16
      1.00e+05      266.740      0.000 = 1.451147e-14
      1.00e+05      271.740      0.000 = 4.542999e-14
      1.00e+05      271.740      0.005 = 4.542999e-14
      1.50e+05      241.740      0.000 = 1.048031e-15
      1.50e+05      266.740      0.000 = 3.265081e-14
      1.50e+05      271.740      0.000 = 1.022175e-13
      1.50e+05      271.740      0.005 = 1.022175e-13
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
