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
pressure = 1.785e+07 Pa = (hydrostatic at depth 2000.00 m)
flowtable:
  (dev stress)   (abs temp) (liq frac) =   (flow)
      1.00e+04      241.740      0.000 = 5.267759e-18
      1.00e+04      266.740      0.000 = 2.123259e-16
      1.00e+04      271.740      0.000 = 5.323971e-15
      1.00e+04      271.740      0.005 = 5.323971e-15
      5.00e+04      241.740      0.000 = 1.316940e-16
      5.00e+04      266.740      0.000 = 5.308148e-15
      5.00e+04      271.740      0.000 = 1.330993e-13
      5.00e+04      271.740      0.005 = 1.330993e-13
      1.00e+05      241.740      0.000 = 5.267759e-16
      1.00e+05      266.740      0.000 = 2.123259e-14
      1.00e+05      271.740      0.000 = 5.323971e-13
      1.00e+05      271.740      0.005 = 5.323971e-13
      1.50e+05      241.740      0.000 = 1.185246e-15
      1.50e+05      266.740      0.000 = 4.777333e-14
      1.50e+05      271.740      0.000 = 1.197893e-12
      1.50e+05      271.740      0.005 = 1.197893e-12
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
