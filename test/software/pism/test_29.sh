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
pressure = 1.785e+07 Pa = (hydrostatic at depth 2000.00 m)
flowtable:
  (dev stress)   (abs temp) (liq frac) =   (flow)
      1.00e+04      241.740      0.000 = 7.324397e-17
      1.00e+04      266.740      0.000 = 5.496298e-15
      1.00e+04      271.740      0.000 = 2.417138e-14
      1.00e+04      271.740      0.005 = 2.417138e-14
      5.00e+04      241.740      0.000 = 2.163601e-16
      5.00e+04      266.740      0.000 = 1.934468e-14
      5.00e+04      271.740      0.000 = 9.044284e-14
      5.00e+04      271.740      0.005 = 9.044284e-14
      1.00e+05      241.740      0.000 = 4.061917e-16
      1.00e+05      266.740      0.000 = 3.397701e-14
      1.00e+05      271.740      0.000 = 1.605747e-13
      1.00e+05      271.740      0.005 = 1.605747e-13
      1.50e+05      241.740      0.000 = 6.689768e-16
      1.50e+05      266.740      0.000 = 4.807048e-14
      1.50e+05      271.740      0.000 = 2.278162e-13
      1.50e+05      271.740      0.005 = 2.278162e-13
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
