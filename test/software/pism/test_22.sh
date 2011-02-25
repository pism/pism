#!/bin/bash

source ../functions.sh

# Test name:
test="Test #22: EISMINT-Ross experiment regression."
# The list of files to delete when done.
files="ross.nc riggs.nc rossComputed.nc ross.txt"
dir=`pwd`

run_test ()
{
    cleanup

    cp ../../../examples/eisross/test/ross.nc.gz .
    cp ../../../examples/eisross/test/riggs.nc.gz .
    gunzip ross.nc.gz
    gunzip riggs.nc.gz

    $MPIDO -n 2 pross -boot_file ross.nc -Mx 147 -My 147 \
        -riggs riggs.nc -o rossComputed.nc > ross.txt

    python <<EOF
from numpy import double, abs
from sys import exit
chi_squared = 1e6
good_chi_squared = 3680.748
rel_tolerance = 0.01

f = open("ross.txt")
for line in f:
  words = line.split(' ')
  if words[0] == "Chi^2":
    chi_squared = double(words[-1])
    break

rel_difference = abs(chi_squared - good_chi_squared) / good_chi_squared
if rel_difference < rel_tolerance:
#  print "Chi^2 compared: %f and %f. Difference: %f%%" % (chi_squared, good_chi_squared, rel_difference*100)
  exit(0)
else:
  print "Chi^2 = %f, should be near %f" % (chi_squared, good_chi_squared)
  exit(1)
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
