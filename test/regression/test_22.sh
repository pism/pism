#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# Test name:
echo "Test #22: EISMINT-Ross experiment regression."
# The list of files to delete when done.
files="ross.nc riggs.nc rossComputed.nc ross.txt"

rm -f $files

cp $PISM_SOURCE_DIR/examples/eisross/test/ross.nc.gz .
cp $PISM_SOURCE_DIR/examples/eisross/test/riggs.nc.gz .
gunzip ross.nc.gz
gunzip riggs.nc.gz

$MPIEXEC -n 2 $PISM_PATH/pross -boot_file ross.nc -Mx 147 -My 147 \
    -riggs riggs.nc -o rossComputed.nc > ross.txt

/usr/bin/env python <<EOF
from numpy import double, abs
from sys import exit
chi_squared = 1e6
good_chi_squared = 3649.427
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
    exit 1
fi

rm -f $files; exit 0

