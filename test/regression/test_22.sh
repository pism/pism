#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #22: Schoof (2003) bed roughness SIA parameterization regression."
# The list of files to delete when done.
files="brout-22.txt"

rm -f $files

$PISM_PATH/bedrough_test > brout-22.txt
# compare results
diff - brout-22.txt  <<END-OF-OUTPUT
BedSmoother TEST
  smoothing domain:  Nx = 2, Ny = 2
  original bed    :  min elev =  -500.000000 m,  max elev =   500.000000 m
  smoothed bed    :  min elev =  -372.992474 m,  max elev =   372.992474 m
  Schoof's theta  :  min      =  0.714730065,    max      =  0.988484365
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

