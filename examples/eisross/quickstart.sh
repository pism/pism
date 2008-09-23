#!/bin/bash

#   Downloads data and runs an EISMINT-Ross example in PISM.

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "ccl3.sh 8" then NN = 8
  NN="$1"
fi
set -e  # exit on error

echo "--------------------------------------------------"
echo "Download the ASCII files from the EISMINT web site:"
for fname in "111by147Grid.dat" "kbc.dat" "inlets.dat"
do
  wget http://homepages.vub.ac.be/~phuybrec/eismint/$fname
done

echo
echo "--------------------------------------------------"
echo "Run eis_ross.py to turn ascii data into NetCDF file:"
./eis_ross.py -o ross.nc


echo
echo "--------------------------------------------------"
echo "Running PISM to compute velocity in Ross ice shelf:"
pismd -ross -bif ross.nc -ssaBC ross.nc -Mx 147 -My 147 -Lz 1000 -Mz 11 \
 -ssa -o rossComputed
 
echo 
echo "--------------------------------------------------"
echo "Done.  View input NetCDF file ross.nc and model output"
echo "rossComputed.nc."
echo 

