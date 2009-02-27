#!/bin/bash
#   Downloads data and runs an EISMINT-Ross example in PISM.  See User's Manual.

NN=1  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "./quickstart.sh 8" then NN = 8
  NN="$1"
fi

set -e  # exit on error

echo "-----  Download the ASCII files from the EISMINT web site:"
for fname in "111by147Grid.dat" "kbc.dat" "inlets.dat"
do
  wget http://homepages.vub.ac.be/~phuybrec/eismint/$fname
done

echo "-----  Run eis_ross.py to turn ascii data into NetCDF file ross.nc:"
./eisross.py -o ross.nc

echo "-----  Also create NetCDF version of RIGGS data; riggs.nc; used later:"
./eisriggs.py -o riggs.nc

echo "-----  Running PISM to compute velocity in Ross ice shelf:"
mpiexec -n $NN pismd -ross -boot_from ross.nc -ssa -ssaBC ross.nc -Mx 147 -My 147 -Mz 3 -Lz 1e3 -o rossComputed.nc
 
echo "-----  Done.  Model output in rossComputed.nc."

