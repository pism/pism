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
  wget -nc http://homepages.vub.ac.be/~phuybrec/eismint/$fname
done

echo "-----  Run eis_ross.py to turn ascii data into NetCDF file ross.nc:"
./eisross.py -o ross.nc

echo "-----  Also create NetCDF version of RIGGS data; riggs.nc; used later:"
./eisriggs.py -o riggs.nc

echo "-----  Running pismd to compute velocity in Ross ice shelf:"
mpiexec -n $NN pismd -ross -boot_from ross.nc -ssa -ssaBC ross.nc -Mx 147 -My 147 -Mz 3 -Lz 1e3 -o rossComputed.nc

exit 

# FIXME:  we want to do this compare to RIGGS data by default, too:

echo "-----  Running pismd to compute velocity in Ross ice shelf, with compare to RIGGS data:"
mpiexec -n $NN pismd -ross -boot_from ross.nc \
 -ssa -ssa_rtol 1e-7 -ssaBC ross.nc -riggs riggs.nc \
 -Mx 147 -My 147 -Mz 3 -Lz 1e3 -o rossComputed.nc
  
echo "-----  Done.  Model output in rossComputed.nc."

