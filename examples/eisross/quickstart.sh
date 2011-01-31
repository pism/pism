#!/bin/bash
#   Downloads data and runs an EISMINT-Ross example in PISM.
#   See User's Manual.

NN=2  # default number of processes
if [ $# -gt 0 ] ; then  # if user says "./quickstart.sh 8" then NN = 8
  NN="$1"
fi

set -e  # exit on error

echo "-----  Download the ASCII files from the EISMINT web site using wget:"
for fname in "111by147Grid.dat" "kbc.dat" "inlets.dat"
do
  wget -nc http://homepages.vub.ac.be/~phuybrec/eismint/$fname
done

echo "-----  Run eisross.py to turn ascii data from EISMINT into NetCDF file"
echo "-----    ross.nc (requires python modules numpy and netCDF3/netCDF4):"
./eisross.py -o ross.nc

echo "-----  Run riggs.py to create NetCDF version riggs.nc of RIGGS data"
echo "-----    (requires python modules numpy and netCDF3 or netCDF4):"
./riggs.py -o riggs.nc

echo "-----  Running 'pross' with $NN processes to compute velocity in"
echo "-----    Ross ice shelf, including comparison to RIGGS data:"
mpiexec -n $NN pross -boot_file ross.nc -Mx 147 -My 147 -ssaBC ross.nc \
        -riggs riggs.nc -o rossComputed.nc

echo "----- Generating figure comparing model vs observed velocity (requires"
echo "-----   python modules numpy, netCDF3/4, pylab, and scikits.delaunay):"
./rossplot.py

echo "-----  Done."

