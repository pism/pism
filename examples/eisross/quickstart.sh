#!/bin/bash
#   Downloads data and runs an EISMINT-Ross example in PISM.
#   See User's Manual.

set -e  # exit on error

NN=2  # default number of processes
if [ $# -gt 0 ] ; then  # if user says "./quickstart.sh 8" then NN = 8
  NN="$1"
fi

MMxx=147
if [ $# -gt 1 ] ; then  # if user says "./quickstart.sh 8 31" then "-Mx 31"
  MMxx="$2"
fi

MMyy=147
if [ $# -gt 2 ] ; then  # if user says "./quickstart.sh 8 31 170" then "-Mx 31 -My 170"
  MMyy="$3"
fi

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

# use FD method normally, but allow user to do 'export METHOD=fem'
if [ -n "${METHOD:+1}" ] ; then
  echo "[METHOD = $METHOD  (already set)]"
else
  METHOD=fd
  echo "[METHOD = fd  (default)]"
fi

echo "-----  Running 'pross' with $NN processes and -ssa_method $METHOD to compute velocity in"
echo "-----    Ross ice shelf, including comparison to RIGGS data:"
set -x
mpiexec -n $NN pross -boot_file ross.nc -Mx $MMxx -My $MMyy -ssa_method $METHOD \
        -riggs riggs.nc -o rossComputed.nc
set +x

echo "----- Generating figure comparing model vs observed velocity (requires"
echo "-----   python modules numpy, netCDF3/4, and pylab:"
./rossplot.py --pism-output rossComputed.nc --riggs riggs_clean.dat

echo "-----  Done."

