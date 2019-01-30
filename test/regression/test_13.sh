#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #13: enthalpy symmetry near the base (pisms -energy enthalpy)."
# The list of files to delete when done.
files="simp_exper-13.nc"

rm -f $files
# run pisms
$MPIEXEC -n 2 $PISM_PATH/pisms -y 10e3 -Lz 5000 -Mx 12 -My 12 -o_size big -energy enthalpy -o simp_exper-13.nc -output.variable_order zyx

/usr/bin/env python <<EOF
from netCDF4 import Dataset as NC
from numpy import abs, arange
from sys import exit

nc = NC("simp_exper-13.nc", 'r')
var = nc.variables['enthalpy']
n = 12; m = 12; tol = 1e-3

for k in [0, 1, 2]:
    v = var[0,k,:,:]	# time,z,y,x
    for i in arange((n-1)/2, dtype=int):
        for j in arange((m-1)/2, dtype=int):
            ii = (n-1) - i
            jj = (m-1) - j

            delta = abs(v[i,j] - v[ii,j])
            if (delta >= tol):
                print("X-symmetry failure at (%d,%d),(%d,%d) level %d (delta = %2.2e)" % (i,j,ii,j,k,delta))
                exit(1)

            delta = abs(v[i,j] - v[i,jj])
            if (delta >= tol):
                print("Y-symmetry failure at (%d,%d),(%d,%d) level %d (delta = %2.2e)" % (i,j,i,jj,k,delta))
                exit(1)
                
            delta = abs(v[i,j] - v[ii,jj])
            if (delta >= tol):
                print("Radial symmetry failure at (%d,%d),(%d,%d) level %d (delta = %2.2e)" % (i,j,ii,jj,k,delta))
                exit(1)
exit(0)
EOF

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

