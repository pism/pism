#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
test="Test #14: enthalpy symmetry near the base (pisms -no_cold)."
# The list of files to delete when done.
files="simp_exper.nc"

rm -f $files
# run pisms
$MPIEXEC -n 2 $PISM_PATH/pisms -y 10e3 -Lz 4100 -Mx 26 -My 26 -o_size big -no_cold

/usr/bin/env python <<EOF
try:
    from netCDF3 import Dataset as NC
except:
    from netCDF4 import Dataset as NC
from numpy import abs, arange
from sys import exit

nc = NC("simp_exper.nc", 'r')
var = nc.variables['enthalpy']
n = 26; m = 26; tol = 1e-3

for k in [0, 1, 2]:
    v = var[0,k,:,:]	# time,z,y,x
    for i in arange((n-1)/2):
        for j in arange((m-1)/2):
            ii = (n-1) - i
            jj = (m-1) - j

            delta = abs(v[i,j] - v[ii,j])
            if (delta >= tol):
                print "X-symmetry failure at (%d,%d),(%d,%d) level %d (delta = %2.2e)" % (i,j,ii,j,k,delta)
                exit(1)

            delta = abs(v[i,j] - v[i,jj])
            if (delta >= tol):
                print "Y-symmetry failure at (%d,%d),(%d,%d) level %d (delta = %2.2e)" % (i,j,i,jj,k,delta)
                exit(1)
                
            delta = abs(v[i,j] - v[ii,jj])
            if (delta >= tol):
                print "Radial symmetry failure at (%d,%d),(%d,%d) level %d (delta = %2.2e)" % (i,j,ii,jj,k,delta)
                exit(1)
exit(0)
EOF

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

