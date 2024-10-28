#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test #10: regridding with different number of processors."

files="foo0-10.nc foo1-10.nc foo2-10.nc foo3-10.nc foo4-10.nc foo5-10.nc"

rm -f $files

set -e -x

NRANGE="1 2 3 4 5"

# Create a file to bootstrap from:
$MPIEXEC -n 1 $PISM_PATH/pism -eisII A -eisII I -Mx 51 -My 101 -y 0 -o foo0-10.nc

# Bootstrap:
for NN in $NRANGE;
do 
    $MPIEXEC -n $NN $PISM_PATH/pism -i foo0-10.nc -bootstrap -Mx 51 -My 101 -Mz 11 -Lz 5000 -y 0 -o foo$NN-10.nc -o_size small
done

set +e

# Compare:
for i in $NRANGE;
do
    for j in $NRANGE;
    do
	if [ $i -le $j ]; then continue; fi
	
	$PISM_PATH/nccmp.py -x -v timestamp foo$i-10.nc foo$j-10.nc
	if [ $? != 0 ];
	then
	    exit 1
	fi
    done
done

rm -f $files; exit 0

