#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test #10: regridding with different number of processors."

files="foo0.nc foo1.nc foo2.nc foo3.nc foo4.nc foo5.nc"

rm -f $files

set -e -x

NRANGE="1 2 3 4 5"

# Create a file to bootstrap from:
$MPIEXEC -n 1 $PISM_PATH/pisms -eisII I -Mx 51 -My 101 -y 0 -o foo0.nc

# Bootstrap:
for NN in $NRANGE;
do 
    $MPIEXEC -n $NN $PISM_PATH/pismr -boot_file foo0.nc -Mx 51 -My 101 -Mz 11 -Lz 5000 -y 0 -o foo$NN.nc -o_size small
done

set +e

# Compare:
for i in $NRANGE;
do
    for j in $NRANGE;
    do
	if [ $i -le $j ]; then continue; fi
	
	$PISM_PATH/nccmp.py foo$i.nc foo$j.nc
	if [ $? != 0 ];
	then
	    exit 1
	fi
    done
done

rm -f $files; exit 0

