#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 2: pismv exact processor independence (SIA only; test G)."
files="foo1.nc foo2.nc foo3.nc foo4.nc foo6.nc"

set -e -x

NRANGE="1 2 3 4 6"

# Create the files:
for NN in $NRANGE;
do 
    $MPIEXEC -n $NN $PISM_PATH/pismv -test G -Mx 40 -My 50 -Mz 60 -y 1 -verbose 1 -o foo$NN.nc 
done

set +e

# Compare:
for i in $NRANGE;
do
    for j in $NRANGE;
    do
	if [ $i -le $j ]; then continue; fi
	
	$PISM_PATH/nccmp.py -x -v rank foo$i.nc foo$j.nc
	if [ $? != 0 ];
	then
	    exit 1
	fi
    done
done

rm -f $files; exit 0
