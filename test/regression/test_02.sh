#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

echo "Test # 2: pismv exact processor independence (SIA only; test G)."
files="foo1-02.nc foo2-02.nc foo3-02.nc foo4-02.nc foo6-02.nc"

set -e -x

NRANGE="1 2 3 4 6"

# Create the files:
for NN in $NRANGE;
do
    $MPIEXEC -n $NN $PISM_PATH/pismv -test G -Mx 30 -My 40 -Mz 20 -y 1 -verbose 1 -o foo$NN-02.nc
done

set +e

# Compare:
for i in $NRANGE;
do
    for j in $NRANGE;
    do
	if [ $i -le $j ]; then continue; fi

	$PISM_PATH/nccmp.py -x -v rank,timestamp foo$i-02.nc foo$j-02.nc
	if [ $? != 0 ];
	then
	    exit 1
	fi
    done
done

rm -f $files; exit 0
