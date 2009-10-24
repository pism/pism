#!/bin/bash

# this script runs pismv to produce a figure for the User's Manual.

N=8

rm -i test-B-errors.nc
for Mx in 31 41 51 61 71 81 91 101 111 121
do
    set +e
    command="mpiexec -np $N pismv -test B -Mx $Mx -My $Mx -Mz 31 -ys 422.45 -y 25000.0 -report_file test-B-errors.nc -verbose 1"
    echo "running $command"
    $command
done
