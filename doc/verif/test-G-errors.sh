#!/bin/bash

# this script runs pismv to produce a figure for the User's Manual.

N=8

rm -i test-G-errors.nc
for Mx in 61 71 81 91 101 111 121 151 181
do
    set +e
    command="mpiexec -np $N pismv -test G -Mx $Mx -My $Mx -Mz $Mx -y 25000.0 -report_file test-G-errors.nc -verbose 1"
    echo "running $command"
    $command
done
