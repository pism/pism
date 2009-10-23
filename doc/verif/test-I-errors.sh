#!/bin/bash

# this script runs pismv to produce a figure for the User's Manual.

N=2

rm -i test-I-errors.nc
for Mx in 51 101 151 201 401 601 801 1001 1501 2001 2501 3073
do
    set +e
    command="mpiexec -np $N pismv -test I -Mx $Mx -My 5 -ssa_rtol 5e-07 -ksp_rtol 1e-12 -report_file test-I-errors.nc -verbose 1"
    echo "running $command"
    $command
done
