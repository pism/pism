#!/bin/bash

# this script runs pismv to produce a figure for the User's Manual.

rm -i test-K-errors.nc
for Mz in 21 41 61 81 101 121 141 161 181 201 221 241 261 281 301 321
do
    set +e
    Mbz=$(( ($Mz - 1) / 4 + 1))
    command="mpiexec -n 2 pismv -test K -Mz $Mz -Mbz $Mbz -Lbz 1000 -Mx 4 -My 4 -y 130000.0 -report_file test-K-errors.nc -verbose 1"
    echo "running $command"
    $command
done
