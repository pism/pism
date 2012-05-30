#!/bin/bash

NN=2

INFILE=startSMALLablate.nc

OPTS="-Lz 3500 -Mx 51 -My 51 -Mz 61 -no_temp -y 500 -extra_vars thk,csurf -extra_times 0:10:500 -ts_times 0:1:500"

mpiexec -n $NN pismr -boot_file $INFILE -gradient haseloff $OPTS -extra_file Hex.nc -ts_file Hts.nc -o H.nc

mpiexec -n $NN pismr -boot_file $INFILE -gradient mahaffy $OPTS -extra_file Mex.nc -ts_file Mts.nc -o M.nc

mpiexec -n $NN pismr -boot_file $INFILE -gradient eta $OPTS -extra_file Eex.nc -ts_file Ets.nc -o E.nc

