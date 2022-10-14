#!/bin/bash

NN=4

INFILE=startSMALLablate.nc
GRID="-Lz 3500 -Mx 51 -My 51 -Mz 61 -y 500"
OPTS="-no_temp -spatial_vars thk,usurf,h_x,h_y -spatial_times 0:10:500  -o_size big -o_order zyx"

# With the bed smoother (default)
for method in mahaffy eta haseloff;
do
    mpiexec -n $NN pismr -i $INFILE -bootstrap $GRID -gradient $method $OPTS -spatial_file ex-$method.nc -o o-$method.nc
done;

# Without the bed smoother
OPTS="$OPTS -bed_smoother_range 0.0"
for method in mahaffy haseloff;
do
    mpiexec -n $NN pismr -i $INFILE -bootstrap $GRID -gradient $method $OPTS -spatial_file ex-no-smoother-$method.nc -o o-no-smoother-$method.nc
done;
