#!/bin/bash

NN=4

INFILE=startSMALLablate.nc
GRID="-Lz 3500 -Mx 51 -My 51 -Mz 61 -y 500"
OPTS="-no_temp -extra_vars thk,usurf,h_x,h_y -extra_times 0:10:500  -o_size big -o_order zyx"

# With the bed smoother (default)
for method in mahaffy eta haseloff new;
do
    mpiexec -n $NN pismr -boot_file $INFILE $GRID -gradient $method $OPTS -extra_file ex-$method.nc -o o-$method.nc
done;

# Without the bed smoother
OPTS="$OPTS -config_override ./no_bed_smoother.nc"
for method in mahaffy haseloff new;
do
    mpiexec -n $NN pismr -boot_file $INFILE $GRID -gradient $method $OPTS -extra_file ex-ns-$method.nc -o o-ns-$method.nc
done;
