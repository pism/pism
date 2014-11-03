#!/bin/bash

./util/FEvoRslab.py 

./util/flowline.py -o fevor-slab-in.nc --expand -d y fevor-slab.nc

mpiexec -n 2 pismr -surface given -boot_file fevor-slab-in.nc -periodicity xy \
 -Mx 201 -My 3 -Lx 10.05 -Mz 11 -Lz 1100 -y 1 \
 -bed_smoother_range 0 \
 -o pism-out.nc

/util/flowline.py -o fevor-slab-out.nc --collapse -d y pism-out.nc
