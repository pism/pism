#!/bin/bash

# Run the same simulation using synchronous and asynchronous output and save profiling
# data for a comparison.

pism_dir=${pism_dir:-$HOME/github/pism/pism}

pism_options="
-i ${pism_dir}/test/regression/pico_split/bedmap2_schmidtko14_50km.nc
-bootstrap
-grid.dx 5km
-grid.dy 5km
-Lz 5000
-Mz 11
-y 100
-output.snapshot.file snapshots.nc
-output.snapshot.times 10days
-y 60days
-atmosphere uniform
-surface simple
-stress_balance sia
-verbose 2
-max_dt 50
"

mpirun -n 7 pism ${pism_options} -profile pism_async.py :\
       -n 1 python3 ${pism_dir}/util/pism_async_writer -d

mpirun -n 7 pism ${pism_options} -profile pism_sync.py
