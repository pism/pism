#!/bin/bash

# Run the same simulation using synchronous and asynchronous output and save profiling
# data for a comparison.

pism_dir=${pism_dir:-$HOME/github/pism/pism}

pism_options="
-i ${pism_dir}/test/regression/pico_split/bedmap2_schmidtko14_50km.nc
-bootstrap
-bed_def lc
-grid.dx 15km
-grid.dy 15km
-Lz 5000
-Mz 11
-y 100
-output.snapshot.times 10days
-y 60days
-atmosphere uniform
-surface simple
-stress_balance sia
-verbose 2
-max_dt 50
"

# run with asynchronous output:
mpirun -n 7 pism ${pism_options} -output.snapshot.file snapshots_async.nc -profile pism_async.py :\
       -n 1 python3 ${pism_dir}/util/pism_async_writer -d

# equivalent run using synchronous output:
mpirun -n 7 pism ${pism_options} -output.snapshot.file snapshots_sync.nc -profile pism_sync.py

# compare snapshot files:
excluded_vars=model_years_per_processor_hour,pism_config,wall_clock_time,mapping
pism_nccmp -x -v ${excluded_vars} snapshots_sync.nc snapshots_async.nc
