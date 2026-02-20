#!/bin/bash

# Run the same simulation using synchronous and asynchronous output and save profiling
# data for a comparison.

pism_dir=${pism_dir:-$HOME/github/pism/pism}

pism_options="
-i ${pism_dir}/test/regression/pico_split/bedmap2_schmidtko14_50km.nc
-bootstrap
-grid.dx 10km
-grid.dy 10km
-Lz 5000
-Mz 11
-y 100
-output.snapshot.times 30days
-y 60days
-atmosphere uniform
-surface simple
-stress_balance sia
-verbose 2
-max_dt 5days
-spatial_times 15days
-spatial_vars thk,rank,mass_fluxes,temp
"
# Note: "rank" should be written *once* (it is time-independent) and "mass_fluxes" are not
# written to snapshot files, so this tests the ability to write disjoint sets of variables
# to different kinds of output files.

# run with asynchronous output:
time mpirun -n 7 pism ${pism_options} \
       -output.snapshot.file snapshots_async.nc \
       -output.spatial.file spatial_async.nc \
       -profile pism_async.py :\
       -n 1 python3 ${pism_dir}/util/pism_async_writer

# equivalent run using synchronous output:
time mpirun -n 7 pism ${pism_options} \
       -output.snapshot.file snapshots_sync.nc \
       -output.spatial.file spatial_sync.nc \
       -output.format netcdf4_parallel \
       -output.compression_level 1 \
       -profile pism_sync.py

# compare snapshot files:
excluded_vars=model_years_per_processor_hour,pism_config,wall_clock_time,mapping
pism_nccmp -x -v ${excluded_vars} snapshots_sync.nc snapshots_async.nc

pism_nccmp -x -v ${excluded_vars} spatial_sync.nc spatial_async.nc
