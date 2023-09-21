#!/bin/bash
# Runs EISMINT-II experiment A, adding new isochrones every 1000 years.

# number of MPI processes to use
N=${N:-8}

# horizontal grid size
M=${M:-201}

mpiexec -n ${N} pismr -eisII A \
        -grid.Lz 5000 \
        -grid.Mx ${M} \
        -grid.My ${M} \
        -grid.Mz 21 \
        -isochrones.deposition_times 1000 \
        -isochrones.bootstrapping.n_layers 0 \
        -output.extra.file ex.nc \
        -output.extra.times 50 \
        -output.extra.vars isochrone_depth,thk \
        -output.file output.nc \
        -output.sizes.medium isochrone_depth,uvel \
        -time.run_length 20e3
