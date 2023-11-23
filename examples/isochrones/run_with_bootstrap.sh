#!/bin/bash

# Runs a flat-bed isothermal SIA setup using EISMINT-II experiment A surface mass balance,
# bootstrapping from the output of a 5000-year long run and creating 5 "bootstrapping"
# isochrone layers. Adds new isochrone layers every 500 years.

# number of MPI processes to use
N=${N:-8}

# horizontal grid size
M=${M:-100}

common_options="
        -stress_balance.sia.flow_law isothermal_glen
        -stress_balance.sia.surface_gradient_method eta
        -energy.enabled no
"

mpiexec -n ${N} pismr -eisII A \
        -grid.Mx ${M} \
        -grid.My ${M} \
        -output.file input.nc \
        -time.end 5e3 \
        ${common_options}

mpiexec -n ${N} pismr \
        -i input.nc \
        -bootstrap \
        -grid.Lz 3500 \
        -isochrones.deposition_times 0:500:20e3 \
        -isochrones.bootstrapping.n_layers 5 \
        -output.extra.file ex_bootstrap.nc \
        -output.file o_bootstrap.nc \
        -output.extra.times 50 \
        -output.extra.vars isochrone_depth,thk \
        -time.start 0 \
        -time.end 20e3 \
        ${common_options}
