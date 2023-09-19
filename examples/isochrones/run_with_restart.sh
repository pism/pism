#!/bin/bash
# Runs a flat-bed isothermal SIA setup using EISMINT-II experiment A surface mass balance,
# starting from zero ice thickness. Adds new isochrones every 1000 years.

# number of MPI processes to use
N=${N:-8}

# horizontal grid size
M=${M:-101}

common_options="
        -output.extra.file ex_with_restart.nc
        -output.extra.times 50
        -output.extra.vars isochrone_depth,thk
        -output.sizes.medium isochrone_depth,uvel
        -stress_balance.sia.flow_law isothermal_glen
        -stress_balance.sia.surface_gradient_method eta
        -isochrones.deposition_times 1000
        -energy.enabled no
"

mpiexec -n ${N} pismr -eisII A \
        -bootstrapping.defaults.geothermal_flux 0 \
        -grid.Lz 3500 \
        -grid.Mx ${M} \
        -grid.My ${M} \
        -grid.Mz 21 \
        -isochrones.bootstrapping.n_layers 0 \
        -output.file intermediate_result.nc \
        -time.end 5e3 \
        ${common_options}

mpiexec -n ${N} pismr \
        -i intermediate_result.nc \
        -output.extra.append \
        -output.file final_result.nc \
        -time.end 10e3 \
        ${common_options}
