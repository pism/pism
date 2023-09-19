#!/bin/bash
# Runs a flat-bed isothermal SIA setup using EISMINT-II experiment A surface mass balance,
# starting from zero ice thickness. Adds new isochrones every 1000 years.

# number of MPI processes to use
N=${N:-8}

# horizontal grid size
M=${M:-101}

mpiexec -n ${N} pismr -eisII A \
        -bootstrapping.defaults.geothermal_flux 0 \
        -energy.enabled no \
        -grid.Lz 3500 \
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
        -stress_balance.sia.flow_law isothermal_glen \
        -stress_balance.sia.surface_gradient_method eta \
        -time.run_length 20e3 \
  ;
