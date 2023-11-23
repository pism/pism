#!/bin/bash

# Run PISM for 50000 years to approximate the steady state.

variables=mass_fluxes,thk,usurf,diffusivity,diffusivity_staggered,flux_mag,wvelsurf,mask,h_x

mpiexec -n 4 \
pismr \
  -bootstrap \
  -constants.ice.density 910 \
  -constants.standard_gravity 9.81 \
  -energy none \
  -flow_law.isothermal_Glen.ice_softness 1e-16.Pa-3.year-1 \
  -grid.Lz 400 \
  -grid.Mz 41 \
  -i input.nc \
  -output.extra.file ex.nc \
  -output.extra.times 1000 \
  -output.extra.vars ${variables} \
  -output.file bedrock_step.nc \
  -output.timeseries.filename ts.nc \
  -output.timeseries.times 1 \
  -stress_balance.ice_free_thickness_standard 0.01 \
  -stress_balance.model sia \
  -stress_balance.sia.Glen_exponent 3 \
  -stress_balance.sia.bed_smoother.range 0 \
  -stress_balance.sia.flow_law isothermal_glen \
  -y 50000 \
  > bedrock_step.log \
  ;
