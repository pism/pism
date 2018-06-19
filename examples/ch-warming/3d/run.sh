#!/bin/bash

python create_inputs.py

run() {
suffix=$1
warming=$2
extra_vars=$3

# bootstrap ice temperature (and so its enthalpy) using the mean-annual surface
# temperature
pismr -bootstrap -i input.nc \
      -regional \
      -regional.no_model_strip 0 \
      -Lz 200 -Mz 201 -grid.ice_vertical_spacing equal \
      -surface given \
      -y 0 \
      -o in.nc -verbose 1

echo "running PISM (${suffix})..."

# bootstrap again and bring in the enthalpy field computed by the run above
pismr -bootstrap -i input.nc -regrid_file in.nc -regrid_vars enthalpy \
      -bootstrapping.defaults.geothermal_flux 0 \
      -regional \
      -regional.no_model_strip 0 \
      -Lz 200 -Mz 201 -grid.ice_vertical_spacing equal \
      -energy.ch_warming.enabled ${warming} \
      -energy.ch_warming.average_channel_spacing 20 \
      -energy.ch_warming.residual_water_fraction 0.005 \
      -energy.ch_warming.temperate_ice_thermal_conductivity_ratio 1.0 \
      -surface given,delta_T \
      -surface_delta_T_file input.nc -surface_delta_T_period 1 \
      -extra_file ex_${suffix}.nc \
      -extra_vars ${extra_vars} \
      -extra_times 10days \
      -y 10 \
      -calendar 360_day \
      -o o_${suffix}.nc -verbose 1

rm in.nc
}

variables=temp,temp_pa,ice_surface_temp

# run without cryo-hydrologic warming
run no_warming False ${variables}

# run with cryo-hydrologic warming
run warming True ${variables},ch_heat_flux,ch_temp,ch_liqfrac

ncpdq -a z,time,y,x -O ex_no_warming.nc ex_no_warming.nc
ncpdq -a z,time,y,x -O ex_warming.nc ex_warming.nc

# compare temperature fields
ncdiff -O -v temp ex_warming.nc ex_no_warming.nc temp_difference.nc
