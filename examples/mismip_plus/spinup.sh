#!/bin/bash

# stop on error
set -e

NN=8
if [ $# -gt 0 ] ; then
  NN="$1" # first arg is number of processes
  echo "Using $NN cores"
fi

python3 preprocess.py
ncgen -o config.nc config.cdl

boot_file=mismip+.nc
climate_file=climate.nc
ocean_file=ocean.nc

spatial_vars="beta,bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag"
pico_spatial_vars="pico_box_mask,pico_shelf_mask,pico_contshelf_mask,pico_ice_rise_mask,pico_T_star,pico_overturning,pico_temperature,pico_salinity,pico_basal_melt_rate"
regrid_vars="litho_temp,enthalpy,tillwat,bmelt,ice_area_specific_volume,thk"

run_length=5000yr
sb="ssa+sia"
resolution="8km"
out=g${resolution}_${sb}_${run_length}.nc
mpirun -np $NN pism -config_override config.nc \
       -geometry.front_retreat.prescribed.file $boot_file \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -input.file $boot_file \
       -input.bootstrap yes \
       -o_size medium \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 100year \
       -output.extra.vars $spatial_vars \
       -output.file state_$out \
       -stress_balance.model $sb \
       -surface.given.file $climate_file \
       -time.run_length $run_length

infile=state_$out

run_length=1000yr
sb="ssa+sia"
resolution="4km"
out=g${resolution}_${sb}_${run_length}.nc
mpirun -np $NN pism -config_override config.nc \
       -stress_balance.model $sb \
       -geometry.front_retreat.prescribed.file $boot_file \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -input.bootstrap yes \
       -input.file $boot_file \
       -input.regrid.file $infile \
       -input.regrid.vars $regrid_vars \
       -o_size medium \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 10year \
       -output.extra.vars $spatial_vars \
       -output.file state_$out \
       -surface.given.file $climate_file \
       -time.run_length $run_length

infile=state_$out

run_length=100yr
sb="ssa+sia"
resolution="2km"
out=g${resolution}_${sb}_${run_length}.nc
mpirun -np $NN pism -config_override config.nc \
       -stress_balance.model $sb \
       -geometry.front_retreat.prescribed.file $boot_file \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -input.bootstrap yes \
       -input.file $boot_file \
       -input.regrid.file $infile \
       -input.regrid.vars $regrid_vars \
       -o_size medium \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 10year \
       -output.extra.vars $spatial_vars \
       -output.file state_$out \
       -surface.given.file $climate_file \
       -time.run_length $run_length

export infile=state_$out
