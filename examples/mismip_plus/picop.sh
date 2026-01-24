#!/bin/bash

# stop on error
set -e

NN=8
if [ $# -gt 0 ] ; then
  NN="$1" # first arg is number of processes
  echo "Using $NN cores"
fi

infile=state_g2km_ssa+sia_100yr.nc
if [ $# -gt 1 ] ; then
  infile="$2" 
  echo "Using state file $infile"
fi

boot_file=mismip+.nc
climate_file=climate.nc
ocean_file=ocean.nc

spatial_vars="beta,bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag"
picop_spatial_vars="picop_basal_melt_rate,picop_grounding_line_slope,picop_grounding_line_elevation"
regrid_vars="litho_temp,enthalpy,tillwat,bmelt,ice_area_specific_volume,thk"

run_length=1s
sb="ssa+sia"
resolution="2km"
out=picop_g${resolution}_${sb}_${run_length}.nc
mpirun -np $NN pism -config_override config.nc \
       -stress_balance.model $sb \
       -geometry.front_retreat.prescribed.file $boot_file \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -input.file $infile \
       -ocean.models picop \
       -ocean.pico.file $ocean_file \
       -o_size medium \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 1s \
       -output.extra.vars $picop_spatial_vars,$spatial_vars \
       -output.file state_$out \
       -surface.given.file $climate_file \
       -time.start 0 \
       -time.run_length $run_length

run_length=100yr
sb="ssa+sia"
resolution="2km"
out=picop_g${resolution}_${sb}_${run_length}.nc
mpirun -np $NN pism -config_override config.nc \
       -stress_balance.model $sb \
       -geometry.front_retreat.prescribed.file $boot_file \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -input.file $infile \
       -ocean.models picop \
       -ocean.pico.file $ocean_file \
       -o_size medium \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 1yr \
       -output.extra.vars $picop_spatial_vars,$spatial_vars \
       -output.file state_$out \
       -surface.given.file $climate_file \
       -time.start 0 \
       -time.run_length $run_length
