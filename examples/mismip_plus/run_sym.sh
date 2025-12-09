#!/bin/bash

# stop on error
set -e

python3 preprocess_sym.py
ncgen -o config.nc config.cdl

boot_file=mismip+_sym.nc
spatial_vars="beta,bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line,tendency_of_ice_mass_due_to_forced_retreat"
pico_spatial_vars="pico_box_mask,pico_shelf_mask,pico_contshelf_mask,pico_ice_rise_mask,pico_T_star,pico_overturning,pico_temperature,pico_salinity,pico_basal_melt_rate"
regrid_vars="litho_temp,enthalpy,tillwat,bmelt,ice_area_specific_volume,thk"

run_length=5000
sb="ssa+sia"
resolution="8km"
out=g${resolution}_${sb}_${run_length}a.nc
# mpirun pism -config_override config.nc \
#        -geometry.front_retreat.prescribed.file $boot_file \
#        -grid.dx $resolution \
#        -grid.dy $resolution \
#        -input.file $boot_file \
#        -input.bootstrap yes \
#        -o_size medium \
#        -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
#        -output.extra.file spatial_$out \
#        -output.extra.times 100year \
#        -output.extra.vars $spatial_vars \
#        -output.file state_$out \
#        -stress_balance.model $sb \
#        -surface.given.file climate_sym.nc \
#        -time.run_length $run_length

infile=state_$out

run_length=1000
sb="ssa+sia"
resolution="4km"
out=g${resolution}_${sb}_${run_length}a.nc
# mpirun pism -config_override config.nc \
#        -stress_balance.model $sb \
#        -geometry.front_retreat.prescribed.file $boot_file \
#        -grid.dx $resolution \
#        -grid.dy $resolution \
#        -input.bootstrap yes \
#        -input.file $boot_file \
#        -input.regrid.file $infile \
#        -input.regrid.vars $regrid_vars \
#        -o_size medium \
#        -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
#        -output.extra.file spatial_$out \
#        -output.extra.times 10year \
#        -output.extra.vars $spatial_vars \
#        -output.file state_$out \
#        -surface.given.file climate_sym.nc \
#        -time.run_length $run_length

infile=state_$out

run_length=100
sb="ssa+sia"
resolution="2km"
out=g${resolution}_${sb}_${run_length}a.nc
# mpirun pism -config_override config.nc \
#        -stress_balance.model $sb \
#        -geometry.front_retreat.prescribed.file $boot_file \
#        -grid.dx $resolution \
#        -grid.dy $resolution \
#        -input.bootstrap yes \
#        -input.file $boot_file \
#        -input.regrid.file $infile \
#        -input.regrid.vars $regrid_vars \
#        -o_size medium \
#        -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
#        -output.extra.file spatial_$out \
#        -output.extra.times 10year \
#        -output.extra.vars $spatial_vars \
#        -output.file state_$out \
#        -surface.given.file climate_sym.nc \
#        -time.run_length $run_length
>>>>>>> Stashed changes

infile=state_$out

N=8
run_length=1s
sb="ssa+sia"
resolution="2km"
out=pico_g${resolution}_${sb}_${run_length}a.nc
# pism -config_override config.nc \
#        -stress_balance.model $sb \
#        -geometry.front_retreat.prescribed.file $boot_file \
#        -grid.dx $resolution \
#        -grid.dy $resolution \
#        -input.bootstrap yes \
#        -input.file $boot_file \
#        -input.regrid.file $infile \
#        -input.regrid.vars $regrid_vars \
#        -ocean.models pico \
#        -ocean.pico.continental_shelf_depth -721 \
#        -ocean.pico.file ocean_sym.nc \
#        -ocean.pico.heat_exchange_coefficent 2e-05 \
#        -ocean.pico.maximum_ice_rise_area 10000.0 \
#        -ocean.pico.number_of_boxes 5 \
#        -ocean.pico.overturning_coefficent 1000000.0 \
#        -o_size medium \
#        -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
#        -output.extra.file spatial_$out \
#        -output.extra.times 1s \
#        -output.extra.vars $pico_spatial_vars,$spatial_vars \
#        -output.file state_$out \
#        -surface.given.file climate_sym.nc \
#        -time.run_length $run_length


run_length=1s
sb="ssa+sia"
resolution="2km"
out=picop_g${resolution}_${sb}_${run_length}a.nc
mpirun -np $N pism -config_override config.nc \
       -stress_balance.model $sb \
       -geometry.front_retreat.prescribed.file $boot_file \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -input.bootstrap yes \
       -input.file $boot_file \
       -input.regrid.file $infile \
       -input.regrid.vars $regrid_vars \
       -ocean.models picop \
       -ocean.pico.continental_shelf_depth -721 \
       -ocean.pico.file ocean_sym.nc \
       -ocean.pico.heat_exchange_coefficent 2e-05 \
       -ocean.pico.maximum_ice_rise_area 10000.0 \
       -ocean.pico.number_of_boxes 5 \
       -ocean.pico.overturning_coefficent 1000000.0 \
       -o_size medium \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 1s \
       -output.extra.vars picop_gammaTS,picop_length_scale,picop_geometric_scale,picop_temperature,picop_salinity,picop_basal_melt_rate,picop_grounding_line_elevation,picop_grounding_line_slope,picop_shelf_base_elevation,$spatial_vars \
       -output.file state_$out \
       -surface.given.file climate_sym.nc \
       -time.run_length $run_length


