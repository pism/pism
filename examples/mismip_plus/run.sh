#!/bin/bash

ncgen -o config.nc config.cdl

spatial_vars="beta,bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line"
regrid_vars="litho_temp,enthalpy,tillwat,bmelt,ice_area_specific_volume,thk"

N=8
run_length=50000
sb="ssa+sia"
resolution="8km"
out=g${resolution}_${sb}_${run_length}a.nc
mpirun pism -config_override config.nc \
       -stress_balance.model $sb \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -input.file mismip+.nc \
       -input.bootstrap yes \
       -o_size medium \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 100year \
       -output.extra.vars $spatial_vars \
       -output.file state_$out \
       -time.run_length $run_length

infile=state_$out

N=8
run_length=1000
sb="ssa+sia"
resolution="4km"
out=g${resolution}_${sb}_${run_length}a.nc
mpirun pism -config_override config.nc \
       -stress_balance.model $sb \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -input.bootstrap yes \
       -input.file mismip+.nc \
       -input.regrid.file $infile \
       -input.regrid.vars $regrid_vars \
       -o_size medium \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 10year \
       -output.extra.vars $spatial_vars \
       -output.file state_$out \
       -time.run_length $run_length

N=8
run_length=500
sb="ssa+sia"
resolution="2km"
out=g${resolution}_${sb}_${run_length}a.nc
mpirun -np $N pism -config_override config.nc \
       -stress_balance.model $sb \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -i $infile \
       -o_size medium \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -time.run_length $run_length




exit




bootfile=boot_$out
cp state_$out $bootfile
ncrename -v u_ssa,ubar -v v_ssa,vbar -O boot_$out

N=8
run_length=1s
sb="prescribed_sliding"
resolution="2km"
out=pico_g${resolution}_${sb}_${run_length}a.nc
mpirun -np $N pism -config_override config.nc \
       -stress_balance.model $sb \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -i $bootfile \
       -stress_balance.prescribed_sliding.file $bootfile \
       -o_size medium \
       -ocean.models pico -ocean.pico.continental_shelf_depth -721 -ocean.pico.file ocean.nc -ocean.pico.heat_exchange_coefficent 2e-05 -ocean.pico.maximum_ice_rise_area 10000.0 -ocean.pico.number_of_boxes 5 -ocean.pico.overturning_coefficent 1000000.0 \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 1s \
       -output.extra.vars pico_temperature,pico_salinity,pico_basal_melt_rate,pico_salinity_box0,pico_temperature_box0,velbase,velsurf,bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line -output.file state_$out \
       -time.run_length $run_length

N=8
run_length=1s
sb="prescribed_sliding"
resolution="2km"
out=picop_g${resolution}_${sb}_${run_length}a.nc
mpirun -np $N pism -config_override config.nc \
       -stress_balance.model $sb \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -i $bootfile \
       -stress_balance.prescribed_sliding.file $bootfile \
       -o_size medium \
       -ocean.models picop -ocean.pico.continental_shelf_depth -721 -ocean.pico.file ocean.nc -ocean.pico.heat_exchange_coefficent 2e-05 -ocean.pico.maximum_ice_rise_area 10000.0 -ocean.pico.number_of_boxes 5 -ocean.pico.overturning_coefficent 1000000.0 \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 1s \
       -output.extra.vars picop_gammaTS,picop_length_scale,picop_geometric_scale,picop_temperature,picop_salinity,picop_basal_melt_rate,picop_grounding_line_elevation,picop_grounding_line_slope,picop_shelf_base_elevation,velbase,velsurf,bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line -output.file state_$out \
       -time.run_length $run_length



bootfile=boot_$out
cp state_$out  $bootfile
ncrename -v uvel,uvel_sigma -v vvel,vvel_sigma $bootfile


N=8
run_length=100
sb="blatter"
resolution="2km"
out=g${resolution}_${sb}_${run_length}a.nc
mpirun -np $N pism -config_override config.nc \
       -stress_balance.model $sb \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -i $bootfile \
       -bp_ksp_monitor -bp_ksp_view_singularvalues -bp_snes_monitor_ratio \
       -bp_pc_type mg \
       -bp_mg_levels_ksp_type richardson \
       -bp_mg_levels_pc_type sor \
       -bp_mg_coarse_ksp_type preonly \
       -bp_mg_coarse_pc_type lu \
       -bp_pc_mg_levels 3 \
       -bp_pc_type mg \
       -bp_snes_ksp_ew 1 \
       -bp_snes_ksp_ew_version 3 \
       -output.extra.file spatial_$out \
       -output.extra.times 1yr \
       -output.extra.vars bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line -output.file state_$out \
       -o_size medium \
       -output.sizes.medium sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -time.run_length $run_length


infile=state_$out
N=8
run_length=300
sb="ssa+sia"
resolution="2km"
out=g${resolution}_${sb}_${run_length}a.nc
mpirun -np $N pism -config_override config.nc \
       -stress_balance.model $sb \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -i $infile \
       -output.extra.file spatial_$out \
       -output.extra.times 1year \
       -output.extra.vars bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line -output.file state_$out \
       -output.sizes.medium sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -time.run_length $run_length

# Initialization 1km
mpirun -np 8 pism -basal_resistance.pseudo_plastic.enabled yes -basal_resistance.pseudo_plastic.q 0.3333333333333333 -basal_resistance.pseudo_plastic.u_threshold 100m/yr -basal_resistance.regularized_coulomb.enabled yes -basal_yield_stress.constant.value 10000000.0 -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 -bp_ksp_monitor  -bp_ksp_type gmres -bp_ksp_view_singularvalues  -bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type lu -bp_mg_levels_ksp_type richardson -bp_mg_levels_pc_type sor -bp_pc_mg_levels 3 -bp_pc_type mg -bp_snes_ksp_ew 1 -bp_snes_ksp_ew_version 3 -bp_snes_monitor_ratio  -constants.ice.density 918.0 -constants.sea_water.density 1028.0 -energy.model none -flow_law.isothermal_Glen.ice_softness 6.338e-25 -geometry.front_retreat.prescribed.file mismip+.nc -geometry.front_retreat.use_cfl  -geometry.part_grid.enabled  -geometry.remove_icebergs  -grid.Lbz 1000 -grid.Lz 6000 -grid.Mbz 11 -grid.Mz 601 -grid.dx 1km -grid.dy 1km -grid.registration corner -hydrology.model null -i state_g2km_5000a.nc -input.bootstrap yes -input.regrid.file state_g2km_5000a.nc -input.regrid.vars litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume,thk -ocean.constant.melt_rate 0.0 -ocean.models constant -output.extra.file spatial_g1km_200a.nc -output.extra.times 1year -output.extra.vars bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line -output.file state_g1km_200a.nc -output.sizes.medium sftgif,velsurf_mag,mask,usurf,bmelt,velbar  -stress_balance.calving_front_stress_bc  -stress_balance.model ssa -stress_balance.sia.flow_law isothermal_glen -stress_balance.sia.max_diffusivity 100000.0 -stress_balance.sia.surface_gradient_method eta -stress_balance.ssa.flow_law isothermal_glen -surface.given.file climate.nc -surface.models given -time.run_length 200 -time_stepping.adaptive_ratio 250 -time_stepping.skip.enabled  -time_stepping.skip.max 100

ncrename -v u_ssa,ubar -v v_ssa,vbar -O state_g1km_200a.nc state_g1km_200a.nc


exit

mkdir -p pico
# PICO run
mpirun -np 8 pism -basal_resistance.pseudo_plastic.enabled yes -basal_resistance.pseudo_plastic.q 0.3333333333333333 -basal_resistance.pseudo_plastic.u_threshold 100m/yr -basal_resistance.regularized_coulomb.enabled yes -basal_yield_stress.constant.value 10000000.0 -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 -bp_ksp_monitor  -bp_ksp_type gmres -bp_ksp_view_singularvalues  -bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type lu -bp_mg_levels_ksp_type richardson -bp_mg_levels_pc_type sor -bp_pc_mg_levels 3 -bp_pc_type mg -bp_snes_ksp_ew 1 -bp_snes_ksp_ew_version 3 -bp_snes_monitor_ratio  -constants.ice.density 918.0 -constants.sea_water.density 1028.0 -energy.model none -flow_law.isothermal_Glen.ice_softness 6.338e-25 -geometry.front_retreat.prescribed.file mismip+.nc -geometry.front_retreat.use_cfl  -geometry.part_grid.enabled  -geometry.remove_icebergs  -grid.Lbz 1000 -grid.Lz 6000 -grid.Mbz 11 -grid.Mz 601 -grid.dx 1km -grid.dy 1km -grid.registration corner -hydrology.model null -i state_g1km_200a.nc -input.bootstrap yes -input.regrid.file state_g1km_200a.nc -input.regrid.vars litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume,thk -ocean.constant.melt_rate 0.0 -ocean.models pico -ocean.pico.continental_shelf_depth -721 -ocean.pico.file ocean.nc -ocean.pico.heat_exchange_coefficent 2e-05 -ocean.pico.maximum_ice_rise_area 10000.0 -ocean.pico.number_of_boxes 5 -ocean.pico.overturning_coefficent 1000000.0 -output.extra.file pico/spatial_g1km.nc -output.extra.times 1s -output.extra.vars pico_basal_melt_rate,pico_salinity_box0,pico_temperature_box0,bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line -output.file pico/state_g1km_1a.nc -output.sizes.medium sftgif,velsurf_mag,mask,usurf,bmelt,velbar  -stress_balance.calving_front_stress_bc  -stress_balance.model prescribed_sliding -stress_balance.prescribed_sliding.file  state_g1km_200a_ubar.nc  -stress_balance.sia.flow_law isothermal_glen -stress_balance.sia.max_diffusivity 100000.0 -stress_balance.sia.surface_gradient_method eta -stress_balance.ssa.flow_law isothermal_glen -surface.given.file climate.nc -surface.models given -time.run_length 1s -time_stepping.adaptive_ratio 250 -time_stepping.skip.enabled  -time_stepping.skip.max 100 -ocean.pico.overturning_coefficent  0.23e6 -pico.heat_exchange_coefficent 1.28e-4 

mkdir -p picop
# PICOP run
mpirun -np 8 pism -basal_resistance.pseudo_plastic.enabled yes -basal_resistance.pseudo_plastic.q 0.3333333333333333 -basal_resistance.pseudo_plastic.u_threshold 100m/yr -basal_resistance.regularized_coulomb.enabled yes -basal_yield_stress.constant.value 10000000.0 -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 -bp_ksp_monitor  -bp_ksp_type gmres -bp_ksp_view_singularvalues  -bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type lu -bp_mg_levels_ksp_type richardson -bp_mg_levels_pc_type sor -bp_pc_mg_levels 3 -bp_pc_type mg -bp_snes_ksp_ew 1 -bp_snes_ksp_ew_version 3 -bp_snes_monitor_ratio  -constants.ice.density 918.0 -constants.sea_water.density 1028.0 -energy.model none -flow_law.isothermal_Glen.ice_softness 6.338e-25 -geometry.front_retreat.prescribed.file mismip+.nc -geometry.front_retreat.use_cfl  -geometry.part_grid.enabled  -geometry.remove_icebergs  -grid.Lbz 1000 -grid.Lz 6000 -grid.Mbz 11 -grid.Mz 601 -grid.dx 1km -grid.dy 1km -grid.registration corner -hydrology.model null -i state_g1km_200a.nc -input.bootstrap yes -input.regrid.file state_g1km_200a.nc -input.regrid.vars litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume,thk -ocean.constant.melt_rate 0.0 -ocean.models picop -ocean.pico.continental_shelf_depth -721 -ocean.pico.file ocean.nc -ocean.pico.heat_exchange_coefficent 2e-05 -ocean.pico.maximum_ice_rise_area 10000.0 -ocean.pico.number_of_boxes 5 -ocean.pico.overturning_coefficent 1000000.0 -output.extra.file picop/spatial_g1km.nc -output.extra.times 1s -output.extra.vars picop_gammaTS,picop_length_scale,picop_geometric_scale,picop_temperature,picop_salinity,picop_basal_melt_rate,picop_grounding_line_elevation,picop_grounding_line_slope,picop_shelf_base_elevation,bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line,picop_temperature,picop_salinity -output.file picop/state_g1km_1a.nc -output.sizes.medium sftgif,velsurf_mag,mask,usurf,bmelt,velbar -stress_balance.calving_front_stress_bc  -stress_balance.model prescribed_sliding -stress_balance.prescribed_sliding.file  state_g1km_200a_ubar.nc  -stress_balance.sia.flow_law isothermal_glen -stress_balance.sia.max_diffusivity 100000.0 -stress_balance.sia.surface_gradient_method eta -stress_balance.ssa.flow_law isothermal_glen -surface.given.file climate.nc -surface.models given -time.run_length 1s -time_stepping.adaptive_ratio 250 -time_stepping.skip.enabled  -time_stepping.skip.max 100 -picop.heat_exchange_parameter  1.69e-4 -ocean.pico.overturning_coefficent  0.23e6 -pico.heat_exchange_coefficent 1.28e-4


mpirun -np 8 pism -basal_resistance.pseudo_plastic.enabled yes -basal_resistance.pseudo_plastic.q 0.3333333333333333 -basal_resistance.pseudo_plastic.u_threshold 100m/yr -basal_resistance.regularized_coulomb.enabled yes -basal_yield_stress.constant.value 10000000.0 -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 -bp_ksp_monitor  -bp_ksp_type gmres -bp_ksp_view_singularvalues  -bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type lu -bp_mg_levels_ksp_type richardson -bp_mg_levels_pc_type sor -bp_pc_mg_levels 3 -bp_pc_type mg -bp_snes_ksp_ew 1 -bp_snes_ksp_ew_version 3 -bp_snes_monitor_ratio  -constants.ice.density 918.0 -constants.sea_water.density 1028.0 -energy.model none -flow_law.isothermal_Glen.ice_softness 6.338e-25 -geometry.front_retreat.prescribed.file mismip+.nc -geometry.front_retreat.use_cfl  -geometry.part_grid.enabled  -geometry.remove_icebergs  -grid.Lbz 1000 -grid.Lz 6000 -grid.Mbz 11 -grid.Mz 601 -grid.dx 1km -grid.dy 1km -grid.registration corner -hydrology.model null -i state_g1km_200a.nc -input.bootstrap yes -input.regrid.file state_g1km_200a.nc -input.regrid.vars litho_temp,enthalpy,age,tillwat,bmelt,ice_area_specific_volume,thk -ocean.constant.melt_rate 0.0 -ocean.models constant -output.extra.file spatial_ho_g1km_1a.nc -output.extra.times monthly -output.extra.vars bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line -output.file state_ho_g1km_1a.nc -output.sizes.medium sftgif,velsurf_mag,mask,usurf,bmelt,velbar  -stress_balance.calving_front_stress_bc  -stress_balance.model blatter -stress_balance.sia.flow_law isothermal_glen -stress_balance.sia.max_diffusivity 100000.0 -stress_balance.sia.surface_gradient_method eta -stress_balance.ssa.flow_law isothermal_glen -surface.given.file climate.nc -surface.models given -time.run_length 1 -time_stepping.adaptive_ratio 250 -time_stepping.skip.enabled  -time_stepping.skip.max 100
