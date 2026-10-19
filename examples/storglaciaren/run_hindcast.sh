#!/bin/bash

set -x -e

ebm=debm_enhanced

# #############################################################################
# Initialization
# #############################################################################

climate_file=carra2_monthly_mean_RGI2000-v7.0-C-08-00658.nc

resolution=50m

outfile=state_g${resolution}_init.nc
mpiexec -n 2 pism  \
        -atmosphere.elevation_change.file ${climate_file}  \
        -atmosphere.elevation_change.periodic yes  \
        -atmosphere.elevation_change.temperature_lapse_rate 6.0  \
        -atmosphere.given.file ${climate_file}  \
        -atmosphere.given.periodic yes  \
        -atmosphere.models given,elevation_change  \
        -basal_resistance.pseudo_plastic.enabled yes  \
        -basal_resistance.pseudo_plastic.q 0.75  \
        -basal_resistance.pseudo_plastic.u_threshold 50m/yr  \
        -basal_yield_stress.model mohr_coulomb  \
        -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025  \
        -basal_yield_stress.mohr_coulomb.till_phi_default  30.0 \
        -bp_ksp_rtol 0.001  \
        -bp_ksp_type fgmres  \
        -bp_mg_coarse_ksp_max_it 50  \
        -bp_mg_coarse_ksp_rtol 0.01  \
        -bp_mg_coarse_ksp_type gmres  \
        -bp_mg_coarse_pc_type gamg  \
        -bp_mg_levels_ksp_type richardson  \
        -bp_mg_levels_pc_type sor  \
        -bp_pc_mg_levels 3  \
        -bp_pc_type mg  \
        -bp_snes_ksp_ew 1  \
        -bp_snes_ksp_ew_version 3  \
        -bp_snes_rtol 0.001  \
        -energy.model enthalpy \
        -geometry.front_retreat.use_cfl yes  \
        -geometry.part_grid.enabled yes  \
        -geometry.remove_icebergs yes  \
        -geometry.update.enabled yes  \
        -grid.Lbz 0  \
        -grid.Lz 250  \
        -grid.Mbz 1  \
        -grid.Mz 251  \
        -grid.dx ${resolution}  \
        -grid.dy ${resolution}  \
        -grid.file grid_RGI2000-v7.0-C-08-00658.nc  \
        -grid.registration center  \
        -hydrology.model null \
        -hydrology.null_diffuse_till_water yes \
        -hydrology.null_diffusion_distance 5e3 \
        -hydrology.null_diffusion_time 500 \
        -input.bootstrap yes  \
        -input.file pism_storglaciaren_1994_g10m.nc  \
        -input.forcing.buffer_size 24  \
        -input.forcing.time_extrapolation yes  \
        -output.compression_level 2  \
        -output.format netcdf4_serial  \
        -output.spatial.file spatial_g${resolution}_init.nc  \
        -output.spatial.stop_missing no  \
        -output.spatial.times 1  \
        -output.spatial.vars velbase_mag,velsurf_mag,tillwat  \
        -output.file $outfile \
        -output.compression_level 2  \
        -output.format netcdf4_serial  \
        -output.runtime.area_scale_factor_log10  0 \
        -output.size medium  \
        -output.sizes.medium sftgif,velsurf_mag,mask,usurf  \
        -stress_balance.blatter.Glen_exponent 3.0  \
        -stress_balance.blatter.Mz 10  \
        -stress_balance.blatter.coarsening_factor 3  \
        -stress_balance.blatter.enhancement_factor 1  \
        -stress_balance.blatter.flow_law gpbld  \
        -stress_balance.blatter.use_eta_transform yes  \
        -stress_balance.calving_front_stress_bc no  \
        -stress_balance.model blatter  \
        -surface.debm_simple.air_temp_all_precip_as_rain 275.15  \
        -surface.debm_simple.air_temp_all_precip_as_snow 273.15  \
        -surface.debm_simple.albedo_input.file ${climate_file}  \
        -surface.debm_simple.albedo_input.periodic yes  \
        -surface.debm_simple.c1 29  \
        -surface.debm_simple.c2 -93  \
        -surface.debm_simple.interpret_precip_as_snow no  \
        -surface.debm_simple.refreeze 0.6  \
        -surface.debm_simple.refreeze_ice_melt no  \
        -surface.debm_simple.std_dev 5.0  \
        -surface.models $ebm,forcing  \
        -surface.force_to_thickness.alpha 0.98  \
        -surface.force_to_thickness.file pism_forcing_mask_all_g10m.nc  \
        -surface.force_to_thickness.ice_free_alpha_factor 10  \
        -time.end 0011-01-01  \
        -time.calendar 365_day \
        -time.reference_date 0001-01-01  \
        -time.start 0001-01-01 \
        -time_stepping.adaptive_ratio 500  \
        -time_stepping.skip.enabled yes  \
        -time_stepping.skip.max 100 

exit

# #############################################################################
# Initialization
# #############################################################################

climate_file=carra2_RGI2000-v7.0-C-08-00658.nc

infile=$outfile
mpiexec -n 2 pism  \
        -atmosphere.elevation_change.file ${climate_file}  \
        -atmosphere.elevation_change.periodic yes  \
        -atmosphere.elevation_change.temperature_lapse_rate 6.0  \
        -atmosphere.given.file ${climate_file}  \
        -atmosphere.given.periodic yes  \
        -atmosphere.models given,elevation_change  \
        -basal_resistance.pseudo_plastic.enabled yes  \
        -basal_resistance.pseudo_plastic.q 0.75  \
        -basal_resistance.pseudo_plastic.u_threshold 10m/yr  \
        -basal_yield_stress.model mohr_coulomb  \
        -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025  \
        -basal_yield_stress.mohr_coulomb.till_phi_default 30  \
        -bp_ksp_rtol 0.001  \
        -bp_ksp_type fgmres  \
        -bp_mg_coarse_ksp_max_it 50  \
        -bp_mg_coarse_ksp_rtol 0.01  \
        -bp_mg_coarse_ksp_type gmres  \
        -bp_mg_coarse_pc_type gamg  \
        -bp_mg_levels_ksp_type richardson  \
        -bp_mg_levels_pc_type sor  \
        -bp_pc_mg_levels 3  \
        -bp_pc_type mg  \
        -bp_snes_ksp_ew 1  \
        -bp_snes_ksp_ew_version 3  \
        -bp_snes_rtol 0.001  \
        -energy.model enthalpy \
        -geometry.front_retreat.use_cfl yes  \
        -geometry.part_grid.enabled yes  \
        -geometry.remove_icebergs yes  \
        -geometry.update.enabled yes  \
        -grid.Lbz 0  \
        -grid.Lz 250  \
        -grid.Mbz 1  \
        -grid.Mz 251  \
        -grid.dx ${resolution}  \
        -grid.dy ${resolution}  \
        -grid.file grid_RGI2000-v7.0-C-08-00658.nc  \
        -input.bootstrap yes  \
        -input.file $infile \
        -input.forcing.buffer_size 24  \
        -input.regrid.file $infile  \
        -input.regrid.vars litho_temp,enthalpy,age,tillwat  \
        -input.forcing.time_extrapolation yes  \
        -output.compression_level 2  \
        -output.format netcdf4_serial  \
        -output.spatial.file spatial_g${resolution}_historical.nc  \
        -output.spatial.stop_missing no  \
        -output.spatial.times monthly  \
        -output.spatial.vars effective_air_temp,effective_precipitation,ice_surface_temp,climatic_mass_balance,pdd_fluxes,mass_fluxes,debm_insolation_driven_melt_flux,debm_insolation_driven_melt_rate,debm_temperature_driven_melt_flux,debm_temperature_driven_melt_rate,debm_offset_melt_flux,debm_offset_melt_rate,snow_depth,surface_albedo,atmosphere_transmissivity,insolation,velsurf_mag,velbase_mag,mask,usurf,ice_mass,thk,sftgif  \
        -output.file state_g${resolution}_historical.nc \
        -output.compression_level 2  \
        -output.format netcdf4_serial  \
        -output.size medium  \
        -output.sizes.medium sftgif,velsurf_mag,mask,usurf  \
        -stress_balance.blatter.Glen_exponent 3.0  \
        -stress_balance.blatter.Mz 10  \
        -stress_balance.blatter.coarsening_factor 3  \
        -stress_balance.blatter.enhancement_factor 2.0  \
        -stress_balance.blatter.flow_law gpbld  \
        -stress_balance.blatter.use_eta_transform yes  \
        -stress_balance.calving_front_stress_bc no  \
        -stress_balance.model blatter  \
        -surface.debm_simple.air_temp_all_precip_as_rain 275.15  \
        -surface.debm_simple.air_temp_all_precip_as_snow 273.15  \
        -surface.debm_simple.albedo_input.file ${climate_file}  \
        -surface.debm_simple.albedo_input.periodic yes  \
        -surface.debm_simple.c1 29  \
        -surface.debm_simple.c2 -93  \
        -surface.debm_simple.interpret_precip_as_snow no  \
        -surface.debm_simple.refreeze 0.6  \
        -surface.debm_simple.refreeze_ice_melt no  \
        -surface.debm_simple.std_dev 5.0  \
        -surface.models $ebm,forcing  \
        -surface.force_to_thickness.alpha 0.98  \
        -surface.force_to_thickness.file pism_forcing_mask_glacier_g10m.nc  \
        -surface.force_to_thickness.ice_free_alpha_factor 10  \
        -time.calendar proleptic_gregorian \
        -time.end 2025-01-01  \
        -time.reference_date 1986-01-01  \
        -time.start 1986-01-01 \
        -time_stepping.adaptive_ratio 500  \
        -time_stepping.skip.enabled yes  \
        -time_stepping.skip.max 100
        

