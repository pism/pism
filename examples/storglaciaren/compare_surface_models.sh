#!/bin/bash

set -x -e

for ebm in pdd debm_simple debm_enhanced; do
    mpiexec -n 1  pism_async_writer -f : -n  6 pism  \
            -atmosphere.elevation_change.file carra2_RGI2000-v7.0-C-08-00658.nc  \
            -atmosphere.elevation_change.periodic yes  \
            -atmosphere.elevation_change.temperature_lapse_rate 6.0  \
            -atmosphere.given.file carra2_RGI2000-v7.0-C-08-00658.nc  \
            -atmosphere.given.periodic yes  \
            -atmosphere.models given,elevation_change  \
            -energy.model none \
            -geometry.front_retreat.use_cfl yes  \
            -geometry.part_grid.enabled yes  \
            -geometry.remove_icebergs yes  \
            -geometry.update.enabled yes  \
            -grid.Lbz 0  \
            -grid.Lz 250  \
            -grid.Mbz 1  \
            -grid.Mz 251  \
            -grid.dx 10m  \
            -grid.dy 10m  \
            -grid.file grid_RGI2000-v7.0-C-08-00658.nc  \
            -grid.registration center  \
            -input.bootstrap yes  \
            -input.file pism_storglaciaren_g10m.nc  \
            -input.forcing.buffer_size 24  \
            -input.forcing.time_extrapolation yes  \
            -output.compression_level 2  \
            -output.format netcdf4_serial  \
            -output.asynchronous \
            -output.spatial.file spatial_g10m_${ebm}.nc  \
            -output.spatial.stop_missing no  \
            -output.spatial.times daily  \
            -output.spatial.vars effective_air_temp,effective_precipitation,ice_surface_temp,climatic_mass_balance,pdd_fluxes,mass_fluxes,debm_insolation_driven_melt_flux,debm_insolation_driven_melt_rate,debm_temperature_driven_melt_flux,debm_temperature_driven_melt_rate,debm_offset_melt_flux,debm_offset_melt_rate,snow_depth,surface_albedo,atmosphere_transmissivity,insolation  \
            -stress_balance.model none  \
            -surface.debm_simple.air_temp_all_precip_as_rain 275.15  \
            -surface.debm_simple.air_temp_all_precip_as_snow 273.15  \
            -surface.debm_simple.albedo_input.file carra2_RGI2000-v7.0-C-08-00658.nc  \
            -surface.debm_simple.albedo_input.periodic yes  \
            -surface.debm_simple.c1 29  \
            -surface.debm_simple.c2 -93  \
            -surface.debm_simple.interpret_precip_as_snow no  \
            -surface.debm_simple.refreeze 0.6  \
            -surface.debm_simple.refreeze_ice_melt no  \
            -surface.debm_simple.std_dev 5.0  \
            -surface.models $ebm  \
            -surface.pdd.factor_ice 0.00879120879120879 \
            -surface.pdd.factor_snow 0.0032967032967033 \
            -surface.pdd.refreeze 0.6  \
            -surface.pdd.refreeze_ice_melt no  \
            -surface.pdd.std_dev.value 5.0  \
            -time.calendar proleptic_gregorian  \
            -time.end 1995-01-01  \
            -time.reference_date 1980-01-01  \
            -time.start 1994-01-01 

  cdo -O timsum -setattribute,surface_melt_flux@units="m" -divdpy -divc,1000 -selvar,surface_melt_flux -seldate,1994-07-04,1994-08-25 spatial_g10m_${ebm}.nc spatial_g10m_${ebm}_SUM.nc
  cdo -O yearmean spatial_g10m_${ebm}.nc spatial_g10m_${ebm}_YM.nc
done
