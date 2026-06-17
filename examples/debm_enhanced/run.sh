#!/bin/bash

for f in bootfile_RGI2000-v7.0-C-01-04374.nc era5_wgs84_RGI2000-v7.0-C-01-04374.nc insolation_RGI2000-v7.0-C-01-04374.nc; do
    wget -nc https://pism-cloud-data.s3.amazonaws.com/debm_enhanced/$f;
done

    
mpirun -np 8 pism -i bootfile_RGI2000-v7.0-C-01-04374.nc \
       -input.bootstrap yes  \
       -grid.dx 500m \
       -grid.dy 500m \
       -energy.model none \
       -stress_balance.model none \
       -atmosphere.models given,elevation_change \
       -atmosphere.elevation_change.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc \
       -atmosphere.elevation_change.temperature_lapse_rate 6 \
       -atmosphere.given.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc \
       -atmosphere.given.periodic yes \
       -surface pdd \
       -time.start 2000-01-01 \
       -time.end 2001-01-01 \
       -time.calendar standard \
       -output.spatial.file pdd.nc \
       -output.spatial.vars mass_fluxes,pdd_fluxes,effective_air_temp,effective_precipitation,climatic_mass_balance \
       -output.spatial.times monthly

mpirun -np 8 pism -i bootfile_RGI2000-v7.0-C-01-04374.nc \
       -input.bootstrap yes  \
       -grid.dx 500m \
       -grid.dy 500m \
       -energy.model none \
       -stress_balance.model none \
       -atmosphere.models given,elevation_change \
       -atmosphere.elevation_change.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc \
       -atmosphere.elevation_change.temperature_lapse_rate 6 \
       -atmosphere.given.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc \
       -atmosphere.given.periodic yes \
       -surface debm_simple \
       -surface.debm_simple.std_dev.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc \
       -surface.debm_simple.std_dev.periodic yes \
       -surface.debm_simple.interpret_precip_as_snow no \
       -surface.debm_simple.c1 25 \
       -surface.debm_simple.c2 -103 \
       -surface.debm_simple.refreeze  0.50 \
       -surface.debm_simple.refreeze_ice_melt no \
       -surface.debm_simple.air_temp_all_precip_as_snow 271.65 \
       -surface.debm_simple.air_temp_all_precip_as_rain 275.0\
       -time.start 2000-01-01 \
       -time.end 2001-01-01 \
       -time.calendar standard \
       -output.spatial.file debm_simple.nc \
       -output.spatial.vars mass_fluxes,pdd_fluxes,effective_air_temp,effective_precipitation,climatic_mass_balance \
       -output.spatial.times monthly


mpirun -np 8 pism -i bootfile_RGI2000-v7.0-C-01-04374.nc \
       -energy.model none \
       -grid.dx 500m \
       -grid.dy 500m \
       -input.bootstrap yes  \
       -input.forcing.time_extrapolation yes \
       -stress_balance.model none \
       -atmosphere.models given \
       -atmosphere.given.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc \
       -atmosphere.given.periodic yes \
       -surface debm_enhanced \
       -surface.debm_enhanced.file insolation_RGI2000-v7.0-C-01-04374.nc \
       -surface.debm_enhanced.periodic yes  \
       -surface.debm_simple.std_dev.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc \
       -surface.debm_simple.std_dev.periodic yes \
       -surface.debm_simple.interpret_precip_as_snow no \
       -surface.debm_simple.c1 25 \
       -surface.debm_simple.c2 -103 \
       -surface.debm_simple.refreeze  0.50 \
       -surface.debm_simple.refreeze_ice_melt no \
       -surface.debm_simple.air_temp_all_precip_as_snow 271.65 \
       -surface.debm_simple.air_temp_all_precip_as_rain 275.0\
       -time.start 2000-01-01 \
       -time.end 2001-01-01 \
       -time.calendar standard \
       -output.spatial.file debm_enhanced.nc \
       -output.spatial.vars mass_fluxes,pdd_fluxes,effective_air_temp,effective_precipitation,climatic_mass_balance,insolation \
       -output.spatial.times monthly
