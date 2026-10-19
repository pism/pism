#!/bin/bash

set -e
set -u
set -x

# Download example files

for f in bootfile_RGI2000-v7.0-C-01-04374.nc era5_wgs84_RGI2000-v7.0-C-01-04374.nc; do
    wget -nc https://pism-cloud-data.s3.amazonaws.com/debm_enhanced/$f;
done

res="200m"
interval="weekly"

common_options="
  -input.file bootfile_RGI2000-v7.0-C-01-04374.nc
  -input.bootstrap yes
  -grid.dx $res
  -grid.dy $res
  -energy.model none
  -stress_balance.model none
  -geometry.update.enabled no
  -atmosphere.models given,elevation_change
  -atmosphere.elevation_change.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc
  -atmosphere.elevation_change.temperature_lapse_rate 6
  -atmosphere.given.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc
  -atmosphere.given.periodic yes
  -time.start 2000-01-01
  -time.end 2001-01-01
  -time.calendar standard
  -output.spatial.vars mass_fluxes,pdd_fluxes,effective_air_temp,effective_precipitation,climatic_mass_balance,insolation,sky_view_factor
  -output.spatial.times $interval
  -output.spatial.stop_missing no
  -output.size none
"

debm_simple_options="
  -surface.debm_simple.interpret_precip_as_snow no
  -surface.debm_simple.c1 25
  -surface.debm_simple.c2 -103
  -surface.debm_simple.refreeze  0.50
  -surface.debm_simple.refreeze_ice_melt no
  -surface.debm_simple.air_temp_all_precip_as_snow 271.65
  -surface.debm_simple.air_temp_all_precip_as_rain 275.0
"

# PDD baseline with default melt factors
mpirun -np 8 pism ${common_options} \
       -surface pdd \
       -output.spatial.file pdd.nc

cdo yearmean pdd.nc pdd_ym.nc 

# dEBM Simple
mpirun -np 8 pism ${common_options} \
       -surface debm_simple ${debm_simple_options} \
       -output.spatial.file debm_simple.nc

cdo yearmean debm_simple.nc debm_simple_ym.nc

# dEBM Enhanced 

for df in 0.0 0.2; do
  mpirun -np 8 pism ${common_options} \
         -surface debm_simple ${debm_simple_options} \
         -surface.debm_simple.insolation_model terrain_shaded \
         -terrain_shaded_insolation.diffuse_fraction ${df} \
         -output.spatial.file debm_enhanced_df_${df}.nc

  cdo yearmean debm_enhanced_df_${df}.nc debm_enhanced_df_${df}_ym.nc
done
