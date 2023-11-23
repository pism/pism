#!/bin/bash

set -e
set -u
set -x

output_file=$1
input_file=$2
climate_file=$2

extra_vars=climatic_mass_balance,surface_accumulation_flux,surface_runoff_flux,air_temp_snapshot,snow_depth,surface_albedo,insolation,atmosphere_transmissivity,debms_insolation_driven_melt_flux,debms_temperature_driven_melt_flux,debms_background_melt_flux,thk,mask

common_options="
  -atmosphere searise_greenland \
  -atmosphere.searise_greenland.file ${climate_file} \
  -bootstrap \
  -energy none \
  -extra_file ${output_file} \
  -extra_times 7days \
  -extra_vars ${extra_vars} \
  -geometry.update.enabled false \
  -i ${input_file} \
  -o_size none \
  -stress_balance none \
  -ys 1-1-1 \
  -ye 2-12-30 \
"

# don't complain about unset variables
set +u

if [[ -n ${present_day} ]];
then
  pismr \
    ${common_options} \
    -surface debm_simple \
    ;
fi

if [[ -n ${paleo} ]];
then
  orbital_parameters=${3:-}
  pismr \
    ${common_options} \
    -surface debm_simple \
    -surface.debm_simple.paleo.enabled \
    -surface.debm_simple.paleo.file ${orbital_parameters} \
    ;
fi

if [[ -n ${albedo} ]];
then
  albedo_file=${3:-}
  pismr \
    ${common_options} \
    -surface debm_simple \
    -surface.debm_simple.albedo_input.file ${albedo_file} \
    ;
fi
