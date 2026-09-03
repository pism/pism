#!/bin/bash
# Idealized debris-covered valley glacier (Verhaegen and Huybrechts, 2026): a
# higher-order (Blatter) isothermal glacier with pseudo-plastic sliding, with the debris
# transport model running as a passive tracer (the effect of debris on melt is not yet
# fed back to the surface mass balance).
#
# The Blatter solver does not cope with the thin, freshly nucleated ice of the first
# years at grid spacings below ~100 m (the linear solve diverges), so the glacier is first
# grown with the hybrid SSA+SIA stress balance (same sliding law and basal strength) for
# SPINUP years and the Blatter run restarts from that state.
#
# Usage: ./run.sh [input file] [run length in years] [output prefix]
# Environment: N (number of MPI processes, default 4), DX (grid spacing in m, default 25),
#              SPINUP (length of the SSA+SIA spin-up in years, default 30; 0 to skip),
#              DEBRIS_INPUT (debris forcing file, default debris_input.nc)

set -euo pipefail

N=${N:-4}
DX=${DX:-25}
SPINUP=${SPINUP:-30}

input=${1:-input.nc}
duration=${2:-300}
prefix=${3:-debris}
debris_input=${DEBRIS_INPUT:-debris_input.nc}

if [ ! -f "${input}" ] || [ ! -f "${debris_input}" ]; then
  python3 create_input.py --dx ${DX} --debris-file "${debris_input}" "${input}"
fi

spatial_vars=thk,usurf,mask,velsurf_mag,uvel,vvel,wvel,climatic_mass_balance
spatial_vars=${spatial_vars},debris_thickness,englacial_debris_concentration,englacial_debris_column_mass
spatial_vars=${spatial_vars},debris_cover_fraction,debris_melt_out_rate,debris_input_rate,debris_removal_rate
spatial_vars=${spatial_vars},debris_surface_velocity,debris_gravitational_flux,ice_melt_enhancement

scalar_vars=ice_volume,ice_area_glacierized,englacial_debris_mass,supraglacial_debris_mass,total_debris_mass
scalar_vars=${scalar_vars},debris_input_mass_flux,debris_melt_out_mass_flux,debris_output_mass_flux
scalar_vars=${scalar_vars},debris_lost_mass_flux,debris_mass_conservation_error

# options shared by the spin-up and the main run
common="
  -grid.Mz 31 -grid.Lz 300
  -stress_balance.sia.surface_gradient_method eta
  -stress_balance.sia.max_diffusivity 100000.0
  -stress_balance.sia.bed_smoother.range 0
  -stress_balance.sia.flow_law isothermal_glen
  -flow_law.isothermal_Glen.ice_softness 1.0e-24
  -basal_resistance.pseudo_plastic.enabled yes
  -basal_resistance.pseudo_plastic.q 0.75
  -basal_resistance.pseudo_plastic.u_threshold 100.m.year-1
  -basal_yield_stress.model mohr_coulomb
  -basal_yield_stress.mohr_coulomb.till_phi_default 35
  -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025
  -time_stepping.adaptive_ratio 50
  -time_stepping.maximum_time_step 1
  -energy none
  -surface given
  -geometry.front_retreat.prescribed.file ${input}
  -debris transport
  -debris.transport.input.file ${debris_input}
  -debris.transport.marginal_length_scale ${DX}
  -debris.ice_melt_enhancement.model verhaegen
  -output.sizes.medium ${spatial_vars}
"

# higher-order ice dynamics (Blatter) with pseudo-plastic sliding
blatter="
  -stress_balance.model blatter
  -stress_balance.blatter.Mz 10
  -stress_balance.blatter.coarsening_factor 3
  -stress_balance.blatter.use_eta_transform yes
  -stress_balance.blatter.flow_law isothermal_glen
  -stress_balance.calving_front_stress_bc no
"

if [ "${SPINUP}" -gt 0 ]; then
  # grow the glacier with the hybrid SSA+SIA model using the same sliding law; the
  # Blatter run then restarts from it, warm-started with the SSA velocities
  mpiexec -n ${N} pism \
    -i "${input}" -bootstrap \
    ${common} \
    -stress_balance.model ssa+sia \
    -output.file "${prefix}_spinup.nc" \
    -time.run_length ${SPINUP} \
    2>&1 | tee "${prefix}_spinup.log"

  start="-i ${prefix}_spinup.nc"
else
  start="-i ${input} -bootstrap"
fi

mpiexec -n ${N} pism \
  ${start} \
  ${common} \
  ${blatter} \
  -output.file "${prefix}.nc" \
  -output.spatial.file "${prefix}_spatial.nc" \
  -output.spatial.times 10 \
  -output.spatial.vars "${spatial_vars}" \
  -output.scalar.file "${prefix}_scalar.nc" \
  -output.scalar.times 1 \
  -output.scalar.variables "${scalar_vars}" \
  -time.run_length ${duration} \
  2>&1 | tee "${prefix}.log"
