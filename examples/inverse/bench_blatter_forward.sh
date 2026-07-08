#!/bin/bash
#
# Benchmark Blatter solver options for the forward PISM run.
# Runs a short 1-month simulation with different solver configs.

set -e

NP=${NP:-6}

STATE=bootfile_g50m_RGI2000-v7.0-C-01-04374.nc
GRID=grid_g50m_RGI2000-v7.0-C-01-04374.nc
CLIMATE=era5_wgs84_RGI2000-v7.0-C-01-04374.nc

COMMON="\
  -atmosphere.given.file ${CLIMATE} \
  -atmosphere.models given \
  -basal_resistance.pseudo_plastic.enabled yes \
  -basal_resistance.pseudo_plastic.q 0.75 \
  -basal_resistance.pseudo_plastic.u_threshold 100m/yr \
  -basal_yield_stress.model mohr_coulomb \
  -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 \
  -basal_yield_stress.mohr_coulomb.till_phi_default 30 \
  -energy.model enthalpy \
  -geometry.front_retreat.use_cfl yes \
  -geometry.part_grid.enabled yes \
  -geometry.remove_icebergs yes \
  -grid.Lbz 0 \
  -grid.Lz 2000 \
  -grid.Mbz 1 \
  -grid.Mz 101 \
  -grid.dx 500m \
  -grid.dy 500m \
  -grid.file ${GRID} \
  -grid.registration center \
  -input.bootstrap yes \
  -input.file ${STATE} \
  -input.forcing.buffer_size 390 \
  -input.forcing.time_extrapolation yes \
  -output.size none \
  -surface.debm_simple.c1 30 \
  -surface.debm_simple.c2 -120 \
  -surface.debm_simple.interpret_precip_as_snow no \
  -surface.force_to_thickness.file ${STATE} \
  -surface.models debm_simple \
  -time.calendar standard \
  -time.end 1978-01-10 \
  -time.reference_date 1978-01-01 \
  -time.start 1978-01-01 \
  -time_stepping.skip.enabled yes \
  -time_stepping.skip.max 100 \
  -stress_balance.model blatter \
  -stress_balance.calving_front_stress_bc yes \
  -stress_balance.blatter.enhancement_factor 2.0 \
  -stress_balance.blatter.flow_law gpbld \
  -stress_balance.blatter.use_eta_transform yes \
  -time_stepping.adaptive_ratio 200 \
  -bp_snes_ksp_ew 1 \
  -bp_snes_ksp_ew_version 3 \
"

declare -A CONFIGS

# A: Current settings (baseline)
CONFIGS[A_baseline]="\
  -stress_balance.blatter.Mz 10 \
  -stress_balance.blatter.coarsening_factor 3 \
  -bp_pc_type mg \
  -bp_pc_mg_levels 3 \
  -bp_ksp_rtol 0.001 \
  -bp_mg_coarse_ksp_type preonly \
  -bp_mg_coarse_pc_type lu \
  -bp_mg_levels_ksp_type richardson \
  -bp_mg_levels_pc_type sor \
  -bp_snes_rtol 0.001 \
"

# B: Chebyshev smoother instead of Richardson
CONFIGS[B_chebyshev]="\
  -stress_balance.blatter.Mz 10 \
  -stress_balance.blatter.coarsening_factor 3 \
  -bp_pc_type mg \
  -bp_pc_mg_levels 3 \
  -bp_ksp_rtol 0.001 \
  -bp_mg_coarse_ksp_type preonly \
  -bp_mg_coarse_pc_type lu \
  -bp_mg_levels_ksp_type chebyshev \
  -bp_mg_levels_pc_type sor \
  -bp_mg_levels_ksp_max_it 3 \
  -bp_snes_rtol 0.001 \
"

# C: Chebyshev + looser tolerances
CONFIGS[C_chebyshev_loose]="\
  -stress_balance.blatter.Mz 10 \
  -stress_balance.blatter.coarsening_factor 3 \
  -bp_pc_type mg \
  -bp_pc_mg_levels 3 \
  -bp_ksp_rtol 0.01 \
  -bp_mg_coarse_ksp_type preonly \
  -bp_mg_coarse_pc_type lu \
  -bp_mg_levels_ksp_type chebyshev \
  -bp_mg_levels_pc_type sor \
  -bp_mg_levels_ksp_max_it 3 \
  -bp_snes_rtol 0.01 \
"

# D: GMRES + Chebyshev MG + loose
CONFIGS[D_gmres_chebyshev]="\
  -stress_balance.blatter.Mz 10 \
  -stress_balance.blatter.coarsening_factor 3 \
  -bp_ksp_type gmres \
  -bp_pc_type mg \
  -bp_pc_mg_levels 3 \
  -bp_ksp_rtol 0.01 \
  -bp_mg_coarse_ksp_type preonly \
  -bp_mg_coarse_pc_type lu \
  -bp_mg_levels_ksp_type chebyshev \
  -bp_mg_levels_pc_type sor \
  -bp_mg_levels_ksp_max_it 3 \
  -bp_snes_rtol 0.01 \
  -bp_snes_linesearch_type bt \
"

# E: Hypre BoomerAMG (algebraic MG, no geometric MG)
CONFIGS[E_hypre]="\
  -stress_balance.blatter.Mz 10 \
  -stress_balance.blatter.coarsening_factor 3 \
  -bp_ksp_type gmres \
  -bp_pc_type hypre \
  -bp_pc_hypre_type boomeramg \
  -bp_ksp_rtol 0.01 \
  -bp_snes_rtol 0.01 \
"

# F: Fewer vertical levels (Mz=7, 2 MG levels)
CONFIGS[F_fewer_levels]="\
  -stress_balance.blatter.Mz 7 \
  -stress_balance.blatter.coarsening_factor 3 \
  -bp_pc_type mg \
  -bp_pc_mg_levels 2 \
  -bp_ksp_rtol 0.01 \
  -bp_mg_coarse_ksp_type preonly \
  -bp_mg_coarse_pc_type lu \
  -bp_mg_levels_ksp_type chebyshev \
  -bp_mg_levels_pc_type sor \
  -bp_mg_levels_ksp_max_it 3 \
  -bp_snes_rtol 0.01 \
"

echo "================================================================"
echo "Blatter forward solver benchmark (10-day run, ${NP} procs)"
echo "================================================================"
echo ""

declare -A TIMES

for config in $(echo "${!CONFIGS[@]}" | tr ' ' '\n' | sort); do
    echo "--- ${config} ---"
    t0=$(date +%s)
    mpirun -np ${NP} pism \
      ${COMMON} ${CONFIGS[$config]} \
      -output.file /tmp/bench_${config}.nc \
      : -n 1 pism_async_writer \
      2>&1 || true
    t1=$(date +%s)
    TIMES[$config]=$((t1 - t0))
    echo "${config} wall time: ${TIMES[$config]}s"
    echo ""
done

echo "================================================================"
echo "Summary (10-day Blatter forward run, ${NP} procs):"
echo "================================================================"
for config in $(echo "${!CONFIGS[@]}" | tr ' ' '\n' | sort); do
    printf "  %-25s %4ds\n" "${config}" "${TIMES[$config]}"
done
baseline=${TIMES[A_baseline]}
if [ "${baseline}" -gt 0 ] 2>/dev/null; then
    echo ""
    echo "Relative to baseline:"
    for config in $(echo "${!CONFIGS[@]}" | tr ' ' '\n' | sort); do
        speedup=$(echo "scale=2; ${baseline} / ${TIMES[$config]}" | bc 2>/dev/null || echo "N/A")
        printf "  %-25s %sx\n" "${config}" "${speedup}"
    done
fi
echo "================================================================"
