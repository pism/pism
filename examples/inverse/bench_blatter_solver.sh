#!/bin/bash
#
# Benchmark Blatter solver options for inversion.
# Runs 1 TAO iteration with old vs new solver settings and compares wall time.

set -e

NP=${NP:-8}
SCRIPTDIR=$(dirname "$0")

STATE=blatter_state_g500m_RGI2000-v7.0-C-01-04374_id_0_1978-01-01_1980-01-01.nc
OBS=obs_RGI2000-v7.0-C-01-04374_0.nc

COMMON="\
  -i ${STATE} \
  -inv_data ${OBS} \
  -inv_ssa tauc \
  -inv_method tikhonov_lmvm \
  -inv_target_misfit 50 \
  -inverse.max_iterations 1 \
  -basal_resistance.pseudo_plastic.enabled yes \
  -basal_resistance.pseudo_plastic.q 0.75 \
  -basal_resistance.pseudo_plastic.u_threshold 100m/yr \
  -basal_yield_stress.model mohr_coulomb \
  -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 \
  -basal_yield_stress.mohr_coulomb.till_phi_default 30 \
  -stress_balance.blatter.Mz 17 \
  -stress_balance.blatter.coarsening_factor 4 \
  -stress_balance.blatter.enhancement_factor 2.0 \
  -stress_balance.blatter.flow_law gpbld \
  -stress_balance.blatter.use_eta_transform yes \
  -bp_pc_type mg \
  -bp_pc_mg_levels 3 \
  -bp_mg_coarse_ksp_type preonly \
  -bp_mg_coarse_pc_type lu \
  -bp_snes_ksp_ew 1 \
  -bp_snes_ksp_ew_version 3 \
"

OLD_OPTS="\
  -bp_ksp_rtol 0.001 \
  -bp_mg_levels_ksp_type richardson \
  -bp_mg_levels_pc_type sor \
  -bp_snes_rtol 0.001 \
"

NEW_OPTS="\
  -bp_ksp_type gmres \
  -bp_ksp_rtol 0.01 \
  -bp_mg_levels_ksp_type chebyshev \
  -bp_mg_levels_pc_type sor \
  -bp_mg_levels_ksp_max_it 3 \
  -bp_snes_rtol 0.01 \
  -bp_snes_max_it 20 \
  -bp_snes_linesearch_type bt \
"

echo "================================================================"
echo "Blatter inversion solver benchmark (1 TAO iteration, ${NP} procs)"
echo "================================================================"
echo ""

echo "--- OLD solver options (rtol=0.001, richardson+sor) ---"
t0=$(date +%s)
mpirun -np ${NP} python ${SCRIPTDIR}/pismi_blatter.py \
  ${COMMON} ${OLD_OPTS} \
  -o /tmp/bench_old.nc
t1=$(date +%s)
OLD_TIME=$((t1 - t0))
echo "OLD wall time: ${OLD_TIME}s"
echo ""

echo "--- NEW solver options (rtol=0.01, chebyshev+sor, gmres) ---"
t0=$(date +%s)
mpirun -np ${NP} python ${SCRIPTDIR}/pismi_blatter.py \
  ${COMMON} ${NEW_OPTS} \
  -o /tmp/bench_new.nc
t1=$(date +%s)
NEW_TIME=$((t1 - t0))
echo "NEW wall time: ${NEW_TIME}s"
echo ""

echo "================================================================"
echo "Summary:"
echo "  OLD: ${OLD_TIME}s"
echo "  NEW: ${NEW_TIME}s"
if [ ${OLD_TIME} -gt 0 ]; then
  SPEEDUP=$(echo "scale=1; ${OLD_TIME} / ${NEW_TIME}" | bc 2>/dev/null || echo "N/A")
  echo "  Speedup: ${SPEEDUP}x"
fi
echo "================================================================"
