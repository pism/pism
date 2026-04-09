#!/bin/bash
#
# Run SSA and Blatter inversions on the Wrangell glacier test case.
#
# Prerequisites:
#   1. Run prepare_wrangell.sh to download data and generate initial states
#   2. Observed velocities in obs_RGI2000-v7.0-C-01-04374.nc
#
# This script:
#   1. Generates synthetic SSA observations from the hybrid_fem state
#   2. Runs the SSA inversion using pismi.py
#   3. Runs the Blatter inversion using pismi_blatter.py

set -ex

NP=${NP:-8}
SCRIPTDIR=$(dirname "$0")

HYBRID_STATE=hybrid_fem_state_g500m_RGI2000-v7.0-C-01-04374_id_0_1978-01-01_1980-01-01.nc
BLATTER_STATE=blatter_state_g500m_RGI2000-v7.0-C-01-04374_id_0_1978-01-01_1980-01-01.nc

OBS=obs_RGI2000-v7.0-C-01-04374_0.nc

COMMON_PHYSICS="\
  -basal_resistance.pseudo_plastic.enabled yes \
  -basal_resistance.pseudo_plastic.q 0.75 \
  -basal_resistance.pseudo_plastic.u_threshold 100m/yr \
  -basal_yield_stress.model mohr_coulomb \
  -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 \
  -basal_yield_stress.mohr_coulomb.till_phi_default 30 \
"

# ── Step 1: Generate synthetic observations from the hybrid_fem state ──

echo "=== Generating synthetic observations ==="
mpirun -np ${NP} python ${SCRIPTDIR}/make_synth_ssa.py \
  -i ${HYBRID_STATE} \
  -o synth_ssa.nc \
  -generate_observed \
  -inv_ssa tauc \
  -design_prior_scale 0.9 \
  ${COMMON_PHYSICS} \
  -stress_balance.ssa.method fem

# ── Step 2: SSA inversion ──────────────────────────────────────────────

echo ""
echo "=== Running SSA inversion ==="
mpirun -np ${NP} python ${SCRIPTDIR}/pismi.py \
  -i ${HYBRID_STATE} \
  -inv_data ${OBS} \
  -o inv_ssa_result.nc \
  -inv_ssa tauc \
  -inv_method tikhonov_lmvm \
  -inv_target_misfit 50 \
  ${COMMON_PHYSICS} \
  -stress_balance.ssa.method fem

# ── Step 3: Generate synthetic observations from the blatter state ──

echo "=== Generating synthetic observations ==="
mpirun -np ${NP} python ${SCRIPTDIR}/make_synth_blatter.py \
  -i ${BLATTER_STATE} \
  -o synth_blatter.nc \
  -generate_observed \
  -inv_ssa tauc \
  -design_prior_scale 0.9 \
  ${COMMON_PHYSICS} \
    -basal_resistance.pseudo_plastic.enabled yes \
    -basal_resistance.pseudo_plastic.q 0.75 \
    -basal_resistance.pseudo_plastic.u_threshold 100m/yr \
    -stress_balance.blatter.Mz 17 \
    -stress_balance.blatter.coarsening_factor 4 \
    -bp_pc_type mg -bp_pc_mg_levels 3 \
    -bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type lu \
    -bp_mg_levels_ksp_type richardson -bp_mg_levels_pc_type sor \
    -bp_snes_rtol 0.001 -bp_ksp_rtol 0.001

# ── Step 3: Blatter inversion ─────────────────────────────────────────

echo ""
echo "=== Running Blatter inversion ==="
mpirun -np ${NP} python ${SCRIPTDIR}/pismi_blatter.py \
  -i ${BLATTER_STATE} \
  -inv_data ${OBS} \
  -o inv_blatter_result.nc \
  -inv_ssa tauc \
  -inv_method tikhonov_lmvm \
  -inv_target_misfit 50 \
  ${COMMON_PHYSICS} \
  -stress_balance.blatter.Mz 10 \
  -stress_balance.blatter.coarsening_factor 3 \
  -stress_balance.blatter.enhancement_factor 2.0 \
  -stress_balance.blatter.flow_law gpbld \
  -stress_balance.blatter.use_eta_transform yes \
  -bp_ksp_type gmres \
  -bp_ksp_rtol 0.01 \
  -bp_pc_type mg \
  -bp_pc_mg_levels 3 \
  -bp_mg_coarse_ksp_type preonly \
  -bp_mg_coarse_pc_type lu \
  -bp_mg_levels_ksp_type chebyshev \
  -bp_mg_levels_pc_type sor \
  -bp_mg_levels_ksp_max_it 3 \
  -bp_snes_rtol 0.01 \
  -bp_snes_max_it 20 \
  -bp_snes_ksp_ew 1 \
  -bp_snes_ksp_ew_version 3 \
  -bp_snes_linesearch_type bt

echo ""
echo "=== Done ==="
echo "Results:"
echo "  SSA:     inv_ssa_result.nc"
echo "  Blatter: inv_blatter_result.nc"
