#!/bin/bash
#
# Benchmark Blatter solver options for inversion.
# Runs 1 TAO iteration with the OLD baseline vs the validated NEW solver recipe
# and compares wall time.
#
#   OLD: direct LU on the coarse level (-bp_mg_coarse_pc_type lu, preonly),
#        plain GMRES outer solve.
#   NEW: flexible outer Krylov + Krylov-wrapped GAMG coarse solve ("gamg_ksp").
#        -bp_ksp_type fgmres is REQUIRED: the GMRES coarse solve makes the MG
#        preconditioner a variable operator, which plain GMRES cannot handle
#        (it diverges with DIVERGED_LINEAR_SOLVE). Everything else (mg levels,
#        smoother, tolerances) is held fixed in COMMON so the timing difference
#        is attributable to the *forward* Blatter coarse solve.
#
# The inverse setup (inverse.*) and the adjoint linear solver (inv_adj_*) are held
# fixed in COMMON, matching a known-good inversion configuration. The adjoint needs
# its own preconditioner (inv_adj_ksp_type gmres / inv_adj_pc_type jacobi) plus the
# functional scaling (velocity_scale / length_scale / huber); without them the
# adjoint KSP diverges (DIVERGED_ITS) and TaoSolve fails. If the adjoint still
# struggles on a harder problem, try inv_adj_pc_type gamg.
#
# Override the input files (the originals are gone) and rank count via env, e.g.
#   NP=8 STATE=my_state.nc OBS=my_obs.nc ./bench_blatter_inv.sh

set -e

NP=${NP:-8}
SCRIPTDIR=$(dirname "$0")

STATE=${STATE:-blatter_state_g500m_RGI2000-v7.0-C-01-04374_id_0_1978-01-01_1980-01-01.nc}
OBS=${OBS:-obs_RGI2000-v7.0-C-01-04374_0.nc}

COMMON="\
  -i ${STATE} \
  -inv_data ${OBS} \
  -inverse.max_iterations 1 \
  -inverse.stress_balance.method tikhonov_blmvm \
  -inverse.design.param exp \
  -inverse.design.cH1 1 \
  -inverse.design.cL2 0 \
  -inverse.tikhonov.penalty_weight 10 \
  -inverse.stress_balance.tauc_min 1e4 \
  -inverse.stress_balance.tauc_max 5e7 \
  -inverse.stress_balance.velocity_scale 1e3 \
  -inverse.stress_balance.length_scale 1e6 \
  -inverse.state_func huber \
  -inverse.huber.delta 1e3 \
  -inverse.use_zeta_fixed_mask yes \
  -inverse.adjoint.method approximate \
  -inv_adj_ksp_type gmres \
  -inv_adj_pc_type jacobi \
  -basal_resistance.pseudo_plastic.enabled yes \
  -basal_resistance.pseudo_plastic.q 0.75 \
  -basal_resistance.pseudo_plastic.u_threshold 100m/yr \
  -basal_yield_stress.model mohr_coulomb \
  -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 \
  -basal_yield_stress.mohr_coulomb.till_phi_default 30 \
  -stress_balance.model blatter \
  -stress_balance.blatter.Mz 10 \
  -stress_balance.blatter.coarsening_factor 3 \
  -stress_balance.blatter.enhancement_factor 2.0 \
  -stress_balance.blatter.flow_law gpbld \
  -stress_balance.blatter.use_eta_transform yes \
  -bp_pc_type mg \
  -bp_pc_mg_levels 3 \
  -bp_mg_levels_ksp_type richardson \
  -bp_mg_levels_pc_type sor \
  -bp_ksp_rtol 1e-3 \
  -bp_snes_rtol 1e-3 \
  -bp_snes_ksp_ew 0 \
  -bp_snes_converged_reason \
"

# OLD: direct LU on the coarse level (original default), plain GMRES outer solve.
OLD_OPTS="\
  -bp_mg_coarse_ksp_type preonly \
  -bp_mg_coarse_pc_type lu \
"

# NEW: flexible outer Krylov + Krylov-wrapped GAMG coarse solve.
NEW_OPTS="\
  -bp_ksp_type fgmres \
  -bp_mg_coarse_pc_type gamg \
  -bp_mg_coarse_ksp_type gmres \
  -bp_mg_coarse_ksp_rtol 1e-2 \
  -bp_mg_coarse_ksp_max_it 50 \
"

echo "================================================================"
echo "Blatter inversion solver benchmark (1 TAO iteration, ${NP} procs)"
echo "================================================================"
echo ""

echo "--- OLD solver options (direct LU coarse, plain GMRES) ---"
t0=$(date +%s)
mpirun -np ${NP} pismi \
  ${COMMON} ${OLD_OPTS} \
  -o bench_old.nc
t1=$(date +%s)
OLD_TIME=$((t1 - t0))
echo "OLD wall time: ${OLD_TIME}s"
echo ""

echo "--- NEW solver options (fgmres + gamg_ksp coarse solve) ---"
t0=$(date +%s)
mpirun -np ${NP} pismi \
  ${COMMON} ${NEW_OPTS} \
  -o bench_new.nc
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
