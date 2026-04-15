#!/bin/bash
#
# Parameter sweep for SSA inversion on the Wrangell glacier test case.
#
# Uses pre-generated observed velocities from make_synth_ssa.sh.
# Loops over penalty weight, H1/L2 regularization, and length scale.

set -ex

NP=${NP:-8}
SCRIPTDIR=$(dirname "$0")

STATE=state_blatter_g500m_RGI2000-v7.0-C-01-04374_id_0_1978-01-01_2023-01-01_0.nc
OBS=obs_RGI2000-v7.0-C-01-04374_0.nc

COMMON_PHYSICS="\
  -bp_ksp_rtol 0.001  \
  -bp_ksp_view_singularvalues   \
  -bp_mg_coarse_ksp_type preonly  \
  -bp_mg_coarse_pc_type lu  \
  -bp_mg_levels_ksp_type chebyshev \
  -bp_mg_levels_pc_type jacobi
  -bp_pc_mg_levels 3  \
  -bp_pc_type mg  \
  -bp_snes_ksp_ew 1  \
  -bp_snes_ksp_ew_version 3  \
  -bp_snes_monitor_ratio   \
  -bp_snes_rtol 0.001  \
  -stress_balance.blatter.Mz 10  \
  -stress_balance.blatter.coarsening_factor 3  \
  -stress_balance.blatter.enhancement_factor 2.0  \
  -stress_balance.blatter.flow_law gpbld  \
  -stress_balance.blatter.use_eta_transform yes  \
  -stress_balance.calving_front_stress_bc yes  \
  -stress_balance.model blatter  \
  -time_stepping.adaptive_ratio 500  \
  -basal_resistance.pseudo_plastic.enabled yes \
  -basal_resistance.pseudo_plastic.q 0.75 \
  -basal_resistance.pseudo_plastic.u_threshold 100m/yr \
  -basal_yield_stress.model mohr_coulomb \
  -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 \
  -basal_yield_stress.mohr_coulomb.till_phi_default 30 \
"

max_iter=100

for penalty in 0.1; do
    for h1 in 0.1; do
        for l2 in 0.1 ; do
            for scale in 1e3; do

                outfile=inv_blatter_it_${max_iter}_p_${penalty}_h1_${h1}_l2_${l2}_ls_${scale}.nc

                # Skip if output already exists
                if [ -f "${outfile}" ]; then
                    echo "=== Skipping ${outfile} (exists) ==="
                    continue
                fi

                echo ""
                echo "=== Blatter inv: penalty=${penalty} h1=${h1} l2=${l2} scale=${scale} ==="
                mpirun -np ${NP} python ${SCRIPTDIR}/pismi_blatter.py \
                    -i ${STATE} \
                    -inv_data ${OBS} \
                    -o ${outfile} \
                    -inv_blatter tauc \
                    -inv_method tikhonov_lmvm \
                    -inverse.stress_balance.length_scale ${scale} \
                    -inverse.design.cH1 ${h1} \
                    -inverse.design.cL2 ${l2} \
                    -inverse.max_iterations ${max_iter} \
                    -inverse.tikhonov.penalty_weight ${penalty} \
                    -inverse.use_zeta_fixed_mask yes \
                    ${COMMON_PHYSICS} \

            done
        done
    done
done

exit
    
for penalty in 0.1 1 10; do
    for h1 in 0 0.1 1 10; do
        for l2 in 0 0.1 1 10; do
            for scale in 1e2 5e2 1e3 5e3; do

                outfile=inv_blatter_it_${max_iter}_p_${penalty}_h1_${h1}_l2_${l2}_ls_${scale}.nc

                # Skip if output already exists
                if [ -f "${outfile}" ]; then
                    echo "=== Skipping ${outfile} (exists) ==="
                    continue
                fi

                echo ""
                echo "=== Blatter inv: penalty=${penalty} h1=${h1} l2=${l2} scale=${scale} ==="
                mpirun -np ${NP} python ${SCRIPTDIR}/pismi_blatter.py \
                    -i ${STATE} \
                    -inv_data ${OBS} \
                    -o ${outfile} \
                    -inv_blatter tauc \
                    -inv_method tikhonov_lmvm \
                    -inverse.stress_balance.length_scale ${scale} \
                    -inverse.design.cH1 ${h1} \
                    -inverse.design.cL2 ${l2} \
                    -inverse.max_iterations ${max_iter} \
                    -inverse.tikhonov.penalty_weight ${penalty} \
                    -inverse.use_zeta_fixed_mask yes \
                    ${COMMON_PHYSICS} \

            done
        done
    done
done

echo ""
echo "=== Sweep complete ==="
echo "Results:"
ls -1 inv_ssa_it${max_iter}_*.nc 2>/dev/null | wc -l
echo "files generated"
