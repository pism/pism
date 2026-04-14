#!/bin/bash
#
# Parameter sweep for SSA inversion on the Wrangell glacier test case.
#
# Uses pre-generated observed velocities from make_synth_ssa.sh.
# Loops over penalty weight, H1/L2 regularization, and length scale.

set -ex

NP=${NP:-8}
SCRIPTDIR=$(dirname "$0")

STATE=hybrid_fem_state_g500m_RGI2000-v7.0-C-01-04374_id_0_1978-01-01_1980-01-01_0.nc
OBS=obs_hybrid_fd_g500m_RGI2000-v7.0-C-01-04374_id_0_1978-01-01_2023-01-01_0.nc

COMMON_PHYSICS="\
  -basal_resistance.pseudo_plastic.enabled yes \
  -basal_resistance.pseudo_plastic.q 0.75 \
  -basal_resistance.pseudo_plastic.u_threshold 100m/yr \
  -basal_yield_stress.model mohr_coulomb \
  -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 \
  -basal_yield_stress.mohr_coulomb.till_phi_default 30 \
"

max_iter=100

for penalty in 0.1 1 10; do
    for h1 in 0.1 1 10; do
        for l2 in 0.1 1 10; do
            for scale in 1e3 5e3 1e4; do

                outfile=inv_ssa_it${max_iter}_p${penalty}_h1${h1}_l2${l2}_ls${scale}.nc

                # Skip if output already exists
                if [ -f "${outfile}" ]; then
                    echo "=== Skipping ${outfile} (exists) ==="
                    continue
                fi

                echo ""
                echo "=== SSA inv: penalty=${penalty} h1=${h1} l2=${l2} scale=${scale} ==="
                mpirun -np ${NP} python ${SCRIPTDIR}/pismi_ssa.py \
                    -i ${STATE} \
                    -inv_data ${OBS} \
                    -o ${outfile} \
                    -inv_ssa tauc \
                    -inv_method tikhonov_lmvm \
                    -inverse.stress_balance.length_scale ${scale} \
                    -inverse.design.cH1 ${h1} \
                    -inverse.design.cL2 ${l2} \
                    -inverse.max_iterations ${max_iter} \
                    -inverse.tikhonov.penalty_weight ${penalty} \
                    -inverse.use_zeta_fixed_mask yes \
                    -remove_sia \
                    ${COMMON_PHYSICS} \
                    -stress_balance.ssa.method fem

            done
        done
    done
done

echo ""
echo "=== Sweep complete ==="
echo "Results:"
ls -1 inv_ssa_it${max_iter}_*.nc 2>/dev/null | wc -l
echo "files generated"
