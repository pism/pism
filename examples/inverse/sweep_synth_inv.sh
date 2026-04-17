#!/bin/bash
#
# Parameter sweep for SSA and Blatter inversions on the Wrangell glacier.
#
# Uses pre-generated observed velocities from make_synth.sh.
#
# Usage:
#   bash sweep_inv.sh debug    # generates scripts that run directly with mpirun
#   bash sweep_inv.sh sbatch   # generates scripts for SLURM submission
#   bash sweep_inv.sh          # defaults to debug

set -e

MODE=${1:-debug}
if [ "$MODE" = "sbatch" ]; then
    NP=${NP:-40}
else
    NP=${NP:-8}
fi
SCRIPTDIR=$(dirname "$0")

start="1978-01-01"
end="2023-01-01"
res=500

SBATCH_HEADER='#!/bin/sh
#SBATCH --partition=t2small
#SBATCH --ntasks=40
#SBATCH --tasks-per-node=40
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=pism.%j
'

DEBUG_HEADER='#!/bin/bash
set -ex
'

RUN_CMD="mpirun -np ${NP}"
if [ "$MODE" = "sbatch" ]; then
    HEADER="$SBATCH_HEADER"
else
    HEADER="$DEBUG_HEADER"
fi

COMMON_PHYSICS="\
  -basal_resistance.pseudo_plastic.enabled yes \
  -basal_resistance.pseudo_plastic.q 0.75 \
  -basal_resistance.pseudo_plastic.u_threshold 100m/yr \
  -basal_yield_stress.model mohr_coulomb \
  -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025 \
  -basal_yield_stress.mohr_coulomb.till_phi_default 30"

BLATTER_PHYSICS="\
  -bp_ksp_rtol 0.001 \
  -bp_mg_coarse_ksp_type preonly \
  -bp_mg_coarse_pc_type lu \
  -bp_mg_levels_ksp_type chebyshev \
  -bp_mg_levels_pc_type jacobi \
  -bp_pc_mg_levels 3 \
  -bp_pc_type mg \
  -bp_snes_ksp_ew 1 \
  -bp_snes_ksp_ew_version 3 \
  -bp_snes_rtol 0.001 \
  -stress_balance.blatter.Mz 10 \
  -stress_balance.blatter.coarsening_factor 3 \
  -stress_balance.blatter.enhancement_factor 2.0 \
  -stress_balance.blatter.flow_law gpbld \
  -stress_balance.blatter.use_eta_transform yes \
  -stress_balance.calving_front_stress_bc yes \
  -time_stepping.adaptive_ratio 10"

SSA_PHYSICS="\
  -stress_balance.ssa.method fem"

max_iter=1000
scriptdir="run_synth_scripts"
mkdir -p ${scriptdir}

count=0

pyscript="pismi.py"
for sb in ssa hybrid blatter; do
    if [ "$sb" = "hybrid" ]; then
        inv_flag="-inv_design tauc"
        sb_physics="${SSA_PHYSICS}"
    elif [ "$sb" = "ssa" ]; then
        inv_flag="-inv_design tauc"
        sb_physics="${SSA_PHYSICS}"
    else
        inv_flag="-inv_design tauc"
        sb_physics="${BLATTER_PHYSICS}"
    fi

    for penalty in 0.1 1 10; do
        for h1 in 0 1 10; do
            for l2 in 0 1 10; do
                for hscale in 1e3; do
                    for vscale in 100; do

                        STATE=state_${sb}_g${res}m_RGI2000-v7.0-C-01-04374_id_0_${start}_${end}_0.nc
                        OBS=synth_obs_${sb}_g${res}m_RGI2000-v7.0-C-01-04374_id_0_${start}_${end}.nc
                        outfile=inv_synth_${sb}_it_${max_iter}_p_${penalty}_h1_${h1}_l2_${l2}_ls_${hscale}_vs_${vscale}.nc
                        jobname=inv_synth_${sb}_p_${penalty}_h1_${h1}_l2_${l2}_ls_${hscale}_vs_${vscale}
                        runscript=${scriptdir}/${jobname}.sh

                        cat > ${runscript} <<EOF
${HEADER}
#SBATCH --job-name=${jobname}

${RUN_CMD} python ${SCRIPTDIR}/${pyscript} \\
  -i ${STATE} \\
  -inv_data ${OBS} \\
  -o ${outfile} \\
  ${inv_flag} \\
  -inv_method tikhonov_lmvm \\
  -inverse.stress_balance.velocity_scale ${vscale} \\
  -inverse.stress_balance.length_scale ${hscale} \\
  -inverse.design.cH1 ${h1} \\
  -inverse.design.cL2 ${l2} \\
  -inverse.max_iterations ${max_iter} \\
  -inverse.tikhonov.penalty_weight ${penalty} \\
  -inverse.use_zeta_fixed_mask yes \\
  ${COMMON_PHYSICS} \\
  ${sb_physics}
EOF
                        chmod +x ${runscript}
                        count=$((count + 1))

                    done
                done
            done
        done
    done
done

echo ""
echo "Generated ${count} scripts in ${scriptdir}/"
echo ""

if [ "$MODE" = "sbatch" ]; then
    echo "Submit all with:"
    echo "  for f in ${scriptdir}/*.sh; do sbatch \$f; done"
else
    echo "Run all with:"
    echo "  for f in ${scriptdir}/*.sh; do bash \$f; done"
    echo ""
    echo "Or run one:"
    echo "  bash ${scriptdir}/<script>.sh"
fi
