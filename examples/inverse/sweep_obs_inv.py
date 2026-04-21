#!/usr/bin/env python3
"""Generate inversion run scripts for real (ITS_LIVE) observations.

Usage:
    python sweep_obs_inv.py debug    # mpirun scripts
    python sweep_obs_inv.py sbatch   # SLURM scripts
"""

import os
import sys
import itertools

MODE = sys.argv[1] if len(sys.argv) > 1 else "debug"
NP = int(os.environ.get("NP", 40 if MODE == "sbatch" else 8))
SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))

start = "1978-01-01"
end = "2023-01-01"
res = 500

OBS = "obs_RGI2000-v7.0-C-01-04374_0.nc"

SBATCH_HEADER = """\
#!/bin/sh
#SBATCH --partition=t2small
#SBATCH --ntasks=40
#SBATCH --tasks-per-node=40
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=pism.%j
"""

DEBUG_HEADER = """\
#!/bin/bash
set -ex
"""

HEADER = SBATCH_HEADER if MODE == "sbatch" else DEBUG_HEADER
RUN_CMD = f"mpirun -np {NP}"

COMMON_PHYSICS = [
    "-basal_resistance.pseudo_plastic.enabled", "yes",
    "-basal_resistance.pseudo_plastic.q", "0.75",
    "-basal_resistance.pseudo_plastic.u_threshold", "100m/yr",
    "-basal_yield_stress.model", "mohr_coulomb",
    "-basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden", "0.025",
    "-basal_yield_stress.mohr_coulomb.till_phi_default", "30",
]

BLATTER_PHYSICS = [
    "-stress_balance.model", "blatter",
    "-bp_ksp_rtol", "0.001",
    "-bp_mg_coarse_ksp_type", "preonly",
    "-bp_mg_coarse_pc_type", "lu",
    "-bp_mg_levels_ksp_type", "chebyshev",
    "-bp_pc_mg_levels", "3",
    "-bp_pc_type", "mg",
    "-bp_snes_ksp_ew", "1",
    "-bp_snes_ksp_ew_version", "3",
    "-bp_snes_rtol", "0.001",
    "-stress_balance.blatter.Mz", "10",
    "-stress_balance.blatter.coarsening_factor", "3",
    "-stress_balance.blatter.enhancement_factor", "2.0",
    "-stress_balance.blatter.flow_law", "gpbld",
    "-stress_balance.blatter.use_eta_transform", "yes",
    "-stress_balance.calving_front_stress_bc", "yes",
    "-time_stepping.adaptive_ratio", "10",
    "-inverse.use_incomplete_adjoint", "no",
    "-inv_adj_pc_type", "jacobi",
]

SSA_PHYSICS = [
    "-stress_balance.model", "ssa",
    "-stress_balance.ssa.method", "fem",
]

max_iter = 250
scriptdir = "run_obs_scripts"
os.makedirs(scriptdir, exist_ok=True)

pyscript = "pismi.py"

solvers = {
    "ssa": {"inv_flag": ["-inv_design", "tauc"], "physics": SSA_PHYSICS},
    "blatter": {"inv_flag": ["-inv_design", "tauc"], "physics": BLATTER_PHYSICS},
}

penalties = [0.01, 1, 10, 100, 1000]
h1_values = [0.01, 1, 10]
l2_values = [0, 1, 10]
hscales = ["5e3", "5e4"]
vscales = [100]

count = 0

for sb, params in solvers.items():
    for penalty, h1, l2, hscale, vscale in itertools.product(
            penalties, h1_values, l2_values, hscales, vscales):

        state = f"state_{sb}_g{res}m_RGI2000-v7.0-C-01-04374_id_0_{start}_{end}_0.nc"
        outfile = f"inv_obs_{sb}_it_{max_iter}_p_{penalty}_h1_{h1}_l2_{l2}_ls_{hscale}_vs_{vscale}.nc"
        jobname = f"inv_obs_{sb}_p_{penalty}_h1_{h1}_l2_{l2}_ls_{hscale}_vs_{vscale}"
        runscript = os.path.join(scriptdir, f"{jobname}.sh")

        cmd_parts = [
            RUN_CMD, "python", os.path.join(SCRIPTDIR, pyscript),
            "-i", state,
            "-inv_data", OBS,
            "-o", outfile,
            *params["inv_flag"],
            "-inverse.stress_balance.velocity_scale", str(vscale),
            "-inverse.stress_balance.length_scale", hscale,
            "-inverse.design.cH1", str(h1),
            "-inverse.design.cL2", str(l2),
            "-inverse.max_iterations", str(max_iter),
            "-inverse.design.param", "exp",
            "-inverse.stress_balance.method", "tikhonov_lmvm",
            "-inverse.tikhonov.atol", "1e-10",
            "-inverse.tikhonov.penalty_weight", str(penalty),
            "-inverse.tikhonov.rtol", "1e-10",
            "-inverse.use_zeta_fixed_mask", "yes",
            "-inv_grounded_ice_tauc",
            "-tao_frtol", "1e-20",
            "-tao_fatol", "1e-20",
            *COMMON_PHYSICS,
            *params["physics"],
        ]

        with open(runscript, "w") as f:
            f.write(HEADER)
            f.write(f"#SBATCH --job-name={jobname}\n\n")
            f.write(" ".join(cmd_parts) + "\n")

        os.chmod(runscript, 0o755)
        count += 1

print(f"\nGenerated {count} scripts in {scriptdir}/\n")

if MODE == "sbatch":
    print(f"Submit all with:\n  for f in {scriptdir}/*.sh; do sbatch $f; done")
else:
    print(f"Run all with:\n  for f in {scriptdir}/*.sh; do bash $f; done")
    print(f"\nOr run one:\n  bash {scriptdir}/<script>.sh")
