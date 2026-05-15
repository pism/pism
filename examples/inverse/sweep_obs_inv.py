#!/usr/bin/env python3
"""Generate inversion run scripts for real (ITS_LIVE) observations.

Usage:
    python sweep_obs_inv.py debug                      # mpirun scripts
    python sweep_obs_inv.py sbatch                     # SLURM scripts
    python sweep_obs_inv.py debug --restart FILE.nc    # restart from FILE.nc
"""

import argparse
import os
import sys
import itertools

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("mode", nargs="?", default="debug",
                    choices=["debug", "sbatch"])
parser.add_argument("--restart", metavar="FILE",
                    help="Restart inversion from zeta_inv in FILE "
                         "(e.g., a previous SSA inversion result)")
args = parser.parse_args()

MODE = args.mode
RESTART_FILE = args.restart
NP = int(os.environ.get("NP", 40 if MODE == "sbatch" else 8))
SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))

state = "state_blatter_g500m_RGI2000-v7.0-C-01-04374_id_0_1980-01-01_1990-01-01_0.nc"
state = "state_hybrid_g500m_RGI2000-v7.0-C-01-04374_id_0_1980-01-01_2020-01-01_0.nc"
OBS = "obs_RGI2000-v7.0-C-01-04374.nc"

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
    "-hydrology.model", "null",
    "-hydrology.null_diffuse_till_water", "",
]

BLATTER_PHYSICS = [
    "-stress_balance.model", "blatter",
    "-bp_ksp_rtol", "1e-4",
    "-bp_mg_coarse_ksp_type", "preonly",
    "-bp_mg_coarse_pc_type", "lu",
    "-bp_mg_levels_ksp_type", "chebyshev",
    "-bp_pc_mg_levels", "3",
    "-bp_pc_type", "mg",
    "-bp_snes_ksp_ew", "0",
    "-bp_snes_ksp_ew_version", "3",
    "-bp_snes_rtol", "1e-4",
    "-stress_balance.blatter.Mz", "10",
    "-stress_balance.blatter.coarsening_factor", "3",
    "-stress_balance.blatter.flow_law", "isothermal_glen",
    "-stress_balance.blatter.use_eta_transform", "yes",
    "-stress_balance.calving_front_stress_bc", "yes",
    "-stress_balance.blatter.enhancement_factor", "2.0"
    "-stress_balance.blatter.flow_law", "gpbld"
    "-time_stepping.adaptive_ratio", "10",
    "-inv_adj_ksp_type", "gmres",
    "-inv_adj_pc_type", "jacobi",
]


HYBRID_PHYSICS = [
    "-stress_balance.model", "ssa+sia",
    "-stress_balance.ssa.method", "fem",
    "-stress_balance.sia.enhancement_factor", "2.0",
    "-stress_balance.sia.flow_law", "gpbld",
    "-stress_balance.sia.max_diffusivity", "100000.0",
    "-stress_balance.sia.surface_gradient_method", "eta",
    "-stress_balance.ssa.flow_law", "gpbld",
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
    "hybrid": {"inv_flag": ["-inv_design", "tauc"], "physics": HYBRID_PHYSICS},
    "blatter": {"inv_flag": ["-inv_design", "tauc"], "physics": BLATTER_PHYSICS},
}


penalties = [1, 10, 100, 1000, 10000, 100000]
h1_values = [1]
l2_values = [0, 0.1]
hscales = ["1e1", "1e2", "1e3"]
vscales = [100]

count = 0

for sb, params in solvers.items():
    for penalty, h1, l2, hscale, vscale in itertools.product(
            penalties, h1_values, l2_values, hscales, vscales):

        outfile = f"inv_obs_{sb}_it_{max_iter}_p_{penalty}_h1_{h1}_l2_{l2}_ls_{hscale}_vs_{vscale}.nc"
        jobname = f"inv_obs_{sb}_it_{max_iter}_p_{penalty}_h1_{h1}_l2_{l2}_ls_{hscale}_vs_{vscale}"
        runscript = os.path.join(scriptdir, f"{jobname}.sh")

        # If restarting, extract zeta_inv from the restart file into a
        # temporary inv_data file so pismi picks it up as the initial guess.
        restart_cmds = ""
        inv_data = OBS
        if RESTART_FILE is not None:
            inv_data = f"inv_data_restart_{sb}.nc"
            restart_cmds = (
                f"# Extract tauc (last time, squeezed to 2D) from restart file,\n"
                f"# rename to tauc_prior, and merge into a copy of the obs file.\n"
                f"# Using tauc (not zeta_inv) avoids NaN issues in ice-free regions.\n"
                f"ncks -O -d time,-1 -v tauc {RESTART_FILE} _tauc_tmp.nc\n"
                f"ncwa -O -a time _tauc_tmp.nc _tauc_tmp.nc\n"
                f"ncrename -v tauc,tauc_prior _tauc_tmp.nc\n"
                f"cp {OBS} {inv_data}\n"
                f"ncks -A -C -v tauc_prior _tauc_tmp.nc {inv_data}\n"
                f"rm -f _tauc_tmp.nc\n"
            )

        cmd_parts = [
            RUN_CMD, "python", os.path.join(SCRIPTDIR, pyscript),
            "-i", state,
            "-inv_data", inv_data,
            "-o", outfile,
            *params["inv_flag"],
            "-inverse.stress_balance.velocity_scale", str(vscale),
            "-inverse.stress_balance.length_scale", hscale,
            "-inverse.design.cH1", str(h1),
            "-inverse.design.cL2", str(l2),
            "-inverse.max_iterations", str(max_iter),
            "-inverse.design.param", "exp",
            "-inverse.stress_balance.method", "tikhonov_lmvm",
            "-inverse.tikhonov.penalty_weight", str(penalty),
            "-inverse.use_zeta_fixed_mask", "yes",
            "-inverse.adjoint.method", "exact",
            "-inv_grounded_ice_tauc",
            *COMMON_PHYSICS,
            *params["physics"],
        ]

        with open(runscript, "w") as f:
            f.write(HEADER)
            f.write(f"#SBATCH --job-name={jobname}\n\n")
            if restart_cmds:
                f.write(restart_cmds)
            f.write(" ".join(cmd_parts) + "\n")

        os.chmod(runscript, 0o755)
        count += 1

print(f"\nGenerated {count} scripts in {scriptdir}/\n")

if MODE == "sbatch":
    print(f"Submit all with:\n  for f in {scriptdir}/*.sh; do sbatch $f; done")
else:
    print(f"Run all with:\n  for f in {scriptdir}/*.sh; do bash $f; done")
    print(f"\nOr run one:\n  bash {scriptdir}/<script>.sh")
