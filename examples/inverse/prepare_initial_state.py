#!/usr/bin/env python3
"""Generate PISM spin-up run scripts for the initial-state preparation step.

Produces one shell script per (stress_balance, hydrology) combination, under
``run_prep_scripts/``. Mirrors the structure of ``sweep_obs_inv.py``.

Usage::

    python prepare_initial_state.py debug                # mpirun scripts
    python prepare_initial_state.py sbatch               # SLURM scripts

The emitted scripts:

1. Fetch the input NetCDFs from the public S3 bucket (``wget -nc`` is a
   no-op if files exist), so a script can be run in isolation without a
   separate setup step.
2. Run PISM with the configured stress balance and hydrology options.
3. Post-process the output: zero out missing values, rename SSA/Blatter
   velocity components to ``u_observed``/``v_observed`` (matching the
   inversion's expected names), strip _FillValue, and copy
   ``pism_config`` into the renamed file.
"""

import argparse
import os


# --- CLI ------------------------------------------------------------------

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("mode", nargs="?", default="debug",
                    choices=["debug", "sbatch"])
parser.add_argument("--start", default="1980-01-01",
                    help="simulation start date (default: %(default)s)")
parser.add_argument("--end", default="2020-01-01",
                    help="simulation end date (default: %(default)s)")
parser.add_argument("--res", default="500m",
                    help="grid resolution with units, e.g. 500m (default: %(default)s)")
parser.add_argument("--rgi-id", default="RGI2000-v7.0-C-01-04374",
                    help="RGI v7 glacier complex ID (default: %(default)s)")
parser.add_argument("--solvers", default="ssa,hybrid,blatter",
                    help="comma-separated list of stress balances to generate "
                         "(default: %(default)s)")
parser.add_argument("--hydrologies", default="null,routing",
                    help="comma-separated list of hydrologies to generate "
                         "(default: %(default)s)")
args = parser.parse_args()


# --- Constants -----------------------------------------------------------

MODE       = args.mode
START      = args.start
END        = args.end
RES        = args.res
RGI_ID     = args.rgi_id
SOLVERS    = [s.strip() for s in args.solvers.split(",") if s.strip()]
HYDROS     = [h.strip() for h in args.hydrologies.split(",") if h.strip()]

NP         = int(os.environ.get("NP", 24 if MODE == "sbatch" else 8))
SCRIPTDIR  = "run_prep_scripts"
S3_BUCKET  = "https://pism-cloud-data.s3.amazonaws.com/inverse"

GRID_FILE  = f"grid_{RGI_ID}.nc"
BOOT_FILE  = f"bootfile_{RGI_ID}.nc"
OBS_FILE   = f"obs_{RGI_ID}_0.nc"
CLIM_FILE  = f"era5_wgs84_{RGI_ID}.nc"
INPUT_FILES = [GRID_FILE, BOOT_FILE, OBS_FILE, CLIM_FILE]


SBATCH_HEADER = """\
#!/bin/sh
#SBATCH --partition=t2small
#SBATCH --ntasks=24
#SBATCH --tasks-per-node=24
#SBATCH --time=4:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=pism.%j
"""

DEBUG_HEADER = """\
#!/bin/bash
set -ex
"""

HEADER  = SBATCH_HEADER if MODE == "sbatch" else DEBUG_HEADER
RUN_CMD = f"mpirun -np {NP}"


# --- Physics flag bundles (one per stress balance / hydrology) -----------

BLATTER_PHYSICS = [
    "-stress_balance.model", "blatter",
    "-bp_ksp_monitor", "",
    "-bp_ksp_rtol", "0.001",
    "-bp_ksp_view_singularvalues", "",
    "-bp_mg_coarse_ksp_type", "preonly",
    "-bp_mg_coarse_pc_type", "lu",
    "-bp_mg_levels_ksp_type", "richardson",
    "-bp_mg_levels_pc_type", "sor",
    "-bp_pc_mg_levels", "3",
    "-bp_pc_type", "mg",
    "-bp_snes_ksp_ew", "1",
    "-bp_snes_ksp_ew_version", "3",
    "-bp_snes_monitor_ratio", "",
    "-bp_snes_rtol", "0.001",
    "-stress_balance.blatter.Mz", "10",
    "-stress_balance.blatter.coarsening_factor", "3",
    "-stress_balance.blatter.use_eta_transform", "yes",
    "-stress_balance.calving_front_stress_bc", "yes",
    "-stress_balance.blatter.enhancement_factor", "2.0",
    "-stress_balance.blatter.flow_law", "gpbld",
    "-time_stepping.adaptive_ratio", "500",
]

HYBRID_PHYSICS = [
    "-stress_balance.model", "ssa+sia",
    "-stress_balance.ssa.method", "fd",
    "-stress_balance.sia.enhancement_factor", "2.0",
    "-stress_balance.sia.flow_law", "gpbld",
    "-stress_balance.sia.max_diffusivity", "100000.0",
    "-stress_balance.sia.surface_gradient_method", "eta",
    "-stress_balance.ssa.flow_law", "gpbld",
]

SSA_PHYSICS = [
    "-stress_balance.model", "ssa",
    "-stress_balance.ssa.method", "fd",
    "-stress_balance.ssa.flow_law", "gpbld",
]

NULL_HYDRO = [
    "-hydrology.model", "null",
    "-hydrology.null_diffuse_till_water", "yes",
]

ROUTING_HYDRO = [
    "-hydrology.model", "routing",
    "-basal_yield_stress.add_transportable_water", "yes",
]


# Solver- and hydrology-keyed lookup tables, plus the per-solver rename
# needed to map SSA/Blatter outputs to the (u_observed, v_observed) names
# the inversion expects.
SOLVERS_TABLE = {
    "ssa":     {"physics": SSA_PHYSICS,
                "rename":  "-v u_ssa,u_observed -v v_ssa,v_observed"},
    "hybrid":  {"physics": HYBRID_PHYSICS,
                "rename":  "-v u_ssa,u_observed -v v_ssa,v_observed"},
    "blatter": {"physics": BLATTER_PHYSICS,
                "rename":  "-v uvelsurf,u_observed -v vvelsurf,v_observed"},
}

HYDRO_TABLE = {
    "null":    NULL_HYDRO,
    "routing": ROUTING_HYDRO,
}


# Common physics applied to every run. Mirrors the body of the original
# prepare_initial_state.sh: forcing, basal physics, geometry, IO, etc.
COMMON_PHYSICS = [
    "-atmosphere.given.file", CLIM_FILE,
    "-atmosphere.models", "given",
    "-basal_resistance.pseudo_plastic.enabled", "yes",
    "-basal_resistance.pseudo_plastic.q", "0.75",
    "-basal_resistance.pseudo_plastic.u_threshold", "100m/yr",
    "-basal_yield_stress.model", "mohr_coulomb",
    "-basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden", "0.025",
    "-basal_yield_stress.mohr_coulomb.till_phi_default", "30",
    "-calving.methods", "float_kill",
    "-energy.model", "enthalpy",
    "-geometry.front_retreat.use_cfl", "yes",
    "-geometry.part_grid.enabled", "yes",
    "-geometry.remove_icebergs", "yes",
    "-grid.Lbz", "0",
    "-grid.Lz", "2000",
    "-grid.Mbz", "1",
    "-grid.Mz", "101",
    "-grid.dx", RES,
    "-grid.dy", RES,
    "-grid.file", GRID_FILE,
    "-stress_balance.sia.bed_smoother.range", RES,
    "-grid.registration", "center",
    "-input.bootstrap", "yes",
    "-input.file", BOOT_FILE,
    "-input.forcing.buffer_size", "390",
    "-input.forcing.time_extrapolation", "yes",
    "-surface.force_to_thickness.file", BOOT_FILE,
    "-surface.force_to_thickness.alpha", "0.99",
    "-surface.force_to_thickness.ice_free_alpha_factor", "10",
    "-surface.models", "pdd,forcing",
    "-time.calendar", "standard",
    "-time.end", END,
    "-time.reference_date", START,
    "-time.start", START,
    "-time_stepping.skip.enabled", "yes",
    "-time_stepping.skip.max", "100",
]

# Output-related flags (kept separate so the spatial-vars list is easy to find
# and the per-script output filenames can be substituted in).
OUTPUT_VARS_SPATIAL = "velsurf_mag,velbase_mag,usurf,thk,velsurf,hardav,bwp,bwat"
OUTPUT_VARS_MEDIUM  = "sftgif,velsurf_mag,mask,usurf,velbase_mag,velsurf,hardav"


def output_flags(ofile, sfile):
    return [
        "-output.file", ofile,
        "-output.spatial.file", sfile,
        "-output.spatial.vars", OUTPUT_VARS_SPATIAL,
        "-output.spatial.times", "monthly",
        "-output.size", "medium",
        "-output.sizes.medium", OUTPUT_VARS_MEDIUM,
        "-output.spatial.stop_missing", "no",
    ]


# --- Emit scripts --------------------------------------------------------

os.makedirs(SCRIPTDIR, exist_ok=True)

DATA_FETCH = "\n".join(
    f"wget -nc {S3_BUCKET}/{f}" for f in INPUT_FILES
) + "\n"


count = 0
for hydro in HYDROS:
    if hydro not in HYDRO_TABLE:
        raise ValueError(f"unknown hydrology '{hydro}'; choose from "
                         f"{list(HYDRO_TABLE)}")
    for sb in SOLVERS:
        if sb not in SOLVERS_TABLE:
            raise ValueError(f"unknown stress balance '{sb}'; choose from "
                             f"{list(SOLVERS_TABLE)}")

        postfix = f"{sb}_{hydro}_g{RES}_{RGI_ID}_id_0_{START}_{END}"
        ofile   = f"state_{postfix}.nc"
        sfile   = f"spatial_{postfix}.nc"
        ofile_0 = f"state_{sb}_g{RES}_{RGI_ID}_id_0_{START}_{END}_0.nc"
        jobname = f"prep_{postfix}"
        runscript = os.path.join(SCRIPTDIR, f"{jobname}.sh")

        cmd_parts = [
            RUN_CMD, "pism",
            *SOLVERS_TABLE[sb]["physics"],
            *HYDRO_TABLE[hydro],
            *COMMON_PHYSICS,
            *output_flags(ofile, sfile),
        ]
        # Filter empty strings (PETSc-style flags with no value, like
        # -bp_ksp_monitor) so the cmdline assembles cleanly.
        run_line = " ".join(p for p in cmd_parts if p != "")

        rename     = SOLVERS_TABLE[sb]["rename"]
        attributes = "-a _FillValue,u_observed,d,, -a _FillValue,v_observed,d,,"
        post = (
            f"\n# --- Post-process: produce the inversion-friendly file ---\n"
            f"cdo setmisstoc,0 {ofile} {ofile_0}\n"
            f"ncrename {rename} {ofile_0}\n"
            f"ncatted {attributes} {ofile_0}\n"
            f"ncks -A -v pism_config {ofile} {ofile_0}\n"
        )

        with open(runscript, "w") as f:
            f.write(HEADER)
            if MODE == "sbatch":
                f.write(f"#SBATCH --job-name={jobname}\n")
            f.write("\n")
            f.write("# --- Fetch inputs (no-op if already present) ---\n")
            f.write(DATA_FETCH)
            f.write("\n")
            f.write(run_line + "\n")
            f.write(post)

        os.chmod(runscript, 0o755)
        count += 1


print(f"\nGenerated {count} scripts in {SCRIPTDIR}/\n")

if MODE == "sbatch":
    print(f"Submit all with:\n  for f in {SCRIPTDIR}/*.sh; do sbatch $f; done")
else:
    print(f"Run all with:\n  for f in {SCRIPTDIR}/*.sh; do bash $f; done")
    print(f"\nOr run one:\n  bash {SCRIPTDIR}/<script>.sh")
