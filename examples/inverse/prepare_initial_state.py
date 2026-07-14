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

from add_latlon import ensure_obs_latlon

# Path to add_latlon.py invoked as a CLI from generated shell scripts: it
# adds CF lat/lon bounds to the PISM target state file at runtime, since
# CDO's remapycon needs cell-corner coords on both source and target.
_ADD_LATLON_CLI = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "add_latlon.py"))


# --- CLI ------------------------------------------------------------------

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("mode", nargs="?", default="debug",
                    choices=["debug", "sbatch"])
parser.add_argument("--start", default="1980-01-01",
                    help="simulation start date (default: %(default)s)")
parser.add_argument("--end", default="1990-01-01",
                    help="simulation end date (default: %(default)s)")
parser.add_argument("--resolution", default="500m",
                    help="grid resolution with units, e.g. 500m (default: %(default)s)")
parser.add_argument("--rgi-id", default="RGI2000-v7.0-C-01-04374",
                    help="RGI v7 glacier complex ID (default: %(default)s)")
parser.add_argument("--solvers", default="ssa,hybrid,blatter",
                    help="comma-separated list of stress balances to generate "
                         "(default: %(default)s)")
parser.add_argument("--ice-thickness", default="frank,maffezzoli",
                    help="comma-separated list of ice thickness datasets to generate "
                         "(default: %(default)s)")
parser.add_argument("--hydrologies", default="null,routing,diffuse",
                    help="comma-separated list of hydrologies to generate "
                         "(default: %(default)s)")
parser.add_argument("--output-dir", default=".",
                    help="directory where PISM output (state, spatial, and "
                         "postprocessed _0 files) is written. The generated "
                         "scripts mkdir -p this directory at runtime. "
                         "(default: %(default)s)")
args = parser.parse_args()


# --- Constants -----------------------------------------------------------

MODE       = args.mode
START      = args.start
END        = args.end
ICE_THICKNESS = [i.strip() for i in args.ice_thickness.split(",") if i.strip()]
RES        = args.resolution
RGI_ID     = args.rgi_id
SOLVERS    = [s.strip() for s in args.solvers.split(",") if s.strip()]
HYDROS     = [h.strip() for h in args.hydrologies.split(",") if h.strip()]
OUTPUT_DIR = args.output_dir.rstrip("/")  # strip trailing slash for clean joins

NP         = int(os.environ.get("NP", 24 if MODE == "sbatch" else 8))
SCRIPTDIR  = "run_prep_scripts"
S3_BUCKET  = "https://pism-cloud-data.s3.amazonaws.com/inverse"

GRID_FILE  = f"grid_{RGI_ID}.nc"
OBS_FILE   = f"obs_{RGI_ID}.nc"
OBS0_FILE   = f"obs_{RGI_ID}_0.nc"
CLIM_FILE  = f"era5_wgs84_{RGI_ID}.nc"

v_max = 1000.0

def bootfile_for(thk):
    """Per-ice-thickness bootstrap file name."""
    return f"bootfile_{thk}_{RGI_ID}.nc"


def input_files_for(thk):
    """All input NetCDFs needed for one (thk) run, in download order."""
    return [GRID_FILE, bootfile_for(thk), OBS_FILE, CLIM_FILE]


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
    "-bp_ksp_type", "fgmres",
    "-bp_ksp_rtol", "0.001",
    "-bp_mg_coarse_ksp_type", "gmres",
    "-bp_mg_coarse_pc_type", "gamg",
    "-bp_mg_coarse_ksp_rtol", "1e-2",
    "-bp_mg_coarse_ksp_max_it", "50",
    "-bp_mg_levels_ksp_type", "richardson",
    "-bp_mg_levels_pc_type", "sor",
    "-bp_pc_mg_levels", "3",
    "-bp_pc_type", "mg",
    "-bp_snes_ksp_ew", "1",
    "-bp_snes_ksp_ew_version", "3",
    "-bp_snes_monitor", "",
    "-bp_snes_rtol", "0.001",
    "-stress_balance.blatter.Mz", "10",
    "-stress_balance.blatter.coarsening_factor", "3",
    "-stress_balance.blatter.use_eta_transform", "yes",
    "-stress_balance.calving_front_stress_bc", "yes",
    "-stress_balance.blatter.flow_law", "isothermal_glen",
    "-time_stepping.adaptive_ratio", "500",
]

HYBRID_PHYSICS = [
    "-stress_balance.model", "ssa+sia",
    "-stress_balance.ssa.method", "fd",
    "-stress_balance.sia.flow_law", "isothermal_glen",
    "-stress_balance.sia.surface_gradient_method", "eta",
    "-stress_balance.ssa.flow_law", "isothermal_glen",
]

SSA_PHYSICS = [
    "-stress_balance.model", "ssa",
    "-stress_balance.ssa.method", "fd",
    "-stress_balance.ssa.flow_law", "isothermal_glen",
]

NULL_HYDRO = [
    "-hydrology.model", "null",
    "-hydrology.null_diffuse_till_water", "no",
]

DIFFUSE_HYDRO = [
    "-hydrology.model", "null",
    "-hydrology.null_diffuse_till_water", "yes",
]

ROUTING_HYDRO = [
    "-hydrology.model", "routing",
    "-basal_yield_stress.add_transportable_water", "yes",
    "-basal_yield_stress.mohr_coulomb.till_log_factor_transportable_water", "10.0",
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
    "diffuse": DIFFUSE_HYDRO,
    "routing": ROUTING_HYDRO,
}


# Common physics applied to every run. Mirrors the body of the original
# prepare_initial_state.sh: forcing, basal physics, geometry, IO, etc.
# Built per-iteration because the bootstrap file (input.file +
# force_to_thickness.file) depends on the current ice-thickness dataset.
def common_physics_for(thk):
    boot = bootfile_for(thk)
    return [
        "-atmosphere.given.file", CLIM_FILE,
        "-atmosphere.models", "given",
        "-basal_resistance.pseudo_plastic.enabled", "yes",
        "-basal_resistance.pseudo_plastic.q", "0.75",
        "-basal_resistance.pseudo_plastic.u_threshold", "100m/yr",
        "-basal_yield_stress.model", "mohr_coulomb",
        "-basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden", "0.025",
        "-basal_yield_stress.mohr_coulomb.till_phi_default", "35",
        "-calving.methods", "float_kill",
        "-energy.model", "none",
        "-geometry.front_retreat.use_cfl", "yes",
        "-geometry.part_grid.enabled", "yes",
        "-geometry.remove_icebergs", "yes",
        "-grid.Lbz", "0",
        "-grid.Lz", "1200",
        "-grid.Mbz", "1",
        "-grid.Mz", "61",
        "-grid.dx", RES,
        "-grid.dy", RES,
        "-grid.file", GRID_FILE,
        "-stress_balance.sia.bed_smoother.range", str(1000),
        "-grid.registration", "center",
        "-input.bootstrap", "yes",
        "-input.file", boot,
        "-input.forcing.buffer_size", "390",
        "-input.forcing.time_extrapolation", "yes",
        "-surface.force_to_thickness.file", boot,
        "-surface.force_to_thickness.alpha", "0.5",
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
OUTPUT_VARS_MEDIUM  = "sftgif,velsurf_mag,mask,usurf,velbase_mag,velsurf,hardav,velsurf"


def output_flags(ofile, sfile):
    return [
        "-output.file", ofile,
        "-output.spatial.file", sfile,
        "-output.spatial.vars", OUTPUT_VARS_SPATIAL,
        "-output.spatial.times", "yearly",
        "-output.size", "medium",
        "-output.sizes.medium", OUTPUT_VARS_MEDIUM,
        "-output.spatial.stop_missing", "no",

    ]


# --- Emit scripts --------------------------------------------------------

os.makedirs(SCRIPTDIR, exist_ok=True)


def data_fetch_for(thk):
    """wget block for the inputs needed by a single (thk) run."""
    return "\n".join(
        f"wget -nc {S3_BUCKET}/{f}" for f in input_files_for(thk)
    ) + "\n"


count = 0
for thk in ICE_THICKNESS:
    for hydro in HYDROS:
        if hydro not in HYDRO_TABLE:
            raise ValueError(f"unknown hydrology '{hydro}'; choose from "
                             f"{list(HYDRO_TABLE)}")
        for sb in SOLVERS:
            if sb not in SOLVERS_TABLE:
                raise ValueError(f"unknown stress balance '{sb}'; choose from "
                                 f"{list(SOLVERS_TABLE)}")

            postfix = f"{sb}_{thk}_{hydro}_g{RES}_{RGI_ID}_id_0_{START}_{END}"
            # Output filenames (paths). When OUTPUT_DIR == "." they reduce to
            # plain filenames in the cwd; otherwise they include the directory.
            ofile   = f"{OUTPUT_DIR}/state_{postfix}.nc"
            ofile_0   = f"{OUTPUT_DIR}/state_{postfix}_0.nc"
            sfile   = f"{OUTPUT_DIR}/spatial_{postfix}.nc"
            jobname = f"prep_{postfix}"
            runscript = os.path.join(SCRIPTDIR, f"{jobname}.sh")

            cmd_parts = [
                RUN_CMD, "pism",
                *SOLVERS_TABLE[sb]["physics"],
                *HYDRO_TABLE[hydro],
                *common_physics_for(thk),
                *output_flags(ofile, sfile),
            ]
            # Filter empty strings (PETSc-style flags with no value, like
            # -bp_ksp_monitor) so the cmdline assembles cleanly.
            run_line = " ".join(p for p in cmd_parts if p != "")

            OBS_LATLON_FILE = OBS_FILE.replace(".nc", "_latlon.nc")
            ensure_obs_latlon(OBS_FILE, OBS_LATLON_FILE)
            ofile_latlon = ofile.replace(".nc", "_latlon.nc")

            post = (
                f"\n# --- Post-process: produce the inversion-friendly file ---\n"
                f"cdo -O setmisstoc,0 {ofile} {ofile_0}\n"
                f"ncks -A -v pism_config {ofile} {ofile_0}\n"
            )
            post_misfit = (
                f"\n# --- Post-process: adjust misfit_weight ---\n"
                f"python {_ADD_LATLON_CLI} {ofile} {ofile_latlon}\n"
                f"cdo remapycon,{ofile_latlon} {OBS_LATLON_FILE} {OBS0_FILE}\n"
                f"ncks -A -v velbase_mag {ofile} {OBS0_FILE}\n"
                f"""ncap2 -O -s "vel_misfit_weight=float(vel_misfit_weight); where(velbase_mag > {v_max}) vel_misfit_weight=0.1" {OBS0_FILE} {OBS0_FILE}\n"""
                )
            with open(runscript, "w") as f:
                f.write(HEADER)
                if MODE == "sbatch":
                    f.write(f"#SBATCH --job-name={jobname}\n")
                f.write("\n")
                f.write("# --- Fetch inputs (no-op if already present) ---\n")
                f.write(data_fetch_for(thk))
                f.write("\n")
                f.write(f"# --- Ensure output directory exists ---\n")
                f.write(f"mkdir -p {OUTPUT_DIR}\n\n")
                f.write(run_line + "\n")
                f.write(post + "\n")
                f.write(post_misfit + "\n")

            os.chmod(runscript, 0o755)
            count += 1


print(f"\nGenerated {count} scripts in {SCRIPTDIR}/\n")

if MODE == "sbatch":
    print(f"Submit all with:\n  for f in {SCRIPTDIR}/*.sh; do sbatch $f; done")
else:
    print(f"Run all with:\n  for f in {SCRIPTDIR}/*.sh; do bash $f; done")
    print(f"\nOr run one:\n  bash {SCRIPTDIR}/<script>.sh")
