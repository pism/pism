#!/usr/bin/env bash
#
# Benchmark Blatter/BP coarse-grid solver options.
#
# Each configuration runs the SAME single ISMIP-HOM Blatter solve (via
# bench_ismiphom_solve.py) and records the solve time. ISMIP-HOM produces real
# 3-D velocities -- experiment C has spatially varying basal sliding -- so the
# BP linear/coarse solve is actually exercised (unlike a frozen-bed dome).
# Afterwards we compare every run's surface velocity against the baseline so you
# can confirm the faster solvers converge to the same answer.
#
# Usage:
#   ./bench_blatter_solvers_local.sh                       # all configs
#   MX=401 MY=401 ./bench_blatter_solvers_local.sh         # bigger coarse grid
#   MZ=10 CF=3 ./bench_blatter_solvers_local.sh            # match Greenland production (10->4->2)
#   NP=8 ./bench_blatter_solvers_local.sh                  # more ranks
#   TEST=A L_KM=10 ./bench_blatter_solvers_local.sh        # a different experiment
#   ONLY="baseline_lu gamg" ./bench_blatter_solvers_local.sh   # subset
#
# MZ / CF / MG_LEVELS set the Blatter vertical grid and MG depth (they mirror the
# production stress_balance.blatter.Mz, .coarsening_factor, and -bp_pc_mg_levels).
# They must be compatible: CF^(MG_LEVELS-1) divides (Mz-1), coarsest grid >= 2.
# Default 9/2/3 (9->5->3); Greenland production is 10/3/3 (10->4->2).
#
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
DRIVER="$HERE/ismiphom_solve.py"

# ---------------------------------------------------------------------------
# Knobs (override from the environment)
# ---------------------------------------------------------------------------
NP="${NP:-4}"
MX="${MX:-201}"
MY="${MY:-201}"
MZ="${MZ:-9}"             # Blatter solver vertical levels (drives MG depth)
CF="${CF:-2}"             # vertical MG coarsening factor
MG_LEVELS="${MG_LEVELS:-3}"  # number of vertical MG levels (independent knob, like production)
TEST="${TEST:-C}"         # ISMIP-HOM experiment (C has sliding -> nonzero BP velocity)
L_KM="${L_KM:-80}"        # domain side length, km
OUTDIR="${OUTDIR:-bench_$(date +%Y%m%d_%H%M%S)}"
ONLY="${ONLY:-}"          # space-separated subset of config names; empty = all
EXTRA="${EXTRA:-}"        # extra options appended to every run, e.g. EXTRA="-log_view"

# Mz, CF, and MG_LEVELS are independent (as in production), but must be mutually
# compatible: PISM's Blatter MG semicoarsens in the vertical only, coarsening
# Mz -> (Mz-1)/CF+1 at each level, so CF^(MG_LEVELS-1) must divide (Mz-1) and the
# coarsest grid must have >= 2 points. Examples: Mz=9,CF=2,levels=3 -> 9->5->3;
# Mz=10,CF=3,levels=3 -> 10->4->2 (production). Validate before running.
{
  d=1
  for (( i = 1; i < MG_LEVELS; i++ )); do d=$(( d * CF )); done   # d = CF^(MG_LEVELS-1)
  coarsest=0
  (( (MZ - 1) % d == 0 )) && coarsest=$(( (MZ - 1) / d + 1 ))
  if (( MG_LEVELS < 2 || coarsest < 2 )); then
    echo "ERROR: Mz=$MZ, CF=$CF, MG_LEVELS=$MG_LEVELS are incompatible." >&2
    echo "       Need CF^(MG_LEVELS-1)=$d to divide (Mz-1)=$(( MZ - 1 )) with coarsest grid >= 2 points." >&2
    echo "       e.g. Mz=9 CF=2 MG_LEVELS=3 (9->5->3), or Mz=10 CF=3 MG_LEVELS=3 (10->4->2)." >&2
    exit 1
  fi
  MG_COARSEST=$coarsest
}

mkdir -p "$OUTDIR/logs"

# Solver options shared by every run. -bp_snes_monitor and the converged reason
# put iteration counts in the logs; -bp_mg_coarse_ksp_monitor shows what the
# coarse solve is doing. The grid (-Mx/-My) and Mz/coarsening are supplied per
# run / by the driver.
COMMON=(
  -bp_ksp_rtol 0.001
  -bp_snes_rtol 0.001
  # Cap the outer (linear) solve so a divergent preconditioner fails in seconds
  # instead of grinding to the default 10000 iterations.
  -bp_ksp_max_it 1000
  -bp_snes_ksp_ew 1
  -bp_snes_ksp_ew_version 3
  # Flexible outer Krylov: required because the "_ksp" coarse configs make the MG
  # preconditioner variable; plain GMRES can diverge (DIVERGED_LINEAR_SOLVE).
  -bp_ksp_type fgmres
  -bp_pc_type mg
  -bp_pc_mg_levels "$MG_LEVELS"
  -bp_mg_levels_ksp_type richardson
  -bp_mg_levels_pc_type sor
  -bp_ksp_monitor
  -bp_snes_monitor
  -bp_snes_converged_reason
  -bp_mg_coarse_ksp_monitor
)

# ---------------------------------------------------------------------------
# Configurations: name -> coarse-grid (and smoother) option string.
# Everything else is held fixed via COMMON so differences are attributable.
# ---------------------------------------------------------------------------
#
# Configs marked [MUMPS] / [hypre] need PETSc built with that package; they will
# fail cleanly (FAILED in the summary) if it is missing -- just drop them from
# ONLY in that case. Everything else uses only built-in PETSc functionality.
names=(
  baseline_lu             # direct LU (built-in serial LU; parallel needs MUMPS)
  mumps                   # explicit MUMPS parallel direct solve            [MUMPS]
  gamg                    # one V-cycle of GAMG (preonly) -- the cheap option
  gamg_ksp                # GAMG wrapped in a Krylov solve (actually converges)
  gamg_cheb_sor           # gamg + Chebyshev/SOR level smoother
  gamg_cheb_sor_ksp       # the above, with the coarse solve Krylov-wrapped
  hypre                   # BoomerAMG (hypre) on the coarse level           [hypre]
  # --- preconditioner-reuse variants (build on gamg_ksp; see opts_for notes) ---
  gamg_ksp_lag            # gamg_ksp, but build the GAMG hierarchy once per solve & reuse
  gamg_ksp_lag2           # gamg_ksp, rebuild GAMG every 2nd Newton step (persists)
  gamg_ksp_reuse          # gamg_ksp, rebuild each step but cheaper GAMG re-setup
  gamg_ksp_lag_reuse      # gamg_ksp, lag within a solve + cheaper re-setup
  hypre_ksp               # BoomerAMG wrapped in a bounded GMRES coarse solve   [hypre]
  hypre_ksp_lag           # hypre_ksp, coarse PC built once per solve & reused   [hypre]
)

# Map a config name to its coarse-grid (and smoother) option string. Implemented
# as a case statement rather than an associative array so it works with the
# macOS system bash (3.2), which predates `declare -A`.
#
# Background: PISM's Blatter multigrid semicoarsens in the vertical only, so the
# coarse level is still the full horizontal grid (large). "preonly" applies the
# coarse PC exactly once; the "_ksp" variants instead run a few Krylov iterations
# so the coarse residual is actually reduced -- usually the difference between a
# weak and a strong preconditioner on stiff (high-velocity) steps.
opts_for () {
  case "$1" in
    # Direct LU. Built-in LU is serial; in parallel PETSc needs MUMPS/SuperLU_dist.
    baseline_lu)
      echo "-bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type lu" ;;
    # Explicit parallel direct solve via MUMPS (the "proper" parallel baseline).
    mumps)
      echo "-bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type lu \
            -bp_mg_coarse_pc_factor_mat_solver_type mumps" ;;
    # PETSc-native algebraic multigrid (GAMG), one V-cycle.
    gamg)
      echo "-bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type gamg" ;;
    # GAMG as a preconditioner inside a bounded GMRES coarse solve.
    gamg_ksp)
      echo "-bp_mg_coarse_pc_type gamg \
            -bp_mg_coarse_ksp_type gmres -bp_mg_coarse_ksp_rtol 1e-2 \
            -bp_mg_coarse_ksp_max_it 50" ;;
    # GAMG coarse + a Chebyshev/SOR level smoother (overrides the richardson/sor
    # smoother from COMMON). SOR (not point-Jacobi) is used because the Blatter
    # operator is strongly anisotropic in the vertical: Chebyshev/Jacobi diverges
    # here, Chebyshev/SOR is stable.
    gamg_cheb_sor)
      echo "-bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type gamg \
            -bp_mg_levels_ksp_type chebyshev -bp_mg_levels_pc_type sor" ;;
    # Same smoother, but let the coarse GAMG run a few Krylov iterations.
    gamg_cheb_sor_ksp)
      echo "-bp_mg_coarse_pc_type gamg \
            -bp_mg_coarse_ksp_type gmres -bp_mg_coarse_ksp_rtol 1e-2 \
            -bp_mg_coarse_ksp_max_it 50 \
            -bp_mg_levels_ksp_type chebyshev -bp_mg_levels_pc_type sor" ;;
    # hypre BoomerAMG on the coarse level (often a strong, scalable AMG).
    hypre)
      echo "-bp_mg_coarse_ksp_type preonly -bp_mg_coarse_pc_type hypre \
            -bp_mg_coarse_pc_hypre_type boomeramg" ;;

    # ---- Preconditioner-reuse variants -------------------------------------
    # These build on gamg_ksp (the production baseline). Background: with Newton
    # + Eisenstat-Walker the Jacobian -- and therefore the coarse GAMG hierarchy
    # -- is rebuilt on *every* SNES iteration, and on stiff (high-velocity) steps
    # that AMG setup can dominate the solve. -snes_lag_preconditioner reuses the
    # hierarchy across Newton steps; EW only rescales the KSP *tolerance*, so it
    # is unaffected. Use -log_view (EXTRA="-log_view") to see PCSetUp vs KSPSolve.
    #
    # -bp_snes_lag_preconditioner -2 : build the PC once at the start of the
    #    SNESSolve and reuse it for the rest of that solve.
    gamg_ksp_lag)
      echo "-bp_mg_coarse_pc_type gamg \
            -bp_mg_coarse_ksp_type gmres -bp_mg_coarse_ksp_rtol 1e-2 \
            -bp_mg_coarse_ksp_max_it 50 \
            -bp_snes_lag_preconditioner -2" ;;
    # Rebuild the hierarchy every 2nd Newton step; the lag counter persists across
    # the retry-ladder / parameter-continuation SNESSolves.
    gamg_ksp_lag2)
      echo "-bp_mg_coarse_pc_type gamg \
            -bp_mg_coarse_ksp_type gmres -bp_mg_coarse_ksp_rtol 1e-2 \
            -bp_mg_coarse_ksp_max_it 50 \
            -bp_snes_lag_preconditioner 2 -bp_snes_lag_preconditioner_persists true" ;;
    # Rebuild every step (no lag), but make each GAMG (re)setup cheaper: reuse the
    # prolongation and coarsen more aggressively. On PETSc < 3.19 replace
    # -..._aggressive_coarsening 1 with -..._square_graph 1.
    gamg_ksp_reuse)
      echo "-bp_mg_coarse_pc_type gamg \
            -bp_mg_coarse_ksp_type gmres -bp_mg_coarse_ksp_rtol 1e-2 \
            -bp_mg_coarse_ksp_max_it 50 \
            -bp_mg_coarse_pc_gamg_reuse_interpolation true \
            -bp_mg_coarse_pc_gamg_threshold 0.05 \
            -bp_mg_coarse_pc_gamg_aggressive_coarsening 1" ;;
    # Combine: lag within a solve AND cheaper re-setup when it does rebuild.
    gamg_ksp_lag_reuse)
      echo "-bp_mg_coarse_pc_type gamg \
            -bp_mg_coarse_ksp_type gmres -bp_mg_coarse_ksp_rtol 1e-2 \
            -bp_mg_coarse_ksp_max_it 50 \
            -bp_snes_lag_preconditioner -2 \
            -bp_mg_coarse_pc_gamg_reuse_interpolation true" ;;
    # hypre BoomerAMG wrapped in a bounded GMRES coarse solve -- a fair AMG-vs-AMG
    # comparison against gamg_ksp.
    hypre_ksp)
      echo "-bp_mg_coarse_pc_type hypre -bp_mg_coarse_pc_hypre_type boomeramg \
            -bp_mg_coarse_ksp_type gmres -bp_mg_coarse_ksp_rtol 1e-2 \
            -bp_mg_coarse_ksp_max_it 50" ;;
    # hypre coarse PC built once per solve and reused across Newton steps.
    hypre_ksp_lag)
      echo "-bp_mg_coarse_pc_type hypre -bp_mg_coarse_pc_hypre_type boomeramg \
            -bp_mg_coarse_ksp_type gmres -bp_mg_coarse_ksp_rtol 1e-2 \
            -bp_mg_coarse_ksp_max_it 50 \
            -bp_snes_lag_preconditioner -2" ;;
    *)
      echo "unknown config '$1'" >&2; return 1 ;;
  esac
}

# ---------------------------------------------------------------------------
# Run loop
# ---------------------------------------------------------------------------
results="$OUTDIR/results.tsv"
printf "config\tstatus\tsolve_s\twall_s\toutput\n" > "$results"

run_one () {
  local name="$1"
  local log="$OUTDIR/logs/${name}.log"
  local out="$OUTDIR/out_${name}.nc"

  echo "==========================================================="
  echo "  $name   (np=$NP, ISMIP-HOM $TEST, L=${L_KM}km, ${MX}x${MY}, Mz=${MZ}, CF=${CF}, ${MG_LEVELS} MG levels -> coarsest Mz=${MG_COARSEST})"
  echo "==========================================================="

  local start end wall
  start=$(date +%s.%N)
  # shellcheck disable=SC2086
  BENCH_TEST="$TEST" BENCH_L_KM="$L_KM" BENCH_MX="$MX" BENCH_MY="$MY" \
  BENCH_MZ="$MZ" BENCH_CF="$CF" BENCH_OUT="$out" \
  mpirun -np "$NP"  python3 "$DRIVER" \
      "${COMMON[@]}" $(opts_for "$name") $EXTRA > "$log" 2>&1
  local rc=$?
  end=$(date +%s.%N)
  wall=$(awk "BEGIN{printf \"%.1f\", $end-$start}")

  # The driver prints the time spent inside the solve itself (excludes startup).
  local solve_s
  solve_s=$(grep -Eo "solve_seconds=[0-9.]+" "$log" | tail -1 | cut -d= -f2)
  [ -z "$solve_s" ] && solve_s="n/a"

  local status
  if [ $rc -eq 0 ] && [ -f "$out" ]; then
    status="OK"
  else
    status="FAILED(rc=$rc)"
    out="-"
  fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$name" "$status" "$solve_s" "$wall" "$out" >> "$results"
  echo ">>> $name: $status  solve=${solve_s}s  wall=${wall}s"
}

for name in "${names[@]}"; do
  if [ -n "$ONLY" ] && [[ " $ONLY " != *" $name "* ]]; then
    continue
  fi
  run_one "$name"
done

# ---------------------------------------------------------------------------
# Compare velocity fields against the baseline. Uses the local python3 (the
# active conda env already has netCDF4).
# ---------------------------------------------------------------------------
cat > "$OUTDIR/compare.py" <<'PY'
import sys, glob, os
import numpy as np
from netCDF4 import Dataset

def velmag(path):
    """Return (label, flat magnitude array) using the first available field."""
    ds = Dataset(path)
    v = ds.variables
    def arr(n):
        return np.ma.filled(v[n][:].astype("f8"), np.nan)
    try:
        if "velsurf_mag" in v:
            a, lab = arr("velsurf_mag"), "velsurf_mag"
        elif "uvelsurf" in v and "vvelsurf" in v:
            a, lab = np.hypot(arr("uvelsurf"), arr("vvelsurf")), "uvelsurf/vvelsurf"
        elif "velbar_mag" in v:
            a, lab = arr("velbar_mag"), "velbar_mag"
        elif "ubar" in v and "vbar" in v:
            a, lab = np.hypot(arr("ubar"), arr("vbar")), "ubar/vbar"
        elif "u_sigma" in v and "v_sigma" in v:
            a, lab = np.hypot(arr("u_sigma"), arr("v_sigma")), "u_sigma/v_sigma"
        else:
            return None, None
        return lab, np.asarray(a).ravel()
    finally:
        ds.close()

files = sorted(sys.argv[1:])
base_path = None
for f in files:
    if "out_baseline_lu" in f:
        base_path = f
if base_path is None and files:
    base_path = files[0]

blab, base = velmag(base_path)
print(f"\nReference: {os.path.basename(base_path)}  (field: {blab})")
bmask = np.isfinite(base)
bmax = np.nanmax(np.abs(base[bmask])) if bmask.any() else float("nan")
print(f"  max |v| = {bmax:.6g} m/yr\n")

hdr = f"{'config':<22}{'field':<18}{'max|v|':>12}{'max abs diff':>15}{'max rel diff':>14}"
print(hdr); print("-"*len(hdr))
for f in files:
    lab, a = velmag(f)
    if a is None:
        print(f"{os.path.basename(f):<22}  (no velocity field found)")
        continue
    m = np.isfinite(a) & bmask
    vmax = np.nanmax(np.abs(a[np.isfinite(a)]))
    if a.shape == base.shape and m.any():
        adiff = np.nanmax(np.abs(a[m]-base[m]))
        rdiff = adiff / bmax if bmax else float("nan")
        ds = f"{adiff:>15.4g}{rdiff:>14.2e}"
    else:
        ds = f"{'shape mismatch':>29}"
    tag = os.path.basename(f).replace("out_","").replace(".nc","")
    print(f"{tag:<22}{lab:<18}{vmax:>12.6g}{ds}")
print()
PY

echo
echo "==========================================================="
echo "  Timing summary ($results)"
echo "==========================================================="
column -t -s $'\t' "$results"

echo
echo "==========================================================="
echo "  Correctness check (velocity vs baseline)"
echo "==========================================================="
shopt -s nullglob
ncfiles=("$OUTDIR"/out_*.nc)
if [ ${#ncfiles[@]} -gt 0 ]; then
  python3 "$OUTDIR/compare.py" "${ncfiles[@]}"
else
  echo "No output files to compare."
fi

echo "All artifacts in: $OUTDIR/"
