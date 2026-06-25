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
#   NP=8 ./bench_blatter_solvers_local.sh                  # more ranks
#   TEST=A L_KM=10 ./bench_blatter_solvers_local.sh        # a different experiment
#   ONLY="baseline_lu gamg" ./bench_blatter_solvers_local.sh   # subset
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
TEST="${TEST:-C}"         # ISMIP-HOM experiment (C has sliding -> nonzero BP velocity)
L_KM="${L_KM:-80}"        # domain side length, km
OUTDIR="${OUTDIR:-bench_$(date +%Y%m%d_%H%M%S)}"
ONLY="${ONLY:-}"          # space-separated subset of config names; empty = all

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
  -bp_pc_mg_levels 3
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
  echo "  $name   (np=$NP, ISMIP-HOM $TEST, L=${L_KM}km, ${MX}x${MY})"
  echo "==========================================================="

  local start end wall
  start=$(date +%s.%N)
  # shellcheck disable=SC2086
  BENCH_TEST="$TEST" BENCH_L_KM="$L_KM" BENCH_MX="$MX" BENCH_MY="$MY" BENCH_OUT="$out" \
  mpirun -np "$NP"  python3 "$DRIVER" \
      "${COMMON[@]}" $(opts_for "$name") > "$log" 2>&1
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
