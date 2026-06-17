#!/usr/bin/env python
"""Benchmark a self-contained terrain horizon-map ray-march two ways.

Both backends reproduce the math of solshade's ``compute_horizon_map`` (per-pixel
maximum terrain elevation angle in ``n_directions`` azimuths) and are invoked
through :func:`xarray.apply_ufunc`:

  * **c**     -- pure C, OpenMP-parallel, compiled to ``libsolstis.dylib`` and
                called via ctypes (``solstis.c`` / ``build.sh``).
  * **metal** -- Apple Metal via ``macmetalpy.RawKernel`` (one thread per pixel).

The script cross-checks both against a scipy reference and benchmarks them.

Usage
-----
    python bench_horizon.py --check        # correctness vs scipy on a small crop
    python bench_horizon.py                # benchmark with defaults
    python bench_horizon.py --crop 0 --n-directions 360 --max-distance 5000

Run from anywhere; paths are resolved relative to this file.
"""
from __future__ import annotations

import argparse
import ctypes
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import xarray as xr

# The C dylib and macmetalpy each link their own libomp; loading both in one
# process trips OpenMP's "multiple runtimes" abort. They don't actually share
# OpenMP state here, so allow the duplicate. Must be set before macmetalpy import.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

# Backends and helpers live in the solstis module (this is the CLI/benchmark
# front-end). The build paths (DYLIB/CSRC/BUILD) are defined there too.
from solstis import (
    load_dem,
    precompute,
    _ensure_built,
    c_horizon,
    ref_horizon,
    make_ufunc,
)
from solstis_metal import metal_horizon, metal_upload

DEM_PATH = "bootfile_RGI2000-v7.0-C-01-04374.nc"

# Authoritative reference: solshade's own compute_horizon_map (the implementation
# we reproduce). Optional -- the cross-check falls back to the scipy reimpl in
# solstis.ref_horizon when solshade isn't installed.
try:
    from solshade.terrain import compute_horizon_map
except Exception:
    compute_horizon_map = None


def solshade_horizon(dem_da, n_dir, max_distance, step):
    """Run solshade's compute_horizon_map on a DataArray; return (n_dir, H, W)."""
    hm = compute_horizon_map(
        dem_da, n_directions=n_dir, max_distance=max_distance, step=step,
        n_jobs=1, progress=False,
    )
    return np.ascontiguousarray(hm.values, dtype=np.float32)


# --------------------------------------------------------------------------- #
# Comparison + timing helpers
# --------------------------------------------------------------------------- #
def maxdiff(a, b):
    """NaN-aware: (max |Δ| over jointly-finite, # NaN-mask mismatches)."""
    a = np.asarray(a, np.float64)
    b = np.asarray(b, np.float64)
    nan_a, nan_b = np.isnan(a), np.isnan(b)
    mism = int(np.count_nonzero(nan_a ^ nan_b))
    both = ~(nan_a | nan_b)
    md = float(np.max(np.abs(a[both] - b[both]))) if np.any(both) else 0.0
    return md, mism


def compare_to_ref(got, ref):
    """Compare a backend result against a reference, NaN-asymmetry-aware.

    Returns (max |Δ| over jointly-finite, kernel_only_nan, ref_only_nan):
      * kernel_only_nan -- pixels where we are NaN but the reference has a value
                           (a real miss; should be 0).
      * ref_only_nan    -- pixels where the reference is NaN but we have a value.
                           Benign: solshade samples in float64 without snapping
                           cardinal-azimuth jitter, so its edge rays go a hair
                           out of bounds and NaN where our float32 kernels keep
                           the (correct) value. These sit on the domain border.
    """
    a = np.asarray(got, np.float64)
    b = np.asarray(ref, np.float64)
    na, nb = np.isnan(a), np.isnan(b)
    both = ~(na | nb)
    md = float(np.max(np.abs(a[both] - b[both]))) if np.any(both) else 0.0
    return md, int(np.count_nonzero(na & ~nb)), int(np.count_nonzero(nb & ~na))


def timeit(fn, warmup, runs):
    for _ in range(warmup):
        fn()
    ts = []
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return np.array(ts)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def center_crop(da, n):
    if not n:
        return da
    H, W = da.shape
    y0 = max(0, (H - n) // 2)
    x0 = max(0, (W - n) // 2)
    return da.isel(y=slice(y0, y0 + n), x=slice(x0, x0 + n))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-directions", type=int, default=64)
    p.add_argument("--max-distance", type=float, default=2000.0)
    p.add_argument("--step", type=float, default=20.0)
    p.add_argument("--crop", type=int, default=256, help="center crop NxN; 0 = full DEM")
    p.add_argument("--runs", type=int, default=5)
    p.add_argument("--warmup", type=int, default=1)
    p.add_argument("--backends", default="c,metal")
    p.add_argument("--atol", type=float, default=0.1, help="cross-check tolerance (deg)")
    p.add_argument("--rebuild", action="store_true")
    p.add_argument("--check", action="store_true", help="run correctness cross-check and exit")
    p.add_argument("--metal-exclude-transfer", action="store_true",
                   help="time only the Metal kernel dispatch (pre-upload inputs)")
    p.add_argument("--metal-batch-mb", type=int, default=512,
                   help="cap per-tile Metal output buffer (dirs are tiled to fit)")
    p.add_argument("--dem", default=DEM_PATH)
    args = p.parse_args(argv)

    backends = [b.strip() for b in args.backends.split(",") if b.strip()]
    if "c" in backends:
        _ensure_built(rebuild=args.rebuild)

    dem_da = load_dem(args.dem)
    res_x = dem_da.rio.transform().a
    res_y = -dem_da.rio.transform().e

    # ----- correctness cross-check -----
    if args.check:
        n = args.crop if args.crop else 128
        n = min(n, 128)
        sub = center_crop(dem_da, n)
        n_dir, maxd = 16, 2000.0
        az, distances, dx_pix, dy_pix = precompute(res_x, res_y, n_dir, maxd, args.step)
        dem_np = np.ascontiguousarray(sub.values, np.float32)
        print(f"[check] crop={sub.shape} n_dir={n_dir} max_distance={maxd} ns={len(distances)}")

        if compute_horizon_map is not None:
            ref = solshade_horizon(sub, n_dir, maxd, args.step)
            ref_name = "solshade"
        else:
            ref = ref_horizon(dem_np, dx_pix, dy_pix, distances)
            ref_name = "scipy-ref"
        print(f"[check] reference: {ref_name}")

        results = {}
        if "c" in backends:
            results["c"] = c_horizon(dem_np, dx_pix, dy_pix, distances)
        if "metal" in backends:
            results["metal"] = metal_horizon(dem_np, dx_pix, dy_pix, distances)

        print(f"{'pair':<18}{'max|Δ| (deg)':>14}{'kernel-NaN':>12}{'ref-NaN':>10}")
        ok = True
        for name in ("c", "metal"):
            if name in results:
                md, kn, rn = compare_to_ref(results[name], ref)
                print(f"{name+' vs '+ref_name:<18}{md:>14.6f}{kn:>12d}{rn:>10d}")
                # Fail only on real misses: angle error, or kernel NaN where the
                # reference has a value. ref-only NaN is the benign edge artifact.
                ok = ok and md <= args.atol and kn == 0
        if "c" in results and "metal" in results:
            md, mism = maxdiff(results["c"], results["metal"])
            print(f"{'c vs metal':<18}{md:>14.6f}{'(NaN mism: '+str(mism)+')':>22}")
        print("[check]", "PASS" if ok else f"FAIL (atol={args.atol} deg)")
        print("  ref-NaN = pixels the reference NaNs at the domain edge (cardinal "
              "azimuths,\n  float64 vs our float32) but our kernels resolve; benign.")
        return 0 if ok else 1

    # ----- benchmark -----
    sub = center_crop(dem_da, args.crop)
    H, W = sub.shape
    az, distances, dx_pix, dy_pix = precompute(
        res_x, res_y, args.n_directions, args.max_distance, args.step)
    ns = len(distances)
    n_dir = args.n_directions
    dem_np = np.ascontiguousarray(sub.values, np.float32)
    samples = H * W * n_dir * (ns - 1)
    out_gb = n_dir * H * W * 4 / 1e9

    nthreads = os.environ.get("OMP_NUM_THREADS", "all")
    print(f"DEM crop: {H}x{W} ({H*W:,} px) | n_dir={n_dir} ns={ns} "
          f"max_distance={args.max_distance} step={args.step}")
    print(f"ray samples: {samples:,.0f} | out buffer: {out_gb:.2f} GB | "
          f"OMP_NUM_THREADS={nthreads}")
    print(f"runs={args.runs} warmup={args.warmup}\n")

    funcs = {}
    if "c" in backends:
        funcs["c"] = lambda: c_horizon(dem_np, dx_pix, dy_pix, distances)
    if "metal" in backends:
        bmb = args.metal_batch_mb
        if args.metal_exclude_transfer:
            pre = metal_upload(dem_np, distances)
            funcs["metal"] = lambda: metal_horizon(dem_np, dx_pix, dy_pix, distances,
                                                   batch_mb=bmb, pre=pre)
        else:
            funcs["metal"] = lambda: metal_horizon(dem_np, dx_pix, dy_pix, distances,
                                                   batch_mb=bmb)

    results = {name: fn() for name, fn in funcs.items()}  # one eager run for agreement
    meds = {}
    print(f"{'backend':<10}{'median s':>12}{'pixels/s':>16}{'samples/s':>16}")
    for name, fn in funcs.items():
        ts = timeit(fn, args.warmup, args.runs)
        med = float(np.median(ts))
        meds[name] = med
        print(f"{name:<10}{med:>12.4f}{H*W/med:>16,.0f}{samples/med:>16.3e}")

    if "c" in meds and "metal" in meds:
        print(f"\nMetal ÷ C speedup: {meds['c']/meds['metal']:.2f}x"
              + (" (kernel only)" if args.metal_exclude_transfer
                 else " (incl. H2D+D2H)"))
        md, mism = maxdiff(results["c"], results["metal"])
        print(f"agreement (c vs metal): max|Δ|={md:.6f} deg, NaN-mismatch={mism}")

    # apply_ufunc overhead (the required entry point), measured once
    print("\n[apply_ufunc] wrapped call (xarray entry point):")
    for name, fn in funcs.items():
        if name == "c":
            backend_func = c_horizon
        else:
            backend_func = lambda d, dx, dy, dist: metal_horizon(
                d, dx, dy, dist, batch_mb=args.metal_batch_mb)
        call = make_ufunc(backend_func, az, dx_pix, dy_pix, distances)
        t0 = time.perf_counter()
        da = call(sub)
        dt = time.perf_counter() - t0
        finite = np.isfinite(da.values)
        lo = float(np.nanmin(da.values)) if finite.any() else float("nan")
        hi = float(np.nanmax(da.values)) if finite.any() else float("nan")
        print(f"  {name:<8} {dt:.4f}s -> DataArray{tuple(da.dims)} {da.shape}, "
              f"range [{lo:.2f}, {hi:.2f}] deg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
