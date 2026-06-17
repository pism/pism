#!/usr/bin/env python
"""Reproduce solshade's compute_slope_aspect_normals with C and Metal, and
benchmark / validate the two backends.

The single solshade function is split into three independent ops, each a local
per-pixel stencil (numpy.gradient edge_order=1, then the GIS slope/aspect/
ENU-normal math):

  * slope   (y, x)      -- angle from horizontal, degrees
  * aspect  (y, x)      -- compass direction of steepest descent, deg cw from N
  * normals (3, y, x)   -- ENU unit normal vectors, bands [east, north, up]

Each has a C (OpenMP, via ctypes) and a Metal (macmetalpy RawKernel) backend in
solstis.py. This script cross-checks both against solshade and times them.

Usage
-----
    python bench_slope.py --check     # correctness vs solshade
    python bench_slope.py             # benchmark with defaults
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import xarray as xr

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

from solstis import (
    load_dem,
    _ensure_built,
    c_slope, c_aspect, c_normals,
)
from solstis_metal import metal_slope, metal_aspect, metal_normals

DEM_PATH = "bootfile_RGI2000-v7.0-C-01-04374.nc"

# Authoritative reference: solshade's own implementation (optional).
try:
    from solshade.terrain import compute_slope_aspect_normals
except Exception:
    compute_slope_aspect_normals = None

# op name -> (C backend, Metal backend, output bands, is-circular-degrees)
OPS = {
    "slope":   (c_slope,   metal_slope,   1, True),
    "aspect":  (c_aspect,  metal_aspect,  1, True),   # circular 0..360
    "normals": (c_normals, metal_normals, 3, False),
}


def resolution(dem_da):
    dy_res, dx_res = dem_da.rio.resolution()
    return abs(dx_res), abs(dy_res)


def center_crop(da, n):
    if not n:
        return da
    H, W = da.shape
    y0 = max(0, (H - n) // 2)
    x0 = max(0, (W - n) // 2)
    return da.isel(y=slice(y0, y0 + n), x=slice(x0, x0 + n))


def angdiff(a, b):
    """Smallest absolute difference between two angles in degrees (circular)."""
    d = np.abs(np.asarray(a, np.float64) - np.asarray(b, np.float64)) % 360.0
    return np.minimum(d, 360.0 - d)


def maxdiff(got, ref, circular=False):
    a = np.asarray(got, np.float64)
    b = np.asarray(ref, np.float64)
    d = angdiff(a, b) if circular else np.abs(a - b)
    return float(np.nanmax(d))


def timeit(fn, warmup, runs):
    for _ in range(warmup):
        fn()
    ts = []
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return np.array(ts)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--crop", type=int, default=512, help="center crop NxN; 0 = full DEM")
    p.add_argument("--runs", type=int, default=5)
    p.add_argument("--warmup", type=int, default=1)
    p.add_argument("--backends", default="c,metal")
    p.add_argument("--ops", default="slope,aspect,normals")
    p.add_argument("--atol", type=float, default=1e-3, help="cross-check tolerance (deg)")
    p.add_argument("--rebuild", action="store_true")
    p.add_argument("--check", action="store_true", help="run correctness cross-check and exit")
    p.add_argument("--dem", default=DEM_PATH)
    args = p.parse_args(argv)

    backends = [b.strip() for b in args.backends.split(",") if b.strip()]
    ops = [o.strip() for o in args.ops.split(",") if o.strip()]
    if "c" in backends:
        _ensure_built(rebuild=args.rebuild)

    dem_da = load_dem(args.dem)
    dx, dy = resolution(dem_da)

    # ----- correctness cross-check -----
    if args.check:
        if compute_slope_aspect_normals is None:
            print("[check] solshade not installed; cannot cross-check.")
            return 1
        n = args.crop if args.crop else 256
        n = min(n, 256)
        sub = center_crop(dem_da, n)
        dem_np = np.ascontiguousarray(sub.values, np.float32)
        sdx, sdy = resolution(sub)
        print(f"[check] reference: solshade | crop={sub.shape} dx={sdx} dy={sdy}")

        slp, asp, nrm = compute_slope_aspect_normals(sub)
        ref = {"slope": slp.values, "aspect": asp.values, "normals": nrm.values}

        print(f"{'op / backend':<24}{'max|Δ|':>14}{'units':>10}")
        ok = True
        for op in ops:
            c_fn, m_fn, _bands, circular = OPS[op]
            unit = "deg" if op != "normals" else "unit"
            for name, fn in (("c", c_fn), ("metal", m_fn)):
                if name not in backends:
                    continue
                got = fn(dem_np, sdx, sdy)
                md = maxdiff(got, ref[op], circular=circular)
                print(f"{op+' '+name:<24}{md:>14.3e}{unit:>10}")
                ok = ok and md <= args.atol
        print("[check]", "PASS" if ok else f"FAIL (atol={args.atol})")
        return 0 if ok else 1

    # ----- benchmark -----
    sub = center_crop(dem_da, args.crop)
    H, W = sub.shape
    sdx, sdy = resolution(sub)
    dem_np = np.ascontiguousarray(sub.values, np.float32)
    nthreads = os.environ.get("OMP_NUM_THREADS", "all")
    print(f"DEM crop: {H}x{W} ({H*W:,} px) | dx={sdx} dy={sdy} | "
          f"OMP_NUM_THREADS={nthreads}")
    print(f"runs={args.runs} warmup={args.warmup}\n")

    print(f"{'op':<10}{'backend':<8}{'median s':>12}{'pixels/s':>16}{'speedup':>10}")
    for op in ops:
        c_fn, m_fn, _bands, _circ = OPS[op]
        funcs = {}
        if "c" in backends:
            funcs["c"] = lambda f=c_fn: f(dem_np, sdx, sdy)
        if "metal" in backends:
            funcs["metal"] = lambda f=m_fn: f(dem_np, sdx, sdy)
        results = {name: fn() for name, fn in funcs.items()}
        meds = {}
        for name, fn in funcs.items():
            ts = timeit(fn, args.warmup, args.runs)
            meds[name] = float(np.median(ts))
        for name in funcs:
            sp = (f"{meds['c']/meds['metal']:.1f}x"
                  if name == "metal" and "c" in meds else "")
            print(f"{op:<10}{name:<8}{meds[name]:>12.5f}{H*W/meds[name]:>16,.0f}{sp:>10}")
        if "c" in results and "metal" in results:
            _c, _m, _b, circ = OPS[op]
            md = maxdiff(results["c"], results["metal"], circular=circ)
            print(f"           agreement c vs metal: max|Δ|={md:.3e}")

    # Show the labelled xarray outputs (slope/aspect (y,x); normals (band,y,x)).
    print("\n[xarray] result DataArrays (C backend):")
    coords2d = {d: sub[d] for d in sub.dims}
    if "slope" in ops:
        s = xr.DataArray(c_slope(dem_np, sdx, sdy), coords=coords2d, dims=sub.dims,
                         attrs={"units": "degrees", "long_name": "slope"})
        print(f"  slope   {tuple(s.dims)} {s.shape}  range [{float(s.min()):.1f}, {float(s.max()):.1f}] deg")
    if "aspect" in ops:
        a = xr.DataArray(c_aspect(dem_np, sdx, sdy), coords=coords2d, dims=sub.dims,
                         attrs={"units": "degrees", "long_name": "aspect"})
        print(f"  aspect  {tuple(a.dims)} {a.shape}  range [{float(a.min()):.1f}, {float(a.max()):.1f}] deg")
    if "normals" in ops:
        nv = c_normals(dem_np, sdx, sdy)
        nrm = xr.DataArray(nv, dims=("band", *sub.dims),
                           coords={"band": ["east", "north", "up"], **coords2d},
                           attrs={"description": "Terrain normal unit vector (ENU)"})
        print(f"  normals {tuple(nrm.dims)} {nrm.shape}  |n| in "
              f"[{float(np.linalg.norm(nv, axis=0).min()):.4f}, "
              f"{float(np.linalg.norm(nv, axis=0).max()):.4f}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
