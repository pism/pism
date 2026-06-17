# Copyright (C) 2026 PISM Authors
#
# This file is part of PISM.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

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
import rioxarray

# The C dylib and macmetalpy each link their own libomp; loading both in one
# process trips OpenMP's "multiple runtimes" abort. They don't actually share
# OpenMP state here, so allow the duplicate. Must be set before macmetalpy import.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

HERE = Path(__file__).resolve().parent
DYLIB = HERE / "libsolstis.dylib"
CSRC = HERE / "solstis.c"
BUILD = HERE / "build.sh"


# DEM loading
# --------------------------------------------------------------------------- #
def load_dem(path: str, variable: str = "surface") -> xr.DataArray:
    import rioxarray  # noqa: F401  (registers the .rio accessor)

    ds = xr.open_dataset(path).squeeze()
    mapping_var = ds.rio.grid_mapping
    crs = ds[mapping_var].attrs["crs_wkt"]
    ds = ds.rio.write_crs(crs)
    return ds[variable]


# --------------------------------------------------------------------------- #
# Shared precompute: ray geometry identical for both backends
# --------------------------------------------------------------------------- #
def precompute(res_x, res_y, n_directions, max_distance, step):
    az = np.linspace(0.0, 360.0, n_directions, endpoint=False)  # float64
    distances = np.arange(0.0, max_distance + step, step)       # includes 0
    rad = np.deg2rad(az)[:, None]
    dx_pix = (np.cos(rad) * distances / res_x).astype(np.float32)   # (n_dir, ns)
    dy_pix = (-np.sin(rad) * distances / res_y).astype(np.float32)
    # Snap sub-micropixel jitter to zero. At cardinal azimuths cos/sin carry
    # ~1e-16 float noise that should be exactly 0; left in, it nudges edge rays
    # infinitesimally out of bounds and float32/float64 then disagree on the
    # boundary NaN classification. 1e-6 px = 1e-4 m is physically negligible.
    dx_pix[np.abs(dx_pix) < 1e-6] = 0.0
    dy_pix[np.abs(dy_pix) < 1e-6] = 0.0
    distances = distances.astype(np.float32)
    return az, distances, np.ascontiguousarray(dx_pix), np.ascontiguousarray(dy_pix)


# --------------------------------------------------------------------------- #
# C backend (ctypes)
# --------------------------------------------------------------------------- #
def _ensure_built(rebuild=False):
    stale = (not DYLIB.exists()) or (DYLIB.stat().st_mtime < CSRC.stat().st_mtime)
    if rebuild or stale:
        print(f"[build] compiling {CSRC.name} -> {DYLIB.name}")
        subprocess.run(["bash", str(BUILD)], check=True)


_LIB = None


def _load_lib():
    global _LIB
    if _LIB is None:
        f32 = np.ctypeslib.ndpointer(dtype=np.float32, flags="C_CONTIGUOUS")
        lib = ctypes.CDLL(str(DYLIB))
        lib.horizon_map.restype = None
        lib.horizon_map.argtypes = [
            f32, ctypes.c_int, ctypes.c_int,   # dem, H, W
            f32, f32,                          # dx_pix, dy_pix
            f32, ctypes.c_int, ctypes.c_int,   # distances, n_dir, ns
            f32,                               # out
        ]
        # slope_map / aspect_map / normals_map share one signature:
        #   (dem, H, W, dx, dy, out)
        for name in ("slope_map", "aspect_map", "normals_map"):
            fn = getattr(lib, name)
            fn.restype = None
            fn.argtypes = [
                f32, ctypes.c_int, ctypes.c_int,       # dem, H, W
                ctypes.c_float, ctypes.c_float,        # dx, dy (pixel size, m)
                f32,                                   # out
            ]
        _LIB = lib
    return _LIB


def c_horizon(dem, dx_pix, dy_pix, distances):
    lib = _load_lib()
    H, W = dem.shape
    n_dir, ns = dx_pix.shape
    dem = np.ascontiguousarray(dem, dtype=np.float32)
    out = np.empty((n_dir, H, W), dtype=np.float32)
    lib.horizon_map(
        dem, H, W,
        np.ascontiguousarray(dx_pix, np.float32),
        np.ascontiguousarray(dy_pix, np.float32),
        np.ascontiguousarray(distances, np.float32),
        n_dir, ns, out,
    )
    return out


def c_slope(dem, dx, dy):
    lib = _load_lib()
    H, W = dem.shape
    dem = np.ascontiguousarray(dem, dtype=np.float32)
    out = np.empty((H, W), dtype=np.float32)
    lib.slope_map(dem, H, W, float(dx), float(dy), out)
    return out


def c_aspect(dem, dx, dy):
    lib = _load_lib()
    H, W = dem.shape
    dem = np.ascontiguousarray(dem, dtype=np.float32)
    out = np.empty((H, W), dtype=np.float32)
    lib.aspect_map(dem, H, W, float(dx), float(dy), out)
    return out


def c_normals(dem, dx, dy):
    lib = _load_lib()
    H, W = dem.shape
    dem = np.ascontiguousarray(dem, dtype=np.float32)
    out = np.empty((3, H, W), dtype=np.float32)
    lib.normals_map(dem, H, W, float(dx), float(dy), out)
    return out


# --------------------------------------------------------------------------- #
# scipy reference (mirrors solshade's map_coordinates + reduction split)
# --------------------------------------------------------------------------- #
def ref_horizon(dem, dx_pix, dy_pix, distances):
    from scipy.ndimage import map_coordinates

    H, W = dem.shape
    n_dir, ns = dx_pix.shape
    yy, xx = np.meshgrid(np.arange(H), np.arange(W), indexing="ij")
    yy = yy.astype(np.float64)
    xx = xx.astype(np.float64)
    out = np.full((n_dir, H, W), np.nan, np.float32)
    for idir in range(n_dir):
        best = np.full((H, W), -np.inf, np.float64)
        elev0 = None
        for k in range(ns):
            cy = yy + dy_pix[idir, k]
            cx = xx + dx_pix[idir, k]
            samp = map_coordinates(
                dem.astype(np.float64), [cy.ravel(), cx.ravel()],
                order=1, mode="constant", cval=np.nan,
            ).reshape(H, W)
            if k == 0:
                elev0 = samp
                continue
            ang = np.degrees(np.arctan2(samp - elev0, distances[k]))
            ang = np.where(np.isnan(samp), -np.inf, ang)
            best = np.fmax(best, ang)
        plane = np.where(np.isnan(elev0) | ~np.isfinite(best), np.nan, best)
        out[idir] = plane.astype(np.float32)
    return out


# --------------------------------------------------------------------------- #
# apply_ufunc wrappers (the required entry point for both backends)
# --------------------------------------------------------------------------- #
def make_ufunc(backend_func, az, dx_pix, dy_pix, distances):
    def f(dem_np):                                   # apply_ufunc hands us dem.values
        out = backend_func(dem_np, dx_pix, dy_pix, distances)   # (n_dir, H, W)
        return np.moveaxis(out, 0, -1)               # -> (y, x, azimuth)

    def call(dem_da):
        try:
            r = xr.apply_ufunc(
                f, dem_da,
                input_core_dims=[["y", "x"]],
                output_core_dims=[["y", "x", "azimuth"]],
                output_dtypes=[np.float32],
                dask="forbidden",
            )
            r = r.assign_coords(azimuth=("azimuth", az))
        except ValueError:
            # Fallback: some xarray versions reject reusing y,x as output core dims.
            r = xr.apply_ufunc(
                f, dem_da,
                input_core_dims=[["y", "x"]],
                output_core_dims=[["y_out", "x_out", "azimuth"]],
                output_dtypes=[np.float32],
                dask="forbidden",
            )
            r = r.rename({"y_out": "y", "x_out": "x"}).assign_coords(
                y=dem_da.y, x=dem_da.x, azimuth=("azimuth", az)
            )
        r = r.transpose("azimuth", "y", "x")
        r.attrs.update(units="degrees", long_name="terrain horizon angle")
        r.name = "horizon_angle"
        return r

    return call
