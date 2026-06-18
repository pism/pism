#!/usr/bin/env python
"""Validate the PISM C++ terrain-horizon kernel against the solshade reference file.

Runs the real kernel (compiled into ./verify_horizon) on an interior crop of the example
DEM and compares the result to horizon(azimuth, y, x) in the reference NetCDF, trying a
set of candidate azimuth-index transforms to identify and confirm the convention.
"""
import subprocess
import sys

import numpy as np
import rioxarray  # noqa: F401
import xarray as xr

BOOT = "bootfile_RGI2000-v7.0-C-01-04374.nc"
REF = "insolation_RGI2000-v7.0-C-01-04374.nc"
N_DIR, MAX_DIST, STEP = 360, 2000.0, 20.0  # pism_compute_insolation defaults used for REF


def load_dem():
    ds = xr.open_dataset(BOOT).squeeze()
    mv = ds.rio.grid_mapping
    ds = ds.rio.write_crs(ds[mv].attrs["crs_wkt"])
    da = ds["surface"]
    # sort to y-ascending, x-ascending -> matches PISM (j north+, i east+)
    da = da.sortby("y").sortby("x")
    dx = abs(float(ds.rio.transform().a))
    dy = abs(float(ds.rio.transform().e))
    return da, dx, dy


def main():
    dem, dx, dy = load_dem()
    My, Mx = dem.shape  # (y, x)
    print(f"DEM {Mx}x{My}, dx={dx} dy={dy}, n_dir={N_DIR} max_dist={MAX_DIST} step={STEP}")

    dem.values.astype(np.float64).ravel(order="C").tofile("dem.bin")

    # interior crop, well away from edges so rays never leave the domain (20-cell reach)
    margin = int(MAX_DIST / dx) + 2
    i0, i1 = Mx // 2 - 50, Mx // 2 + 50
    j0, j1 = My // 2 - 50, My // 2 + 50
    assert i0 > margin and j0 > margin and i1 < Mx - margin and j1 < My - margin

    out = subprocess.run(
        ["./verify_horizon", "dem.bin", str(Mx), str(My), str(dx), str(dy),
         str(N_DIR), str(MAX_DIST), str(STEP), str(i0), str(i1), str(j0), str(j1)],
        capture_output=True, check=True,
    ).stdout
    mine = np.frombuffer(out, dtype=np.float64).reshape(N_DIR, j1 - j0, i1 - i0)
    mine = np.degrees(mine)  # kernel returns radians; file is in degrees

    # reference, sorted to the same y/x ordering, same crop
    ref = xr.open_dataset(REF)["horizon"].sortby("y").sortby("x")
    ref = ref.isel(y=slice(j0, j1), x=slice(i0, i1)).values  # (azimuth, j, i)

    print(f"crop cells: {(i1 - i0)}x{(j1 - j0)}; finite ref frac "
          f"{np.isfinite(ref).mean():.3f}")

    k = np.arange(N_DIR)
    candidates = {
        "identity            theta=k": k,
        "90 - k (CCW-from-E) ": (90 - k) % N_DIR,
        "90 + k              ": (90 + k) % N_DIR,
        "180 - k             ": (180 - k) % N_DIR,
        "270 - k             ": (270 - k) % N_DIR,
        "k + 180             ": (k + 180) % N_DIR,
        "-k                  ": (-k) % N_DIR,
    }

    print(f"\n{'azimuth transform':<24}{'mean|Δ| deg':>14}{'max|Δ| deg':>14}")
    best = None
    for name, idx in candidates.items():
        diff = np.abs(mine - ref[idx])
        finite = np.isfinite(diff)
        mad = float(np.mean(diff[finite]))
        mx = float(np.max(diff[finite]))
        print(f"{name:<24}{mad:>14.4f}{mx:>14.4f}")
        if best is None or mad < best[1]:
            best = (name, mad, mx)

    print(f"\nbest match: '{best[0].strip()}'  mean|Δ|={best[1]:.4f} deg, "
          f"max|Δ|={best[2]:.4f} deg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
