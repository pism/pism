#!/usr/bin/env python
"""Add CF-compliant 2D lat/lon coordinates (with cell-corner bounds) to a
NetCDF file whose grid is described by `x`, `y`, and a CF `grid_mapping`
variable.

CDO's `remapcon` / `remapycon` can only compute lat/lon cell corners for a
short list of projections (polar_stereographic, lambert_*, sinusoidal,
rotated_latitude_longitude). Files in other projections (e.g.
transverse_mercator) abort with "Source/Target grid cell corner
coordinates missing!". This script pre-computes the bounds with pyproj
and writes them as CF curvilinear-grid bounds, which CDO accepts for any
projection.

Usage:
    python add_latlon.py SRC.nc DST.nc

Behavior:
    - Reads x, y, and the projection (preferred: crs_wkt; fallback: CF
      grid_mapping attrs) from SRC.
    - Computes lat/lon at cell centers and 4-corner lat/lon bounds.
    - Writes DST containing every variable from SRC except any pre-
      existing lat/lon/*_bnds (those are dropped to keep the output
      unambiguous), plus fresh `lat`, `lon`, `lat_bnds`, `lon_bnds`.
    - Tags every (y, x)-shaped data variable with `coordinates = "lat lon"`.
    - Caches: if DST is newer than SRC and already has lat_bnds + lon_bnds,
      does nothing.
"""
import os
import sys

import netCDF4 as nc
import numpy as np
import pyproj


_LATLON_SKIP = {
    "lat", "lon", "latitude", "longitude",
    "lat_bnds", "lon_bnds", "latitude_bnds", "longitude_bnds",
}


def ensure_obs_latlon(in_file: str, out_file: str) -> None:
    if (os.path.exists(out_file)
            and os.path.getmtime(out_file) >= os.path.getmtime(in_file)):
        try:
            with nc.Dataset(out_file) as cached:
                if {"lat", "lon", "lat_bnds", "lon_bnds"} <= set(cached.variables):
                    return
        except OSError:
            pass

    with nc.Dataset(in_file) as src:
        x = src.variables["x"][:]
        y = src.variables["y"][:]
        nx, ny = x.size, y.size
        dx = float(x[1] - x[0])
        dy = float(y[1] - y[0])

        crs_var = None
        for v in src.variables.values():
            if "grid_mapping_name" in v.ncattrs():
                crs_var = v
                break
        if crs_var is None:
            raise RuntimeError(f"No grid_mapping variable found in {in_file}")

        if "crs_wkt" in crs_var.ncattrs():
            crs = pyproj.CRS.from_wkt(crs_var.getncattr("crs_wkt"))
        else:
            crs = pyproj.CRS.from_cf({
                k: crs_var.getncattr(k) for k in crs_var.ncattrs()
            })
        transformer = pyproj.Transformer.from_crs(crs, "EPSG:4326",
                                                  always_xy=True)

        xx, yy = np.meshgrid(x, y)
        lon_c, lat_c = transformer.transform(xx, yy)

        x_corners = np.concatenate(([x[0] - dx / 2],
                                    x[:-1] + dx / 2,
                                    [x[-1] + dx / 2]))
        y_corners = np.concatenate(([y[0] - dy / 2],
                                    y[:-1] + dy / 2,
                                    [y[-1] + dy / 2]))
        xc, yc = np.meshgrid(x_corners, y_corners)
        lon_corner, lat_corner = transformer.transform(xc, yc)

        lat_bnds = np.stack([lat_corner[:-1, :-1], lat_corner[:-1, 1:],
                             lat_corner[1:, 1:], lat_corner[1:, :-1]],
                            axis=-1)
        lon_bnds = np.stack([lon_corner[:-1, :-1], lon_corner[:-1, 1:],
                             lon_corner[1:, 1:], lon_corner[1:, :-1]],
                            axis=-1)

        with nc.Dataset(out_file, "w") as dst:
            for attr in src.ncattrs():
                dst.setncattr(attr, src.getncattr(attr))
            for dname, dim in src.dimensions.items():
                dst.createDimension(dname,
                                    len(dim) if not dim.isunlimited() else None)
            if "nv" not in dst.dimensions:
                dst.createDimension("nv", 4)

            for vname, v in src.variables.items():
                if vname in _LATLON_SKIP:
                    continue
                new_v = dst.createVariable(
                    vname, v.dtype, v.dimensions,
                    zlib=True, complevel=4,
                    fill_value=v.getncattr("_FillValue")
                    if "_FillValue" in v.ncattrs() else None,
                )
                for attr in v.ncattrs():
                    if attr == "_FillValue":
                        continue
                    new_v.setncattr(attr, v.getncattr(attr))
                new_v[:] = v[:]
                if (v.dimensions[-2:] == ("y", "x")
                        and vname not in ("x", "y")):
                    new_v.coordinates = "lat lon"

            lat_v = dst.createVariable("lat", "f8", ("y", "x"),
                                       zlib=True, complevel=4)
            lat_v.units = "degrees_north"
            lat_v.standard_name = "latitude"
            lat_v.long_name = "latitude"
            lat_v.bounds = "lat_bnds"
            lat_v[:] = lat_c

            lon_v = dst.createVariable("lon", "f8", ("y", "x"),
                                       zlib=True, complevel=4)
            lon_v.units = "degrees_east"
            lon_v.standard_name = "longitude"
            lon_v.long_name = "longitude"
            lon_v.bounds = "lon_bnds"
            lon_v[:] = lon_c

            lat_b = dst.createVariable("lat_bnds", "f8", ("y", "x", "nv"),
                                       zlib=True, complevel=4)
            lat_b[:] = lat_bnds
            lon_b = dst.createVariable("lon_bnds", "f8", ("y", "x", "nv"),
                                       zlib=True, complevel=4)
            lon_b[:] = lon_bnds


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} SRC.nc DST.nc", file=sys.stderr)
        sys.exit(2)
    ensure_obs_latlon(sys.argv[1], sys.argv[2])
