#!/usr/bin/env python3
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

"""Download the Tarfala / Storglaciären radio-echo-sounding (RES) dataset and turn it
into a PISM-bootstrappable NetCDF file with a CF grid mapping.

The dataset (Wang et al., 2026; doi:10.17043/tarfala-storglaciaren-res-3, ODC-By v1.0)
provides 10 m gridded ice thickness and bed elevation GeoTIFFs on SWEREF99 TM
(EPSG:3006). This script downloads the archive, reads the bed and thickness rasters,
aligns them onto a common grid, and writes a NetCDF file with::

    bed        (standard_name bedrock_altitude,           m)
    thickness  (standard_name land_ice_thickness,         m)
    surface    (standard_name surface_altitude,           m) = bed + thickness

a `spatial_ref` grid-mapping variable carrying the EPSG:3006 CRS, and CF-style
projection coordinates `x`, `y`. PISM reads the projection from the grid-mapping
variable's `crs_wkt`.

The RES survey only covers the glacier itself, so outside it the bed (and hence the
surface) is filled with the public Copernicus GLO-30 (COP30) DEM, fetched from AWS Open
Data and resampled onto the 10 m grid (disable with `--no-fill`). Note the small vertical
datum mismatch between the RES product (geoid SWEN17_RH2000) and COP30 (EGM2008).

Examples
--------
    ./create_storglaciaren_res.py
    ./create_storglaciaren_res.py --bed raw -o sg_res_raw.nc
"""
from __future__ import annotations

import argparse
import sys
import time
import zipfile
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlopen

import numpy as np
import rioxarray  # noqa: F401  -- registers the `.rio` accessor
import xarray as xr
from pyproj import Transformer
from rasterio.enums import Resampling
from rioxarray.merge import merge_arrays

DATA_URL = (
    "https://bolin.su.se/data/s3/upload-next/next-s3-uploads/"
    "19d8b371-74cf-45dc-aa20-e3b038fd3bb6/tarfala-storglaciaren-res-3.zip"
)
ARCHIVE = "tarfala-storglaciaren-res-3.zip"
CRS = "EPSG:3006"  # SWEREF99 TM

BED_FILES = {"filtered": "SGBED_Filtered.tif", "raw": "SGBED.tif"}
THICKNESS_FILE = "SGIT.tif"

DEFAULT_OUTPUT = "storglaciaren_res.nc"

# Copernicus GLO-30 (COP30) global DEM, public on AWS Open Data (no authentication).
# Tiles are 1 deg x 1 deg, named by their SW corner; elevations are orthometric (EGM2008).
COP30_BUCKET = "https://copernicus-dem-30m.s3.amazonaws.com"


def download(url: str, dest: Path) -> None:
    """Stream *url* to *dest* unless *dest* already exists."""
    if dest.exists():
        print(f"  using cached archive {dest}")
        return
    print(f"  downloading {url}")
    with urlopen(url) as response, open(dest, "wb") as handle:  # noqa: S310 (trusted URL)
        chunk = response.read(1 << 20)
        while chunk:
            handle.write(chunk)
            chunk = response.read(1 << 20)
    print(f"  saved {dest} ({dest.stat().st_size / 1e6:.1f} MB)")


def extract(archive: Path, workdir: Path) -> None:
    """Extract *archive* into *workdir*."""
    with zipfile.ZipFile(archive) as zf:
        zf.extractall(workdir)


def find_file(workdir: Path, name: str) -> Path:
    """Return the (single) path to *name* anywhere under *workdir*."""
    matches = list(workdir.rglob(name))
    if not matches:
        raise FileNotFoundError(f"{name} not found in {workdir}")
    return matches[0]


def load_raster(path: Path) -> xr.DataArray:
    """Open a single-band GeoTIFF as a 2-D (y, x) DataArray with nodata masked to NaN."""
    da = rioxarray.open_rasterio(path, masked=True)
    if "band" in da.dims:
        da = da.squeeze("band", drop=True)
    return da


def cop30_tile_name(lat: int, lon: int) -> str:
    """COP30 tile name for the 1x1 degree tile whose SW corner is (lat, lon)."""
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return (f"Copernicus_DSM_COG_10_{ns}{abs(lat):02d}_00_"
            f"{ew}{abs(lon):03d}_00_DEM")


def cop30_tiles_for(reference: xr.DataArray):
    """Return the COP30 tile names covering the geographic extent of *reference*."""
    xmin, ymin, xmax, ymax = reference.rio.bounds()
    t = Transformer.from_crs(reference.rio.crs, "EPSG:4326", always_xy=True)
    lons, lats = t.transform([xmin, xmin, xmax, xmax], [ymin, ymax, ymin, ymax])
    tiles = []
    for lat in range(int(np.floor(min(lats))), int(np.floor(max(lats))) + 1):
        for lon in range(int(np.floor(min(lons))), int(np.floor(max(lons))) + 1):
            tiles.append(cop30_tile_name(lat, lon))
    return tiles


def fetch_cop30(reference: xr.DataArray, workdir: Path) -> xr.DataArray:
    """Download the COP30 tiles covering *reference* and reproject them onto its grid."""
    arrays = []
    for name in cop30_tiles_for(reference):
        dest = workdir / f"{name}.tif"
        url = f"{COP30_BUCKET}/{name}/{name}.tif"
        try:
            download(url, dest)
        except HTTPError as e:
            # ocean-only tiles do not exist; skip them
            print(f"  COP30 tile {name} unavailable (HTTP {e.code}), skipping")
            continue
        arrays.append(load_raster(dest))

    if not arrays:
        raise RuntimeError("no COP30 tiles could be downloaded for this domain")

    dem = arrays[0] if len(arrays) == 1 else merge_arrays(arrays)
    # reproject + bilinearly resample onto the reference (EPSG:3006, 10 m) grid
    return dem.rio.reproject_match(reference, resampling=Resampling.bilinear)


def build_dataset(bed: xr.DataArray, thickness: xr.DataArray,
                  fill_dem: xr.DataArray | None = None) -> xr.Dataset:
    """Assemble the PISM bootstrapping dataset from bed and thickness rasters."""
    # Put thickness on the bed grid if the two rasters do not already share it.
    if bed.rio.transform() != thickness.rio.transform() or bed.shape != thickness.shape:
        print("  regridding thickness onto the bed grid (reproject_match)")
        thickness = thickness.rio.reproject_match(bed)

    # No ice (or no measurement) -> zero thickness; clip the few negative grid artifacts.
    thk = thickness.fillna(0.0)
    thk = thk.where(thk > 0.0, 0.0)

    # Off the RES survey there is no ice, so the bedrock equals the bare-ground surface:
    # fill the missing bed with the external DEM (already on this grid). On-glacier bed and
    # thickness are kept, so surface = bed + thickness is the RES ice surface there and the
    # DEM elevation off-glacier.
    if fill_dem is not None:
        bed = bed.fillna(fill_dem)

    surface = bed + thk

    ds = xr.Dataset(
        {
            "bed": bed.astype("f4"),
            "thickness": thk.astype("f4"),
            "surface": surface.astype("f4"),
        }
    )

    ds["bed"].attrs.update(units="m", long_name="bedrock surface elevation",
                           standard_name="bedrock_altitude")
    ds["thickness"].attrs.update(units="m", long_name="land ice thickness",
                                 standard_name="land_ice_thickness")
    ds["surface"].attrs.update(units="m", long_name="ice upper surface elevation",
                               standard_name="surface_altitude")

    for axis, std in (("x", "projection_x_coordinate"),
                      ("y", "projection_y_coordinate")):
        ds[axis].attrs.update(units="m", standard_name=std,
                              long_name={"x": "easting", "y": "northing"}[axis])

    ds = ds.rio.write_crs(CRS).rio.set_spatial_dims(x_dim="x", y_dim="y")

    # Clean encoding: NaN-aware fill value for the floating-point fields, compression.
    for name in ("bed", "thickness", "surface"):
        ds[name].encoding.update(_FillValue=np.float32(-9999.0), zlib=True, complevel=4,
                                 dtype="float32")
    for name in ("x", "y"):
        ds[name].encoding.update(_FillValue=None, dtype="float64")

    ds.attrs.update(
        Conventions="CF-1.8",
        title="Storglaciaren radio-echo-sounding bed and ice thickness",
        source=DATA_URL,
        reference="Wang et al. (2026), doi:10.17043/tarfala-storglaciaren-res-3",
        license="ODC-By v1.0",
        history=time.asctime() + ": " + " ".join(sys.argv),
    )
    return ds


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-o", "--output", type=Path, default=Path(DEFAULT_OUTPUT),
                        help="output NetCDF file")
    parser.add_argument("--bed", choices=sorted(BED_FILES), default="filtered",
                        help="which bed raster to use (filtered = Gaussian-smoothed)")
    parser.add_argument("--url", default=DATA_URL, help="dataset archive URL")
    parser.add_argument("--workdir", type=Path, default=Path("res_data"),
                        help="directory for the downloaded/extracted files")
    parser.add_argument("--no-fill", dest="fill", action="store_false",
                        help="do not fill the off-survey bed/surface with the COP30 DEM")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    args.workdir.mkdir(parents=True, exist_ok=True)

    archive = args.workdir / ARCHIVE
    download(args.url, archive)
    extract(archive, args.workdir)

    bed_path = find_file(args.workdir, BED_FILES[args.bed])
    thk_path = find_file(args.workdir, THICKNESS_FILE)
    print(f"  bed:       {bed_path.name}")
    print(f"  thickness: {thk_path.name}")

    bed = load_raster(bed_path)
    thickness = load_raster(thk_path)

    fill_dem = None
    if args.fill:
        print("  filling off-survey area with the Copernicus GLO-30 (COP30) DEM")
        fill_dem = fetch_cop30(bed, args.workdir)

    ds = build_dataset(bed, thickness, fill_dem=fill_dem)

    ny, nx = ds["bed"].shape
    valid = float(np.isfinite(ds["bed"].values).mean())
    print(f"  grid: {nx} x {ny} @ {abs(ds.rio.resolution()[0]):.0f} m, CRS {CRS}")
    print(f"  bed coverage: {100 * valid:.1f}% of cells have a value"
          + ("" if args.fill else " (rest are off-survey nodata; use COP30 fill to complete)"))

    ds.to_netcdf(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
