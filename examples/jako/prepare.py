# Copyright (C) 2026 the PISM authors

from __future__ import annotations

import hashlib
import logging
import os
from pathlib import Path
from typing import Iterable
from urllib.parse import urlparse

import boto3
from boto3.s3.transfer import TransferConfig
from botocore.config import Config
from tqdm import tqdm

import geopandas as gpd
import numpy as np
import rioxarray  # noqa: F401
import xarray as xr

def create_domain(
    x_bnds: list | np.ndarray,
    y_bnds: list | np.ndarray,
    resolution: float | None = None,
    x_dim: str = "x",
    y_dim: str = "y",
    crs: str = "EPSG:3413",
) -> xr.Dataset:
    """
    Create an xarray.Dataset representing a domain with specified x and y boundaries.

    Parameters
    ----------
    x_bnds : list or numpy.ndarray
        A list or array containing the minimum and maximum x-coordinate boundaries.
    y_bnds : list or numpy.ndarray
        A list or array containing the minimum and maximum y-coordinate boundaries.
    resolution : float or None, optional
        The resolution of the grid, by default None.
    x_dim : str, optional
        The name of the x dimension, by default "x".
    y_dim : str, optional
        The name of the y dimension, by default "y".
    crs : str, optional
        The coordinate reference system (CRS) for the domain, by default "EPSG:3413".

    Returns
    -------
    xarray.Dataset
        An xarray.Dataset containing the domain information, including coordinates,
        boundary data, and mapping attributes.

    Notes
    -----
    The dataset includes:
    - `x` and `y` coordinates with associated metadata.
    - A `mapping` DataArray with polar stereographic projection attributes.
    - A `domain` DataArray with a reference to the `mapping`.
    - `x_bnds` and `y_bnds` DataArrays representing the boundaries of the domain.

    Examples
    --------
    >>> x_bnds = [0, 1000]
    >>> y_bnds = [0, 2000]
    >>> ds = create_domain(x_bnds, y_bnds)
    >>> print(ds)
    """

    if resolution is not None:
        # np.arange(start, stop) includes start but not stop
        xb = np.arange(x_bnds[0], x_bnds[1] + resolution, resolution)
        yb = np.arange(y_bnds[0], y_bnds[1] + resolution, resolution)
        x = np.arange(x_bnds[0] + resolution / 2, x_bnds[1] + resolution - resolution / 2, resolution)
        y = np.arange(y_bnds[0] + resolution / 2, y_bnds[1] + resolution - resolution / 2, resolution)
        x_bounds = np.stack([xb[:-1], xb[1:]]).T
        y_bounds = np.stack([yb[:-1], yb[1:]]).T
    else:
        x = np.array([0])
        y = np.array([0])
        x_bounds = np.array([[x_bnds[0], x_bnds[1]]])
        y_bounds = np.array([[y_bnds[0], y_bnds[1]]])

    x_bnds_dim = f"{x_dim}_bnds"
    y_bnds_dim = f"{y_dim}_bnds"
    coords = {
        x_dim: (
            [x_dim],
            x,
            {
                "axis": "X",
                "bounds": x_bnds_dim,
                "units": "m",
                "standard_name": "projection_x_coordinate",
                "long_name": "X-coordinate in projected coordinate system",
            },
        ),
        y_dim: (
            [y_dim],
            y,
            {
                "axis": "Y",
                "bounds": y_bnds_dim,
                "units": "m",
                "standard_name": "projection_y_coordinate",
                "long_name": "Y-coordinate in projected coordinate system",
            },
        ),
    }
    ds = xr.Dataset(
        {
            "domain": xr.DataArray(
                data=0,
                dims=[y_dim, x_dim],
                coords={x_dim: coords[x_dim], y_dim: coords[y_dim]},
                attrs={
                    "dimensions": f"{x_dim} {y_dim}",
                },
            ),
            x_bnds_dim: xr.DataArray(
                data=x_bounds,
                dims=[x_dim, "nv2"],
                coords={x_dim: coords[x_dim]},
            ),
            y_bnds_dim: xr.DataArray(
                data=y_bounds,
                dims=[y_dim, "nv2"],
                coords={y_dim: coords[y_dim]},
            ),
        },
        attrs={"Conventions": "CF-1.8"},
    ).rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim)
    ds.rio.write_crs(crs, inplace=True).rio.write_coordinate_system(inplace=True)
    for var in list(ds.data_vars) + list(ds.coords):
        ds[var].encoding["_FillValue"] = None

    return ds

def download_from_s3(s3_uri: str, dest: str | Path) -> Path:
    """
    Download a file from AWS S3.

    Parameters
    ----------
    s3_uri : str
        URI of S3 object to download.
    dest : str or Path
        Path to the downloaded file.

    Returns
    -------
    Path
        Path to the downloaded file.
    """
    dest = Path(dest)

    parsed_url = urlparse(s3_uri)
    bucket = parsed_url.netloc
    prefix = parsed_url.path.lstrip("/")

    s3 = boto3.client("s3")
    head = s3.head_object(Bucket=bucket, Key=prefix)
    total_size = head["ContentLength"]

    with tqdm(total=total_size, unit="B", unit_scale=True, desc=dest.name) as pbar:
        s3.download_file(bucket, prefix, str(dest), Callback=pbar.update)

    return dest



x_min = -199925 + 4500
x_max = -50025 - 9000
y_min = -2314025 + 4500
y_max = -2224025 - 4500

bucket = "pism-cloud-data"

obs_file = "GreenlandObsISMIP7-v1.3.nc"

if not Path(obs_file).exists():
    obs = download_from_s3(f"s3://{bucket}/obs/{obs_file}", ".")
else:
    obs  = xr.open_dataset(obs_file)
crs = obs.mapping.crs_wkt
obs.rio.write_crs(crs, inplace=True).rio.write_grid_mapping("mapping", inplace=True).rio.write_coordinate_system(inplace=True)

outline_file = "Greenland_Basins_PS_v1.4.2_all_with_shelves.gpkg"
outline_local = Path(".") / outline_file
outline_s3 =  f"""s3://{bucket}/basins/{outline_file}"""
if not outline_local.exists():
    download_from_s3(outline_s3, outline_local)
    
prefix = "kitp/input"
regrid_file = "g900m_id_BAYES-MEDIAN_2007-01-01_2010-01-01.nc"
version = "v4"

if not Path(regrid_file).exists():
    download_from_s3(f"s3://{bucket}/{prefix}/{version}/{regrid_file}", regrid_file)

regrid = xr.open_dataset(regrid_file)
regrid = regrid.sel({"x": slice(x_min - 1500, x_max + 1500), "y": slice(y_min - 1500, y_max + 1500)})

x = regrid["x"].values
y = regrid["y"].values

# True for cells outside the active model domain (i.e., in the 1500m
# buffer we tacked onto each side via .sel above).
in_x_band = (x < x_min) | (x > x_max)        # shape (nx,)
in_y_band = (y < y_min) | (y > y_max)        # shape (ny,)

# 2D mask, (y, x) order to match PISM conventions: 1 in the buffer,
# 0 in the active interior.
no_model = np.zeros((len(y), len(x)), dtype=np.int8)
no_model[in_y_band, :] = 1
no_model[:, in_x_band] = 1

regrid["no_model_mask"] = (("y", "x"), no_model, {
    "long_name": "regional no-model boundary mask "
    "(1 = held fixed at boundary state, 0 = active interior)",
    "units": "1",
})
regrid.to_netcdf("state_jako.nc")
    
basins = gpd.read_file(outline_local).to_crs(crs)
jib = basins[basins["NAME"] == "JAKOBSHAVN_ISBRAE"]


if obs.y[0] > obs.y[-1]:
    obs = obs.sortby("y")
obs_jib = obs.sel({"x": slice(x_min - 1500, x_max + 1500), "y": slice(y_min - 1500, y_max + 1500)})

# coordinate metadata on x/y
for c, axis, stdname in (("x", "X", "projection_x_coordinate"), ("y", "Y", "projection_y_coordinate")):
    attrs = {
        "units": "m",
        "axis": axis,
        "standard_name": stdname,
        "long_name": f"{c}-coordinate in projected coordinate system",
    }
    if c in obs_jib.coords:
        obs_jib[c].attrs.update(attrs)

zeta_mask = obs_jib["bed"].rio.clip(jib.geometry, drop=False)
outside_geometry = zeta_mask.isnull()
no_ice = (obs_jib["thickness"] == 0) | obs_jib["bed"].isnull()
obs_jib["zeta_fixed_mask"] = (
    xr.where(outside_geometry | no_ice, 1, 0).fillna(0).astype(int)
)
  
def fix_xy_attrs(ds):
    """Ensure x/y coordinates have proper CF attributes and float64 dtype."""
    for c, axis, stdname in (("x", "X", "projection_x_coordinate"),
                              ("y", "Y", "projection_y_coordinate")):
        if c in ds.coords:
            ds[c].attrs.update({
                "units": "m",
                "axis": axis,
                "standard_name": stdname,
                "long_name": f"{c}-coordinate in projected coordinate system",
            })
            ds[c].encoding["dtype"] = "float64"
    return ds

obs_jib = obs_jib.fillna(0).isel(vel_time=1)
boot_ds = fix_xy_attrs(obs_jib[["bed", "thickness", "surface_grimp"]])
boot_ds["bed"].attrs.update({"standard_name": "bedrock_altitude", "units": "m"})
boot_ds["surface_grimp"].attrs.update({"standard_name": "surface_altitude", "units": "m"})
encoding = {v: {"_FillValue": None} for v in (boot_ds.data_vars and boot_ds.coords)}
boot_ds.to_netcdf("boot_jako.nc", encoding=encoding)

obs_ds = fix_xy_attrs(
    obs_jib[["vx_timeseries", "vy_timeseries", "zeta_fixed_mask"]]
    .rename_vars({"vx_timeseries": "u_observed", "vy_timeseries": "v_observed"})
    .fillna(0)
)
encoding = {v: {"_FillValue": None} for v in (obs_ds.data_vars and obs_ds.coords)}

obs_ds.to_netcdf("obs_jako.nc", encoding=encoding)

grid_file = Path("grid_jako.nc")
grid_ds = create_domain([x_min + 75, x_max -75], [y_min + 75, y_max-75])
grid_ds.to_netcdf(grid_file)
