#!/usr/bin/env python3
#
# Copyright (C) 2025 Andy Aschwanden
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

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import matplotlib.pylab as plt
import numpy as np
import xarray as xr
from pathlib import Path
import re
from functools import partial
from collections import OrderedDict

def preprocess_config(
    ds,
    regexp: str = "id_(.+?)_",
    dim: str = "exp_id",
    drop_vars: list[str] | None = None,
    drop_dims: list[str] = ["nv4"],
) -> xr.Dataset:
    """
    Add experiment identifier to the dataset.

    This function processes the dataset by extracting an experiment identifier from the filename
    using a regular expression, adding it as a new dimension, and optionally dropping specified
    variables and dimensions from the dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset to be processed.
    regexp : str, optional
        The regular expression pattern to extract the experiment identifier from the filename, by default "id_(.+?)_".
    dim : str, optional
        The name of the new dimension to be added to the dataset, by default "exp_id".
    drop_vars : list[str]| None, optional
        A list of variable names to be dropped from the dataset, by default None.
    drop_dims : list[str], optional
        A list of dimension names to be dropped from the dataset, by default ["nv4"].

    Returns
    -------
    xarray.Dataset
        The processed dataset with the experiment identifier added as a new dimension, and specified variables and dimensions dropped.

    Raises
    ------
    AssertionError
        If the regular expression does not match any part of the filename.
    """
    m_id_re = re.search(regexp, ds.encoding["source"])
    ds = ds.expand_dims(dim)
    assert m_id_re is not None
    m_id: str | int
    try:
        m_id = int(m_id_re.group(1))
    except:
        m_id = str(m_id_re.group(1))
    ds[dim] = [m_id]

    p_config = ds["pism_config"]
    p_run_stats = ds["run_stats"]

    # List of suffixes to exclude
    suffixes_to_exclude = ["_doc", "_type", "_units", "_option", "_choices"]

    # Filter the dictionary
    config = {
        k: v
        for k, v in p_config.attrs.items()
        if not any(k.endswith(suffix) for suffix in suffixes_to_exclude)
    }
    stats = p_run_stats
    config_sorted = OrderedDict(sorted(config.items()))

    pc_keys = np.array(list(config_sorted.keys()))
    pc_vals = np.array(list(config_sorted.values()))
    rs_keys = np.array(list(stats.attrs.keys()))
    rs_vals = np.array(list(stats.attrs.values()))

    pism_config = xr.DataArray(
        pc_vals.reshape(-1, 1),
        dims=["pism_config_axis", "exp_id"],
        coords={"pism_config_axis": pc_keys, "exp_id": [m_id]},
        name="pism_config",
    )
    run_stats = xr.DataArray(
        rs_vals.reshape(-1, 1),
        dims=["run_stats_axis", "exp_id"],
        coords={"run_stats_axis": rs_keys, "exp_id": [m_id]},
        name="run_stats",
    )
    ds = xr.merge(
        [
            ds.drop_vars(["pism_config", "run_stats"], errors="ignore"),
            pism_config,
            run_stats,
        ]
    )
    return ds.drop_vars(drop_vars, errors="ignore").drop_dims(
        drop_dims, errors="ignore"
    )

if __name__ == "__main__":
    __spec__ = None  # type: ignore

    # set up the option parser
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "FILES",
        help="""Ensemble netCDF files.""",
        nargs="*",
    )

    options, unknown = parser.parse_known_args()
    m_files = options.FILES
    ds = xr.open_mfdataset(m_files, preprocess=partial(preprocess_config, regexp = "sg_(.+?)_"), decode_cf=True,
                           decode_timedelta=True, parallel=True).squeeze()

    m_var = "velsurf_mag"
    fig = ds[m_var].plot(col="exp_id", vmin=0, vmax=25, cmap="viridis").fig
    res = float(ds.sel({"pism_config_axis": "grid.dx"})["pism_config"].isel({"exp_id": -1}))
    fig.savefig(f"sg_stress_balance_{m_var}_{res}m.png")
    m_var = "velbase_mag"
    fig = ds[m_var].plot(col="exp_id", vmin=0, vmax=10, cmap="viridis").fig
    res = float(ds.sel({"pism_config_axis": "grid.dx"})["pism_config"].isel({"exp_id": -1}))
    fig.savefig(f"sg_stress_balance_{m_var}_{res}m.png")
