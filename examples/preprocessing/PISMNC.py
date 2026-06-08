#!/usr/bin/env python3
"""PISMNC: small helpers for building PISM-compatible NetCDF files with xarray.

This module replaces the old ``class PISMDataset(netCDF4.Dataset)`` shim. The
public surface is now a set of free functions that operate on an
``xarray.Dataset``. A typical preprocessing script looks like::

    import xarray as xr
    from PISMNC import (create_dimensions, append_time,
                        define_2d_field, write_2d_field, write_timeseries)

    ds = xr.Dataset()
    create_dimensions(ds, x, y, time_dependent=True)
    write_2d_field(ds, "topg", topg, attrs={"units": "m", ...})
    ds.to_netcdf(out_path, unlimited_dims=["time"])

The legacy class ``PISMDataset`` is kept as a thin compatibility shim that
delegates to these functions; new code should call the functions directly.
"""

import numpy as np

try:
    import xarray as xr
except ImportError as exc:
    raise SystemExit("xarray is not installed!") from exc


def create_time(ds, use_bounds=False, length=None, units=None):
    """Add a ``time`` coordinate (and optionally ``time_bounds``) to ``ds``.

    ``length`` controls the initial length of the time axis (``None`` makes
    it empty so it can be appended to via :func:`append_time`).
    """
    n = 0 if length is None else int(length)
    t_attrs = {
        "axis": "T",
        "long_name": "time",
        "units": units or "seconds since 1-1-1",
    }
    if use_bounds:
        t_attrs["bounds"] = "time_bounds"

    ds.coords["time"] = ("time", np.zeros(n, dtype="f8"), t_attrs)

    if use_bounds:
        ds["time_bounds"] = (("time", "n_bounds"),
                             np.zeros((n, 2), dtype="f8"))


def create_dimensions(ds, x, y, time_dependent=False, use_time_bounds=False):
    """Create PISM-compatible x/y (and optionally time) coordinates on ``ds``."""
    if time_dependent and "time" not in ds.variables:
        create_time(ds, use_bounds=use_time_bounds)

    ds.coords["x"] = ("x", np.asarray(x, dtype="f8"), {
        "axis": "X",
        "long_name": "X-coordinate in Cartesian system",
        "units": "m",
        "standard_name": "projection_x_coordinate",
    })
    ds.coords["y"] = ("y", np.asarray(y, dtype="f8"), {
        "axis": "Y",
        "long_name": "Y-coordinate in Cartesian system",
        "units": "m",
        "standard_name": "projection_y_coordinate",
    })


def append_time(ds, value, bounds=None):
    """Append one time step (and optional bounds pair) to ``ds``.

    Existing time-dependent variables are extended along ``time`` to match
    the new length; the new record is filled with zeros. This mirrors
    netCDF4's unlimited-time behaviour.

    xarray won't let us resize a coord while variables that share it have
    a different size, so we snapshot every time-dependent variable, drop
    them, resize ``time``, then re-add them padded to the new length.
    """
    if "time" not in ds.dims:
        return

    captured = {}
    for name in list(ds.variables):
        if name == "time":
            continue
        da = ds[name]
        if "time" in da.dims:
            captured[name] = (da.dims,
                              np.asarray(da.values).copy(),
                              dict(da.attrs),
                              dict(da.encoding))
            del ds[name]

    new_time = np.append(np.asarray(ds["time"].values), float(value))
    t_attrs = dict(ds["time"].attrs)
    ds.coords["time"] = ("time", new_time.astype("f8"), t_attrs)
    new_n = len(new_time)

    for name, (dims, arr, attrs, encoding) in captured.items():
        old_n = arr.shape[dims.index("time")] if "time" in dims else 0
        if name == "time_bounds":
            new_arr = np.zeros((new_n, 2), dtype="f8")
            if old_n > 0:
                new_arr[:old_n, :] = arr
            if bounds is not None:
                new_arr[new_n - 1, :] = np.asarray(bounds, dtype="f8")
        else:
            new_shape = tuple(new_n if d == "time" else ds.sizes[d]
                              for d in dims)
            new_arr = np.zeros(new_shape, dtype=arr.dtype)
            if old_n > 0:
                slc = tuple(slice(None, old_n) if d == "time" else slice(None)
                            for d in dims)
                new_arr[slc] = arr
        ds[name] = (dims, new_arr)
        ds[name].attrs.update(attrs)
        ds[name].encoding.update(encoding)


def set_attrs(ds, var_name, attrs):
    """Apply ``attrs`` (dict) to the named variable, skipping ``_FillValue``.

    ``_FillValue`` belongs to ``encoding`` in xarray; callers that need it
    should pass it through ``encoding={var_name: {"_FillValue": fv}}`` at
    write time, or set ``ds[var_name].encoding["_FillValue"] = fv`` directly.
    """
    if not attrs:
        return
    for name, value in attrs.items():
        if name == "_FillValue":
            ds[var_name].encoding["_FillValue"] = value
            continue
        ds[var_name].attrs[name] = value


def _default_dims(time_dependent):
    return ("time", "y", "x") if time_dependent else ("y", "x")


def define_2d_field(ds, var_name, time_dependent=False, dims=None,
                    nc_type="f8", attrs=None, data=None):
    """Define (and optionally initialise) a 2D-or-3D field on ``ds``.

    If ``data`` is ``None``, the variable is filled with zeros of the right
    shape; otherwise ``data`` is used directly (must broadcast to the dim
    sizes already present on ``ds``).
    """
    if dims is None:
        dims = _default_dims(time_dependent)

    if var_name in ds.variables:
        if attrs:
            set_attrs(ds, var_name, attrs)
        return ds[var_name]

    shape = tuple(ds.sizes[d] for d in dims)
    if data is None:
        data = np.zeros(shape, dtype=np.dtype(nc_type))
    else:
        data = np.asarray(data, dtype=np.dtype(nc_type))

    ds[var_name] = (dims, data)
    if attrs:
        set_attrs(ds, var_name, attrs)
    return ds[var_name]


def define_timeseries(ds, var_name, attrs=None):
    """Define a 1D time-series variable on ``ds``."""
    if var_name not in ds.variables:
        n = ds.sizes.get("time", 0)
        ds[var_name] = (("time",), np.zeros(n, dtype="f8"))
    if attrs:
        set_attrs(ds, var_name, attrs)
    return ds[var_name]


def write_2d_field(ds, var_name, data, time_dependent=False, attrs=None,
                   dims=None, nc_type="f8"):
    """Write a 2D map (optionally time-stamped) into ``ds``.

    For ``time_dependent=True``, the data is written into the *last* time
    record; this mirrors the original :class:`PISMDataset` behaviour.
    """
    if dims is None:
        dims = _default_dims(time_dependent)

    if var_name not in ds.variables:
        define_2d_field(ds, var_name, time_dependent=time_dependent,
                        dims=dims, nc_type=nc_type, attrs=attrs)
    elif attrs:
        set_attrs(ds, var_name, attrs)

    saved_attrs = dict(ds[var_name].attrs)
    saved_encoding = dict(ds[var_name].encoding)

    if time_dependent:
        last = ds.sizes["time"] - 1
        arr = np.asarray(ds[var_name].values).copy()
        arr[last, :, :] = data
        ds[var_name] = (dims, arr)
    else:
        ds[var_name] = (dims, np.asarray(data))

    ds[var_name].attrs.update(saved_attrs)
    ds[var_name].encoding.update(saved_encoding)
    return ds[var_name]


def write_timeseries(ds, var_name, data, attrs=None):
    """Write a (possibly multi-dim) time-series array into ``ds``.

    If the variable already exists, its dims are preserved so a 2D
    ``time_bounds(time, n_bounds)`` survives a write of a 2D array.
    """
    arr = np.asarray(data)
    if var_name in ds.variables:
        dims = ds[var_name].dims
        saved_attrs = dict(ds[var_name].attrs)
        saved_encoding = dict(ds[var_name].encoding)
    else:
        dims = ("time",) if arr.ndim == 1 else ("time", "n_bounds")
        saved_attrs = {}
        saved_encoding = {}

    ds[var_name] = (dims, arr.astype("f8") if arr.dtype.kind == "f" else arr)
    ds[var_name].attrs.update(saved_attrs)
    ds[var_name].encoding.update(saved_encoding)
    if attrs:
        set_attrs(ds, var_name, attrs)
    return ds[var_name]


def write(ds, var_name, data, time_dependent=False, attrs=None):
    """Dispatch to :func:`write_2d_field` or :func:`write_timeseries` based
    on the rank of ``data``."""
    arr = np.asarray(data)
    if arr.ndim == 1:
        return write_timeseries(ds, var_name, arr, attrs=attrs)
    if arr.ndim == 2:
        return write_2d_field(ds, var_name, arr, time_dependent, attrs=attrs)
    return None


# ---------------------------------------------------------------------------
# Backwards-compat shim. Existing scripts that did
#     from PISMNC import PISMDataset
#     nc = PISMDataset("foo.nc", "w")
#     nc.create_dimensions(...)
# can keep working unchanged; new scripts should use the functions above.
# ---------------------------------------------------------------------------
class PISMDataset:
    """Compatibility shim mimicking the netCDF4-backed PISMDataset class."""

    def __init__(self, filename, mode="w", **kwargs):
        self._filename = filename
        self._mode = mode
        if mode in ("r", "a") and filename:
            try:
                self.ds = xr.open_dataset(filename, decode_times=False,
                                          decode_cf=False).load()
                self.ds.close()
            except FileNotFoundError:
                self.ds = xr.Dataset()
        else:
            self.ds = xr.Dataset()

    # -- ergonomic shims (forwarders) ----------------------------------
    def create_time(self, *a, **k):       return create_time(self.ds, *a, **k)
    def create_dimensions(self, *a, **k): return create_dimensions(self.ds, *a, **k)
    def append_time(self, *a, **k):       return append_time(self.ds, *a, **k)
    def set_attrs(self, *a, **k):         return set_attrs(self.ds, *a, **k)
    def define_2d_field(self, *a, **k):   return define_2d_field(self.ds, *a, **k)
    def define_timeseries(self, *a, **k): return define_timeseries(self.ds, *a, **k)
    def write_2d_field(self, *a, **k):    return write_2d_field(self.ds, *a, **k)
    def write_timeseries(self, *a, **k):  return write_timeseries(self.ds, *a, **k)
    def write(self, *a, **k):             return write(self.ds, *a, **k)

    def sync(self):
        if self._mode in ("w", "a") and self._filename:
            self.ds.to_netcdf(
                self._filename, mode="w",
                unlimited_dims=["time"] if "time" in self.ds.dims else None,
            )

    def close(self):
        if self._mode in ("w", "a") and self._filename:
            self.ds.to_netcdf(
                self._filename, mode="w",
                unlimited_dims=["time"] if "time" in self.ds.dims else None,
            )

    # Allow `with PISMDataset(path, "w") as nc:` usage.
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    # Support the netCDF4-style ``nc.foo = value`` global-attribute assignment.
    def __setattr__(self, name, value):
        if name in ("_filename", "_mode", "ds"):
            object.__setattr__(self, name, value)
        else:
            self.ds.attrs[name] = value

    def __getattr__(self, name):
        # Forward netCDF4-style structural accessors to the underlying
        # xr.Dataset, plus reads of global attributes (so legacy code that
        # does e.g. nc.history works).
        if name == "variables":
            return self.__dict__["ds"].variables
        if name == "dimensions":
            return self.__dict__["ds"].dims
        if name == "ncattrs":
            return lambda: list(self.__dict__["ds"].attrs)
        if "ds" in self.__dict__ and name in self.__dict__["ds"].attrs:
            return self.__dict__["ds"].attrs[name]
        raise AttributeError(name)


if __name__ == "__main__":
    # Produce a NetCDF file for testing
    nx, ny = 101, 201
    x = np.linspace(-100, 100, nx)
    y = np.linspace(-100, 100, ny)
    xx, yy = np.meshgrid(x, y)

    ds = xr.Dataset()
    create_dimensions(ds, x, y, time_dependent=True, use_time_bounds=True)
    define_2d_field(ds, "xx", time_dependent=True,
                    attrs={"long_name": "xx", "comment": "test variable",
                           "valid_range": (-200.0, 200.0)})
    define_2d_field(ds, "yy", time_dependent=True)

    for t in [0, 1, 2, 3]:
        append_time(ds, t, (t - 1, t))
        write(ds, "xx", xx + t, time_dependent=True)
        write(ds, "yy", yy + 2 * t, time_dependent=True)

    ds.to_netcdf("foo.nc", unlimited_dims=["time"])
