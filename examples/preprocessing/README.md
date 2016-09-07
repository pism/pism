Ross flow model example
=================

This directory contains scripts demonstrating preparing datasets for use with PISM.

- `pism_python.py` uses Python to create a PISM-readable NetCDF file,
- `pism_matlab.m` uses MATLAB,
- `flowlineslab.py` creates an input file for a hypothetical flowline setup,

The Python module in `PISMNC.py` extends `netCDF4.Dataset` from
[`netcdf4-python`](http://unidata.github.io/netcdf4-python/) with
methods that make it easier to create NetCDF files PISM can use.
