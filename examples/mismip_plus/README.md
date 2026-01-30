## PICOP examples

This directory contains scripts for PICOP.

Here we use a MISMIP+ like geometry, but because PISM does not support axial-symmetric
coordinates, our geometry is mirrored at x=0. Due to a current limitation in PICO, only the
left side (x<0) gives the expected melt rates.

## Basic usage

First make sure you have the necessary python packages installed:

    $ pip install xarray h5netcdf cf_xarray pint_xarray

Then do

    $ ./spinup.py

This will give you an initial state.