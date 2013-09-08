# Greenland SeaRISE-like simulation example

The first section of the PISM _User's Manual_ documents this example.

This directory contains scripts for a Greenland ice dynamics simulation which
uses the SeaRISE data (see ) and which is similar to the PISM submitted runs
for the actual SeaRISE process.  The submitted runs are described, in the
context of evaluating and comparing all SeaRISE simulations, in

  * S. Nowicki and others (2013). _Insights into spatial sensitivities of ice_
    _mass response to environmental change from the SeaRISE ice sheet modeling_
    _project: II. Greenland._ J. Geophys. Res.: Earth Surface 118 (2),
    doi:10.1002/jgrf.20076, 1025--1044.

  * R. Bindshadler and others (2013). _Ice-sheet model sensitivities to_
    _environmental forcing and their use in projecting future sea-level_
    _(The SeaRISE Project)._ J. Glaciol. 59 (214), 195--224.

## Basic usage

First do

    $ ./preprocess.sh

This downloads the file `Greenland_5km_v1.1.nc` from the SeaRISE site and it
builds PISM-readable NetCDF file `pism_Greenland_5km_v1.1.nc`.  It also
generates additional climate-forcing files and it converts the configuration
overrides file `searise_config.cdl` to get `searise_config.nc`.

To run the spinup in the default (coarsest) grid resolution on 8 processes, for
example, and run it in the background and keep a copy of the stdout information
on the run in the file `out.spinup`, do

    $ ./spinup.sh 8 &> out.spinup &

It might be wise to see what will run as the spinup first; for that do

    $ PISM_DO=echo ./spinup.sh 8

There are also experiment and post-processing scripts which serve as examples,
but without maintained documentation.

## Changing configuration constants

Edit `searise_config.cdl` to set your preferred value and rerun `preprocess.sh`
and then `spinup.sh` or `experiments.sh` as desired.
