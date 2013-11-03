## "Standard" Greenland spinup examples

This directory contains scripts for standard Greenland ice dynamics simulations.
The first section of the PISM _User's Manual_ documents these examples.  There
is one major run script `spinup.sh` but it takes options to determine major (and
minor) spinup choices.

These runs are motivated primarily by the "constant-climate" and "paleo-climate"
runs which are evaluated in

  * A. Aschwanden, G. Adalgeirsdottir, and C. Khroulev (2013). _Hindcasting to_
    _measure ice sheet model sensitivity to initial states._ The Cryosphere 7,
    doi:10.5194/tc-7-1083-2013, 1083--1093.

However, the runs here use 5km SeaRISE data including present-day bedrock
topography, ice thickness, and basal heat flux; see references at bottom.  Also,
in contrast to the above _Hindcasting_ paper, runs here use publicly-available
RACMO atmospheric/ocean circulation model outputs, including modeled present-day
2m air temperature, precipitation, and surface mass balance, as distributed by
the SeaRISE project; see reference at bottom.

## Basic usage

First do

    $ ./preprocess.sh

This downloads the file `Greenland_5km_v1.1.nc` from the SeaRISE site and it
builds PISM-readable NetCDF file `pism_Greenland_5km_v1.1.nc`.  It also
generates additional climate-forcing files `pism_dT.nc` and `pism_dSL.nc`
and it converts the configuration overrides file `std_config.cdl` to get
`std_config.nc`.

FIXME

To run the spinup in the default (coarsest) grid resolution on 8 processes, for
example, and run it in the background and keep a copy of the stdout information
on the run in the file `out.spinup`, do

    $ ./spinup.sh 8 &> out.spinup &

It might be wise to see what will run as the spinup first; for that do

    $ PISM_DO=echo ./spinup.sh 8

There are also experiment and post-processing scripts which serve as examples,
but without maintained documentation.

## Changing configuration constants

Edit `std_config.cdl` to set your preferred value and rerun `preprocess.sh`
and then `spinup.sh` as desired.

# SeaRISE Greenland background information

The SeaRISE assessment process included PISM Greenland simulations.  These runs
are described in these papers, which evaluate and compare all SeaRISE
simulations:

  * S. Nowicki and others (2013). _Insights into spatial sensitivities of ice_
    _mass response to environmental change from the SeaRISE ice sheet modeling_
    _project: II. Greenland._ J. Geophys. Res.: Earth Surface 118 (2),
    doi:10.1002/jgrf.20076, 1025--1044.

  * R. Bindshadler and others (2013). _Ice-sheet model sensitivities to_
    _environmental forcing and their use in projecting future sea-level_
    _(The SeaRISE Project)._ J. Glaciol. 59 (214), 195--224.

The reference for RACMO model outputs is:

  * J. Ettema, M.R. van den Broeke, E. van Meigaard, W.J. van de Berg,
    J.L. Bamber, J.E. Box, and R.C. Bales (2009). _Higher surface mass balance_
    _of the Greenland ice sheet revealed by high-resolution climate modeling._
    Geophys. Res. Lett. 36 (L12501), doi:10.1029/2009GL038110.

Scripts for the actual PISM submitted runs to SeaRISE can be found in the
`stable0.5` branch of PISM.  See

    https://github.com/pism/pism/tree/stable0.5/examples/searise-greenland

