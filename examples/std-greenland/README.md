## "Standard" Greenland spinup examples

This directory contains scripts for standard Greenland ice dynamics simulations.
The first section of the PISM _User's Manual_ documents these examples.

These runs are motivated by (but not identical to) the "constant-climate" and
"paleo-climate" initial states which are evaluated in

  * A. Aschwanden, G. Adalgeirsdottir, and C. Khroulev (2013). _Hindcasting to_
    _measure ice sheet model sensitivity to initial states._ The Cryosphere 7,
    doi:10.5194/tc-7-1083-2013, 1083--1093.

The runs here use 5km SeaRISE data including present-day bedrock
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
generates additional climate-forcing files `pism_dT.nc` and `pism_dSL.nc`.

Then do

    $ ./spinup.sh

*In fact this does no actual simulation.*  Instead it provides a usage/help
message.  Specifically, it shows two basic usages, one with constant-climate
and SIA-only choices, and one with paleo-climate and SIA/SSA-hybrid choices.
The first of these examples can be run in a few minutes of computer time,
which should help you get started.

# SeaRISE Greenland background information

The SeaRISE assessment process included PISM Greenland simulations among many
ice sheet model results.  These runs are described in the following papers,
which evaluate and compare all SeaRISE simulations:

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

Note that the RACMO surface mass balance is missing where there is no
(present-day) ice.

For further references on the SeaRISE data, run

    $ ncdump -h Greenland_5km_v1.1.nc |less

after running `preprocess.sh` above, and look at the `source` attributes of
various fields.

Scripts for the actual PISM submitted runs to SeaRISE are more complex than
the runs shown here.  They can be found in the `stable0.5` branch of PISM.  See

    https://github.com/pism/pism/tree/stable0.5/examples/searise-greenland

# Basic Parameter study

An example for a basic parameter study is given by running the script `param.sh`
and either `runsequential.sh` or `runparallel.sh`.  See the User's Manual.  Also
see scripts in `advanced/`.

# Visualization

Results can be visualized using, e.g. `im-plot.py` from
[pypismtools](https://github.com/pism/pypismtools):

    im-plot.py -v velsurf_mag --colorbar_label --inner_titles "q=0.1 f=0.01,q=0.25 f=0.01,q=0.8 f=0.01,q=0.1 f=0.02,q=0.25 f=0.02,q=0.1 f=0.02,q=0.8 f=0.0,q=0.1 f=0.05,q=0.25 f=0.05,q=0.8 f=0.05" -o nosgl.pdf g20km_*0.01_*.nc g20km_*0.02_*.nc g20km_*0.05_*.nc

    im-plot.py -v velsurf_mag --colorbar_label --inner_titles "q=0.1 f=0.01,q=0.25 f=0.01,q=0.8 f=0.01,q=0.1 f=0.02,q=0.25 f=0.02,q=0.1 f=0.02,q=0.8 f=0.0,q=0.1 f=0.05,q=0.25 f=0.05,q=0.8 f=0.05" -o sgl.pdf g20km_*0.01.nc g20km_*0.02.nc g20km_*0.05.nc

