# Nbreen hydrology example

This directory contains draft scripts for a hydrology-only simulation using
ice geometry data for the Nordenskioldbreen glacier on Svalbard.  This example
is a basic test of mass-conservative subglacial hydrology models with cavity
evolution (i.e. cavitation by sliding and cavity closure by creep).  The basal
melt rate input is constant and the sliding is artificially set to 50 m/a
everywhere.

This example is related to a manuscript in preparation by Ward van Pelt and
Ed Bueler.

In the future there may be an additional ice dynamics simulation coupled to the
hydrology submodel.

## Basic usage

First do

    ./preprocess.sh

This builds PISM-readable NetCDF file `pismnbreen.nc` and it generates a
synthetic summer event input file `fakesummerevent.nc` (which has a large
pulse of water put into the subglacial system in a few weeks in the middle of
the year).  It also converts the `nbreen_config.cdl` to get `nbreen_config.nc`.

To see how to run examples, type

    ./run.sh

This will generate a usage message.


## Changing configuration constants

Edit `nbreen_config.cdl` to set your preferred value and rerun `preprocess.sh`.
