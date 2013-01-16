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

First build PISM-readable NetCDF file `pismnbreen.nc` and also convert the `.cdl`
to get `nbreen_config.nc` using the `preprocess.sh` script.  It also generates a
synthetic summer event input file `fakesummerevent.nc` which has a large
pulse of water put into the subglacial system in a few weeks in the middle of
the year.

    ./preprocess.sh

The following run uses the model in `src/base/hydrology/PISMDistributedHydrology.cc`:

    ./run.sh 4 500 5 dist &> out.y5_500m &

That is, it is a 4 processor run on 500m grid for 5 model years.  For the model
in `src/base/hydrology/PISMLakesHydrology.cc` use:

    ./run.sh 4 500 5 lakes &> out.y5_500m_lakes &

This is much faster, and it is also mass-conserving, but it lacks cavity
evolution.  Finally here is a high resolution run using the summer event input
file:

    ./run.sh 4 125 1 event &> out.y1_125m_event &


## Changing configuration constants

Edit `nbreen_config.cdl` to set your preferred value and rerun `preprocess.sh`.
