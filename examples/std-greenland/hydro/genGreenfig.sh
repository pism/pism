#!/bin/bash

# these figures show "inputs" to the hydrology model, i.e. fields
#   velsurf_mag, velbase_mag, bmelt
# and surfvelmag from Greenland_5km_v1.1.nc
# run as:
#   $ ./genGreenfig.sh g2km-init Greenland_5km_v1.1
# to generate figure files

RESULTROOT=$1
GREENROOT=$2

for varname in bmelt velsurf_mag velbase_mag; do
   ./basemapfigs.py $RESULTROOT $varname
done

./basemapfigs.py $GREENROOT surfvelmag

