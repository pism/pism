#!/bin/bash

# Copyright (C) 2011-2012, 2015, 2020, 2021 the PISM authors

# executes all preprocessing actions appropriate to Jakobshavn case,
# with a particular selection of terminus rectangle

# assumes local presence of "regional-tools" script pism_regional.py
# scripts called here download large files if they are not already present
# result is file jako.nc

set -x
set -e

./preprocess.sh

python3 pism_regional.py -i gr1km.nc -o jakomask.nc -x 360,382 -y 1135,1176 -b 50

ncks -O -d x,299,918 -d y,970,1394 gr1km.nc jako.nc
# rename x and y to avoid type clash (float in gr1km.nc and double in jako.nc)
ncrename -O -v x,x1 -v y,y1 jako.nc jako.nc
ncks -A -d x,299,918 -d y,970,1394 jakomask.nc jako.nc
ncap2 -A \
      -s 'land_ice_area_fraction_retreat=0 * thk' \
      -s 'where(thk > 0 || topg > 0) land_ice_area_fraction_retreat=1' \
      -s 'land_ice_area_fraction_retreat@standard_name=""' \
      -s 'land_ice_area_fraction_retreat@units="1"' \
      jako.nc jako.nc
