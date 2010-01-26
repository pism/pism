#!/bin/bash

#   Downloads
#       grid20-EISMINT, suaq20-EISMINT, specmap.017, sum89-92-ss09-50yr.stp
# from
#       http://homepages.vub.ac.be/~phuybrec/eismint/
# using wget.  Then applies scripts eisgreen.py, eiscore.py, and fill_missing.
# See PISM User's Manual.

#set -e  # exit on error

for fname in "grid20-EISMINT" "suaq20-EISMINT" "specmap.017" "sum89-92-ss09-50yr.stp"; do
  echo "PREPROCESS.SH: getting $fname from http://homepages.vub.ac.be/~phuybrec/eismint/ ..."
  wget -nc http://homepages.vub.ac.be/~phuybrec/eismint/$fname
done

echo ""
echo "PREPROCESS.SH: running eisgreen.py to create eis_green20.nc ..."
./eisgreen.py

echo ""
echo "PREPROCESS.SH: running util/fill_missing.py to fill missing values in topg"
echo "  variable in data (eis_green20.nc) ..."
fill_missing.py -f eis_green20.nc -v topg -o eis_green_smoothed.nc
echo "  eis_green_smoothed.nc written ..."

echo ""
echo "PREPROCESS.SH: running eiscore.py to create grip_dT.nc and specmap_dSL.nc ..."
./eiscore.py

