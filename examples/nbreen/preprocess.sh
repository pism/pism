#!/bin/bash

# depends on NCO (ncrename, ncap, ncatted, ncpdq, ncks)

set -e  # exit on error

DATA=nbreen_input.nc
PISMDATA=pismnbreen.nc

echo "making PISM-readable file $PISMDATA by copying $DATA and adjusting metadata ..."
rm -f $PISMDATA
cp $DATA $PISMDATA

# following use NCO (http://nco.sourceforge.net/)
ncap2 -O -s "climatic_mass_balance=0.0*topg" $PISMDATA $PISMDATA
ncatted -O -a units,climatic_mass_balance,m,c,"m year-1" $PISMDATA
ncatted -O -a standard_name,climatic_mass_balance,m,c,"land_ice_surface_specific_mass_balance" $PISMDATA
ncap2 -O -s "ice_surface_temp=0.0*topg+260.0" $PISMDATA $PISMDATA
ncatted -O -a units,ice_surface_temp,m,c,"K" $PISMDATA
ncatted -O -a standard_name,ice_surface_temp,d,c, $PISMDATA
# Phi0(i,j) = (3.0/900.0)*(900.0-min(usurf(i,j),900.0))/p.spera;
ncap2 -O -s "bmelt=(3.0/900.0)*(900.0-usurf)" $PISMDATA $PISMDATA
ncatted -O -a units,bmelt,m,c,"m year-1" $PISMDATA
echo "  PISM-readable file $PISMDATA created."

