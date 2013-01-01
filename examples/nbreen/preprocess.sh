#!/bin/bash

# depends on NCO (http://nco.sourceforge.net/)

set -e  # exit on error

DATA=nbreen_input.nc
PISMDATA=pismnbreen.nc

echo "making PISM-readable file $PISMDATA by copying $DATA and adjusting metadata ..."
rm -f $PISMDATA
cp $DATA $PISMDATA

ncap2 -O -s "climatic_mass_balance=0.0*topg" $PISMDATA $PISMDATA
ncatted -O -a units,climatic_mass_balance,m,c,"m year-1" $PISMDATA
ncatted -O -a standard_name,climatic_mass_balance,m,c,"land_ice_surface_specific_mass_balance" $PISMDATA
ncap2 -O -s "ice_surface_temp=0.0*topg+260.0" $PISMDATA $PISMDATA
ncatted -O -a units,ice_surface_temp,m,c,"K" $PISMDATA
ncatted -O -a standard_name,ice_surface_temp,d,c, $PISMDATA
# Phi0(i,j) = (3.0/900.0)*(900.0-min(usurf(i,j),900.0))/p.spera;
ncap2 -O -s "bmelt=(3.0/900.0)*(900.0-usurf)" $PISMDATA $PISMDATA
ncap2 -O -s "where(usurf>900.0) bmelt=0.0" $PISMDATA $PISMDATA
# outline(i,j)>0.5, Phi(i,j) = Phi0(i,j); end
ncap2 -O -s "where(nbreen<=0.5) bmelt=0.0" $PISMDATA $PISMDATA
ncatted -O -a units,bmelt,m,c,"m year-1" $PISMDATA

CONF=nbreen_config
echo "creating PISM-readable config override file $CONF.nc ..."
rm -f $CONF.nc
ncgen -o $CONF.nc $CONF.cdl

