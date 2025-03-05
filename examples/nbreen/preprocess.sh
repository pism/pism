#!/bin/bash

# depends on NCO (http://nco.sourceforge.net/)

set -e  # exit on error

DATA=nbreen_input.nc
PISMDATA=pismnbreen.nc

wget -nc https://github.com/pism/example-inputs/raw/main/nbreen/${DATA}

echo "making PISM-readable file $PISMDATA by copying $DATA and adjusting metadata ..."
rm -f $PISMDATA
# copy but skip "outline" which is (thk>0) and thus redundant
ncks -O -v outline -x $DATA $PISMDATA

# surface climate for "-surface given" just has placeholders
ncap2 -O -s "climatic_mass_balance=0.0*topg" $PISMDATA $PISMDATA
ncatted -O -a units,climatic_mass_balance,o,c,"kg m^-2 year^-1" $PISMDATA
ncatted -O -a standard_name,climatic_mass_balance,o,c,"land_ice_surface_specific_mass_balance_flux" $PISMDATA
ncap2 -O -s "ice_surface_temp=0.0*topg+260.0" $PISMDATA $PISMDATA
ncatted -O -a units,ice_surface_temp,o,c,"kelvin" $PISMDATA
ncatted -O -a standard_name,ice_surface_temp,d,c, $PISMDATA

# edit mask name
ncatted -O -a standard_name,nbreen,d,c, $PISMDATA
ncatted -O -a long_name,nbreen,o,c,"mask for Nordenskioldbreen glacier on Svalbard" $PISMDATA

# input into subglacial aquifer is linear in surface elevation where nbreen==1, otherwise zero
ncap2 -O -s "basal_melt_rate_grounded=(3.0/900.0)*(900.0-usurf)" $PISMDATA $PISMDATA
ncap2 -O -s "where(usurf>900.0) basal_melt_rate_grounded=0.0" $PISMDATA $PISMDATA
ncap2 -O -s "where(nbreen<=0.5) basal_melt_rate_grounded=0.0" $PISMDATA $PISMDATA
ncatted -O -a units,basal_melt_rate_grounded,m,c,"m year^-1" $PISMDATA
ncatted -O -a standard_name,basal_melt_rate_grounded,d,c, $PISMDATA
ncatted -O -a long_name,basal_melt_rate_grounded,o,c,"basal melt rate as input to subglacial hydrology" $PISMDATA

# 50 m/a basal speed; only the magnitude affects "-hydrology distributed"
ncap2 -O -s "u_bc=0.0*topg+50.0" $PISMDATA $PISMDATA
ncatted -O -a units,u_bc,o,c,"m year^-1" $PISMDATA
ncatted -O -a long_name,u_bc,o,c,"x-component of prescribed sliding velocity" $PISMDATA
ncatted -O -a standard_name,u_bc,d,c, $PISMDATA
ncap2 -O -s "v_bc=0.0*topg" $PISMDATA $PISMDATA
ncatted -O -a units,v_bc,o,c,"m year^-1" $PISMDATA
ncatted -O -a long_name,v_bc,o,c,"y-component of prescribed sliding velocity" $PISMDATA
ncatted -O -a standard_name,v_bc,d,c, $PISMDATA
ncap2 -O -s "vel_bc_mask=0.0*topg+1.0" $PISMDATA $PISMDATA
ncatted -O -a units,vel_bc_mask,d,c, $PISMDATA
ncatted -O -a long_name,vel_bc_mask,o,c,"equals one where (u_bc,v_bc) should be applied as sliding seen by hydrology" $PISMDATA
ncatted -O -a standard_name,vel_bc_mask,d,c, $PISMDATA

INTOBED=fakesummerevent.nc
echo "calling fake-inputtobed.py to create PISM-readable -input_to_bed file $INTOBED ..."
./fake-inputtobed.py

