#!/bin/bash

# depends on NCO (http://nco.sourceforge.net/)

set -e  # exit on error

DATA=nbreen_input.nc
PISMDATA=pismnbreen.nc

echo "making PISM-readable file $PISMDATA by copying $DATA and adjusting metadata ..."
rm -f $PISMDATA
# copy but skip "outline" which is (thk>0) and thus redundant
ncks -O -v outline -x $DATA $PISMDATA

# surface climate for "-surface given" just has placeholders
ncap2 -O -s "climatic_mass_balance=0.0*topg" $PISMDATA $PISMDATA
ncatted -O -a units,climatic_mass_balance,o,c,"m year-1" $PISMDATA
ncatted -O -a standard_name,climatic_mass_balance,o,c,"land_ice_surface_specific_mass_balance" $PISMDATA
ncap2 -O -s "ice_surface_temp=0.0*topg+260.0" $PISMDATA $PISMDATA
ncatted -O -a units,ice_surface_temp,o,c,"K" $PISMDATA
ncatted -O -a standard_name,ice_surface_temp,d,c, $PISMDATA

# edit mask name
ncatted -O -a standard_name,nbreen,d,c, $PISMDATA
ncatted -O -a long_name,nbreen,o,c,"mask for Nordenskioldbreen glacier on Svalbard" $PISMDATA

# input into subglacial aquifer is linear in surface elevation where nbreen==1, otherwise zero
ncap2 -O -s "bmelt=(3.0/900.0)*(900.0-usurf)" $PISMDATA $PISMDATA
ncap2 -O -s "where(usurf>900.0) bmelt=0.0" $PISMDATA $PISMDATA
ncap2 -O -s "where(nbreen<=0.5) bmelt=0.0" $PISMDATA $PISMDATA
ncatted -O -a units,bmelt,m,c,"m year-1" $PISMDATA
ncatted -O -a standard_name,bmelt,d,c, $PISMDATA
ncatted -O -a long_name,bmelt,o,c,"basal melt rate as input to subglacial hydrology" $PISMDATA

# 50 m/a basal speed; only the magnitude affects "-hydrology distributed"
ncap2 -O -s "u_ssa_bc=0.0*topg+50.0" $PISMDATA $PISMDATA
ncatted -O -a units,u_ssa_bc,o,c,"m year-1" $PISMDATA
ncatted -O -a long_name,u_ssa_bc,o,c,"x-component of prescribed sliding velocity" $PISMDATA
ncatted -O -a standard_name,u_ssa_bc,d,c, $PISMDATA
ncap2 -O -s "v_ssa_bc=0.0*topg" $PISMDATA $PISMDATA
ncatted -O -a units,v_ssa_bc,o,c,"m year-1" $PISMDATA
ncatted -O -a long_name,v_ssa_bc,o,c,"y-component of prescribed sliding velocity" $PISMDATA
ncatted -O -a standard_name,v_ssa_bc,d,c, $PISMDATA
ncap2 -O -s "bcflag=0.0*topg+1.0" $PISMDATA $PISMDATA
ncatted -O -a units,bcflag,d,c, $PISMDATA
ncatted -O -a long_name,bcflag,o,c,"equals one where (u_ssa_bc,v_ssa_bc) should be applied as sliding seen by hydrology" $PISMDATA
ncatted -O -a standard_name,bcflag,d,c, $PISMDATA

CONF=nbreen_config
echo "creating PISM-readable config override file $CONF.nc ..."
rm -f $CONF.nc
ncgen -o $CONF.nc $CONF.cdl

INTOBED=fakesummerevent.nc
echo "calling fake-inputtobed.py to create PISM-readable -input_to_bed file $INTOBED ..."
./fake-inputtobed.py

