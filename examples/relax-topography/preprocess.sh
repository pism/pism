#!/bin/bash

# Copyright (C) 2009-2014, 2016, 2023 Ed Bueler and Andy Aschwanden

# downloads SeaRISE "Present Day Greenland" master dataset NetCDF file, adjusts metadata,
# saves under new name, ready for PISM

# depends on wget and NCO (ncrename, ncap, ncatted, ncpdq, ncks)

set -e  # exit on error

echo "# =================================================================================="
echo "# preprocessing SeaRISE-Greenland data for relaxation experiment"
echo "# =================================================================================="
echo

# get file; see page http://wiki2.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATAURL=https://github.com/pism/example-inputs/raw/main/std-greenland/
DATANAME=Greenland_5km_v$DATAVERSION.nc

echo "fetching master file ... "
wget -nc ${DATAURL}${DATANAME}
echo "  ... done."
echo

PISMVERSION=pism_$DATANAME
echo -n "creating bootstrapable $PISMVERSION from $DATANAME ..."
ncks -O $DATANAME $PISMVERSION  # just copies over, but preserves history and global attrs

# adjust metadata; uses NCO (http://nco.sourceforge.net/)
# convert from water equivalent rate to "kg m-2 year-1"; assumes water density 1000.0 kg m-3
ncap2 -O -s "precipitation=presprcp*1000.0" $PISMVERSION $PISMVERSION
ncatted -O -a units,airtemp2m,c,c,"degC" $PISMVERSION
ncatted -O -a units,precipitation,m,c,"kg m-2 year-1" $PISMVERSION
ncatted -O -a long_name,precipitation,c,c,"ice-equivalent mean annual precipitation rate" $PISMVERSION
# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -a standard_name,bheatflx,d,, $PISMVERSION
ncrename -O -v smb,climatic_mass_balance -v airtemp2m,ice_surface_temp $PISMVERSION $PISMVERSION
ncks -O -v lat,lon,bheatflx,topg,thk,precipitation,mapping,climatic_mass_balance,ice_surface_temp  $PISMVERSION $PISMVERSION
# convert from [m/year] ice equivalent to [kg m-2 year-1]; assumes ice density of 910.0 [kg m-3]
ncap2 -O -s "climatic_mass_balance=climatic_mass_balance*910.0" $PISMVERSION $PISMVERSION
ncatted -O -a units,climatic_mass_balance,m,c,"kg m-2 year-1" $PISMVERSION
echo "done."
echo

echo "create target ice thickness file"
ncap2 -O -s "thk=thk*0.; climatic_mass_balance[\$time,\$y1,\$x1]=-500.f*910.0" $PISMVERSION target_$PISMVERSION
ncatted -a units,climatic_mass_balance,o,c,"kg m-2 year-1" target_$PISMVERSION
echo "done."
echo

