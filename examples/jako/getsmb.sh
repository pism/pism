#!/bin/bash

# Copyright (C) 2011-2012 the PISM authors

# downloads the SeaRISE 5km Greenland data set, from which we will only
# grab the "smb" = Surface Mass Balance and "airtemp2m" fields which were
# produced by the regional atmosphere model RACMO

# depends on wget and NCO (ncrename, ncap, ncwa, ncks)

set -e  # exit on error

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATAURL=http://websrv.cs.umt.edu/isis/images/a/a5/
DATANAME=Greenland_5km_v$DATAVERSION.nc
echo "fetching 5km SeaRISE data file which contains surface mass balance ... "
wget -nc ${DATAURL}${DATANAME}
echo "  ... done"
echo

CLIMATEFILE=g5km_climate.nc
echo "creating PISM-readable climate file $CLIMATEFILE from airtemp2m and smb in data file ..."
ncks -O -v mapping,smb,airtemp2m $DATANAME $CLIMATEFILE
ncrename -O -v airtemp2m,ice_surface_temp $CLIMATEFILE
ncatted -O -a units,ice_surface_temp,a,c,"Celsius" $CLIMATEFILE
ncap -O -s "climatic_mass_balance=(1000.0/910.0)*smb" $CLIMATEFILE $CLIMATEFILE
ncatted -O -a standard_name,climatic_mass_balance,a,c,"land_ice_surface_specific_mass_balance" $CLIMATEFILE
ncatted -O -a units,climatic_mass_balance,a,c,"meters/year" $CLIMATEFILE
ncks -O -x -v smb $CLIMATEFILE $CLIMATEFILE
echo "... done with creating climate file $CLIMATEFILE"
echo

