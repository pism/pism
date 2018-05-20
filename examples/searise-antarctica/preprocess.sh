#!/bin/bash

# Copyright (C) 2009-2018  PISM authors

# downloads development version of SeaRISE "Present Day Antarctica" master
# dataset NetCDF file, adjusts metadata, breaks up, saves under new names,
# ready for PISM

# depends on wget and NCO (ncrename, ncap2, ncatted, ncpdq, ncks)

set -e  # exit on error

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica
DATAVERSION=dev1.0
DATANAME=Antarctica_5km_$DATAVERSION.nc
echo
echo "downloading master $DATANAME ..."
wget -nc http://websrv.cs.umt.edu/isis/images/4/4d/$DATANAME

echo "making PISM-readable file by copying parts of $DATANAME"
echo "  and adjusting metadata ..."
PISMVERSION=pism_Antarctica_5km.nc
cp $DATANAME $PISMVERSION

# following use NCO (http://nco.sourceforge.net/)
# rename dimensions
ncrename -O -v x1,x -v y1,y -d x1,x -d y1,y $PISMVERSION
ncrename -O -v time,t -d time,t $PISMVERSION
# fix polar stereographic parameter
ncatted -O -a standard_parallel,mapping,m,d,-71.0 $PISMVERSION
# rename usurf for convenience
ncrename -O -v usrf,usurf $PISMVERSION
# fix surface temperature name and make K
ncap2 -O -s "air_temp=temp+273.15" $PISMVERSION $PISMVERSION
ncatted -O -a units,air_temp,m,c,"K" $PISMVERSION
# choose Van de Berg et al version of accumulation; will treat as
# ice-equivalent snow rate and convert from an ice-equivalent
# thickness rate ("m year-1") to "kg m-2 year-1" by assuming ice
# density of 910 kg m-3
ncap2 -O -s "precipitation=accr*910.0" $PISMVERSION $PISMVERSION
ncatted -O -a units,precipitation,m,c,"kg m-2 year-1" $PISMVERSION
# use bheatflx_shapiro as the default bheatflx data
ncrename -O -v bheatflx_shapiro,bheatflx $PISMVERSION
ncatted -O -a units,bheatflx,m,c,"W m-2" $PISMVERSION
# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -O -a standard_name,bheatflx,d,, $PISMVERSION
# keep only the fields we actually use at bootstrapping
ncks -O -v x,y,lat,lon,bheatflx,topg,thk,precipitation,air_temp,mapping \
     $PISMVERSION $PISMVERSION
# remove the time dimension
ncwa -O -a t $PISMVERSION $PISMVERSION
# remove the time (t) coordinate variable
ncks -O -x -v t $PISMVERSION $PISMVERSION
echo "  PISM-readable file $PISMVERSION created; only has fields"
echo "    used in bootstrapping."

echo "now run spin-up script 'antspin-coarse.sh'"
echo
