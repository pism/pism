#!/bin/bash

# downloads development version of SeaRISE "Present Day Antarctica" master
# dataset NetCDF file, adjusts metadata, breaks up, saves under new names,
# ready for PISM

# depends on wget and NCO (ncrename, ncap, ncatted, ncpdq, ncks)

set -e  # exit on error

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica
DATAVERSION=dev1.0
DATANAME=Antarctica_5km_$DATAVERSION.nc
echo
echo "downloading master $DATANAME ..."
wget -nc http://websrv.cs.umt.edu/isis/images/4/4d/$DATANAME

echo "making PISM-readable files by copying parts of $DATANAME"
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
ncap -O -s "temp_ma=temp+273.15" $PISMVERSION $PISMVERSION
ncatted -O -a units,temp_ma,a,c,"K" $PISMVERSION
# choose Van de Berg et al version of accumulation; will treat as ice-equivalent snow rate
ncrename -O -v accr,precip $PISMVERSION
ncatted -O -a units,precip,m,c,"m a-1" $PISMVERSION
# use bheatflx_shapiro as the default bheatflx data and 
ncrename -O -v bheatflx_shapiro,bheatflx $PISMVERSION
ncatted -O -a units,bheatflx,m,c,"W m-2" $PISMVERSION
# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -a standard_name,bheatflx,d,, $PISMVERSION
# keep only the fields we actually use at bootstrapping
ncks -O -v x,y,lat,lon,bheatflx,topg,thk,precip,temp_ma,mapping \
  $PISMVERSION $PISMVERSION
echo "  PISM-readable file $PISMVERSION created; only has fields"
echo "    used in bootstrapping."

# extract time series into files suitable for -dTforcing and -dSLforcing;
TEMPSERIES=pism_dT.nc
SLSERIES=pism_dSL.nc
ncks -O -v temptimes,temp_time_series $DATANAME $TEMPSERIES
ncrename -O -d temptimes,t -v temptimes,t -v temp_time_series,delta_T $TEMPSERIES
ncpdq -O --rdr=-t $TEMPSERIES $TEMPSERIES  # reverse time dimension so that
ncap -O -s "t=-t" $TEMPSERIES $TEMPSERIES  #   times follow same convention as PISM
ncatted -O -a units,t,a,c,"years since 1-1-1" $TEMPSERIES
echo "  PISM-readable paleo-temperature file $TEMPSERIES; for option -dTforcing"

ncks -O -v sealeveltimes,sealevel_time_series $DATANAME $SLSERIES
ncrename -O -d sealeveltimes,t -v sealeveltimes,t -v sealevel_time_series,delta_sea_level $SLSERIES
ncpdq -O --rdr=-t $SLSERIES $SLSERIES  # reverse time dimension so that
ncap -O -s "t=-t" $SLSERIES $SLSERIES  #   times follow same convention as PISM
ncatted -O -a units,t,a,c,"years since 1-1-1" $SLSERIES
echo "  PISM-readable paleo-sea-level file $SLSERIES; for option -dSLforcing"
echo

echo "now run spin-up script 'antspinup.sh'"
echo

