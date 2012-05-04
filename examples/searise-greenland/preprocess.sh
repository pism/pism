#!/bin/bash

# Copyright (C) 2009-2012 Ed Bueler and Andy Aschwanden

# PISM SeaRISE Greenland
#
# downloads SeaRISE "Present Day Greenland" master dataset NetCDF file, adjusts metadata,
# saves under new name, ready for PISM

# depends on wget and NCO (ncrename, ncap, ncatted, ncpdq, ncks)

set -e  # exit on error

echo "# =================================================================================="
echo "# PISM SeaRISE Greenland: preprocessing"
echo "# =================================================================================="
echo

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATAURL=http://websrv.cs.umt.edu/isis/images/a/a5/
DATANAME=Greenland_5km_v$DATAVERSION.nc

echo "fetching master file ... "
wget -nc ${DATAURL}${DATANAME}
echo "  ... done."
echo

PISMVERSION=pism_$DATANAME
echo -n "creating bootstrapable $PISMVERSION from $DATANAME ..."
ncks -O $DATANAME $PISMVERSION  # just copies over, but preserves history and global attrs

# adjust metadata; uses NCO (http://nco.sourceforge.net/)
# convert from water equiv to ice thickness change rate; assumes ice density 910.0 kg m-3
ncap -O -s "precipitation=presprcp*(1000.0/910.0)" $PISMVERSION $PISMVERSION
ncatted -O -a units,precipitation,a,c,"m a-1" $PISMVERSION
ncatted -O -a long_name,precipitation,a,c,"ice-equivalent mean annual precipitation rate" $PISMVERSION
# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -a standard_name,bheatflx,d,, $PISMVERSION
ncks -O -v lat,lon,bheatflx,topg,thk,precipitation,mapping \
  $PISMVERSION $PISMVERSION
echo "done."
echo

# extract time series into files suitable for -atmosphere ...,delta_T and -ocean ...,delta_SL
# compare resulting files to grip_dT.nc and specmap_dSL.nc in pism-dev/examples/eisgreen/
TEMPSERIES=pism_dT.nc
SLSERIES=pism_dSL.nc
echo -n "creating paleo-temperature file $TEMPSERIES from $DATANAME for option -atmosphere ...,delta_T ... "
ncks -O -v oisotopestimes,temp_time_series $DATANAME $TEMPSERIES
ncrename -O -d oisotopestimes,time -v oisotopestimes,time -v temp_time_series,delta_T $TEMPSERIES
ncpdq -O --rdr=-time $TEMPSERIES $TEMPSERIES  # reverse time dimension
ncap -O -s "time=-time" $TEMPSERIES $TEMPSERIES  # make times follow same convention as PISM
ncatted -O -a units,time,a,c,"years since 1-1-1" $TEMPSERIES
ncatted -O -a units,delta_T,m,c,"Kelvin" $TEMPSERIES
echo "done."
echo
echo -n "creating paleo-sea-level file $SLSERIES from $DATANAME for option -ocean ...,delta_SL ... "
ncks -O -v sealeveltimes,sealevel_time_series $DATANAME $SLSERIES
ncrename -O -d sealeveltimes,time -v sealeveltimes,time -v sealevel_time_series,delta_SL $SLSERIES
ncpdq -O --rdr=-time $SLSERIES $SLSERIES  # reverse time dimension
ncap -O -s "time=-time" $SLSERIES $SLSERIES  # make times follow same convention as PISM
ncatted -O -a units,time,a,c,"years since 1-1-1" $SLSERIES
echo "done."
echo

# generate config file
CONFIG=searise_config
echo -n "generating config file $CONFIG.nc from ascii file $CONFIG.cdl ... "
ncgen -o ${CONFIG}.nc ${CONFIG}.cdl
echo "done."
echo

echo "fetching anomaly files ..."
echo "   (these files can be re-generated from scripts in subdirectory future_forcing/) "
URL=http://www.pism-docs.org/download
PRECIP=ar4_precip_anomaly.nc
TEMP=ar4_temp_anomaly.nc
wget -nc ${URL}/$PRECIP
wget -nc ${URL}/$TEMP
echo "  ... done."
echo

ANOMALY=ar4_anomaly
echo -n "creating scaled anomaly files ... "
ncks -O $PRECIP ${ANOMALY}.nc
ncks -A $TEMP ${ANOMALY}.nc

ncks -O ${ANOMALY}.nc ${ANOMALY}_scalefactor_1.0.nc

ncap2 -O -s "precipitation_anomaly = 1.5 * precipitation_anomaly; air_temp_anomaly = 1.5 * air_temp_anomaly" ${ANOMALY}.nc ${ANOMALY}_scalefactor_1.5.nc
ncap2 -O -s "precipitation_anomaly = 2.0 * precipitation_anomaly; air_temp_anomaly = 2.0 * air_temp_anomaly" ${ANOMALY}.nc ${ANOMALY}_scalefactor_2.0.nc

echo "done."
echo

