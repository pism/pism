#!/bin/bash

# Copyright (C) 2009-2014, 2016, 2017, 2018, 2019, 2020, 2023, 2024 The PISM Authors

# Downloads SeaRISE "Present Day Greenland" master dataset NetCDF file, adjusts
# metadata, and saves under new name ready for PISM.  See README.md.

# depends on wget and NCO (ncrename, ncap2, ncatted, ncpdq, ncks; see
# http://nco.sourceforge.net/)

set -e  # exit on error

echo "# =================================================================================="
echo "# PISM std Greenland example: preprocessing"
echo "# =================================================================================="
echo

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATAURL=https://github.com/pism/example-inputs/raw/main/std-greenland/
DATANAME=Greenland_5km_v$DATAVERSION.nc

echo "fetching master file ... "
wget -nc ${DATAURL}${DATANAME}   # -nc is "no clobber"
echo "  ... done."
echo

PISMVERSION=pism_$DATANAME
echo -n "creating bootstrapable $PISMVERSION from $DATANAME ... "
# copy the vars we want, and preserve history and global attrs
ncks -O -v mapping,lat,lon,bheatflx,topg,thk,presprcp,smb,airtemp2m $DATANAME $PISMVERSION
# convert from water equivalent thickness rate ("m year-1") to "kg m-2 year-1".
# Assumes water density of 1000.0 kg m-3
ncap2 -O -s "precipitation=presprcp*1000.0" $PISMVERSION $PISMVERSION
ncatted -O -a units,precipitation,m,c,"kg m-2 year-1" $PISMVERSION
ncatted -O -a long_name,precipitation,c,c,"mean annual precipitation rate" $PISMVERSION
# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -a standard_name,bheatflx,d,, $PISMVERSION
# use pism-recognized name for 2m air temp
ncrename -O -v airtemp2m,ice_surface_temp $PISMVERSION
ncatted -O -a units,ice_surface_temp,c,c,"degree_Celsius" $PISMVERSION
# use pism-recognized name and standard_name for surface mass balance, after
# converting from liquid water equivalent thickness per year to [kg m-2 year-1]
ncap2 -O -s "climatic_mass_balance=1000.0*smb" $PISMVERSION $PISMVERSION
# Note: The RACMO field smb has value 0 as a missing value, unfortunately,
# everywhere the ice thickness is zero.  Here we replace with 1 m a-1 ablation.
# This is a *choice* of the model of surface mass balance in thk==0 areas.
ncap2 -O -s "where(thk <= 0.0){climatic_mass_balance=-1000.0;}" $PISMVERSION $PISMVERSION
ncatted -O -a standard_name,climatic_mass_balance,m,c,"land_ice_surface_specific_mass_balance_flux" $PISMVERSION
ncatted -O -a units,climatic_mass_balance,m,c,"kg m-2 year-1" $PISMVERSION
# de-clutter by only keeping vars we want
ncks -O -v mapping,lat,lon,bheatflx,topg,thk,precipitation,ice_surface_temp,climatic_mass_balance \
  $PISMVERSION $PISMVERSION
# add projection information
ncatted -O -a proj,global,c,c,"+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" $PISMVERSION
ncap2 -A \
      -s 'land_ice_area_fraction_retreat=0 * thk' \
      -s 'where(thk > 0 || topg > 0) land_ice_area_fraction_retreat=1' \
      -s 'land_ice_area_fraction_retreat@units="1"' \
      $PISMVERSION $PISMVERSION
ncatted -a standard_name,land_ice_area_fraction_retreat,d,, $PISMVERSION

echo "done."
echo

# create a NCAP2 script that will add time bounds:
script=$(mktemp add_time_bounds_XXXX.txt)

cat > ${script} <<"EOF"
time@bounds="time_bnds";
defdim("nv",2);
time_bnds=array(0.0f,0.0f,/$time,$nv/);
time_bnds(:,1)=time;
time_bnds(1:,0)=time(:-2);
time_bnds(0,0)=2*time(0)-time(1);
EOF

# extract paleo-climate time series into files suitable for option
# -atmosphere ...,delta_T
TEMPSERIES=pism_dT.nc
echo -n "creating paleo-temperature file $TEMPSERIES from $DATANAME ... "
ncks -O -v oisotopestimes,temp_time_series $DATANAME $TEMPSERIES
ncrename -O -d oisotopestimes,time      $TEMPSERIES
ncrename -O -v temp_time_series,delta_T $TEMPSERIES
ncrename -O -v oisotopestimes,time      $TEMPSERIES
# reverse time dimension
ncpdq -O --rdr=-time $TEMPSERIES $TEMPSERIES
# make times follow same convention as PISM
ncap2 -O -s "time=-time" $TEMPSERIES $TEMPSERIES
ncatted -O -a units,time,m,c,"common_years since 1-1-1" $TEMPSERIES
ncatted -O -a calendar,time,c,c,"365_day" $TEMPSERIES
ncatted -O -a units,delta_T,m,c,"kelvin" $TEMPSERIES
ncap2 -O -S ${script} $TEMPSERIES  $TEMPSERIES
echo "done."
echo

# extract paleo-climate time series into files suitable for option
# -sea_level ...,delta_SL
SLSERIES=pism_dSL.nc
echo -n "creating paleo-sea-level file $SLSERIES from $DATANAME ... "
ncks -O -v sealeveltimes,sealevel_time_series $DATANAME $SLSERIES
ncrename -O -d sealeveltimes,time $SLSERIES
ncrename -O -v sealeveltimes,time $SLSERIES
ncrename -O -v sealevel_time_series,delta_SL $SLSERIES
# reverse time dimension
ncpdq -O --rdr=-time $SLSERIES $SLSERIES
# make times follow same convention as PISM
ncap2 -O -s "time=-time" $SLSERIES $SLSERIES
ncatted -O -a units,time,m,c,"common_years since 1-1-1" $SLSERIES
ncatted -O -a calendar,time,c,c,"365_day" $SLSERIES
ncap2 -O -S ${script} $SLSERIES  $SLSERIES
echo "done."
echo

rm ${script}
