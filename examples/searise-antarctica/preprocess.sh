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
ncatted -O -a standard_name,bheatflx,d,, $PISMVERSION
# keep only the fields we actually use at bootstrapping
ncks -O -v x,y,lat,lon,bheatflx,topg,thk,precip,temp_ma,mapping \
  $PISMVERSION $PISMVERSION
echo "  PISM-readable file $PISMVERSION created; only has fields"
echo "    used in bootstrapping."

# extract time series into files suitable for -atmosphere ...,delta_T and -ocean ...,delta_SL
TEMPSERIES=pism_dT.nc
SLSERIES=pism_dSL.nc
ncks -O -v temptimes,temp_time_series $DATANAME $TEMPSERIES
ncrename -O -d temptimes,t -v temptimes,t -v temp_time_series,delta_T $TEMPSERIES
ncpdq -O --rdr=-t $TEMPSERIES $TEMPSERIES  # reverse time dimension so that
ncap -O -s "t=-t" $TEMPSERIES $TEMPSERIES  #   times follow same convention as PISM
ncatted -O -a units,t,a,c,"years since 1-1-1" $TEMPSERIES
echo "  PISM-readable paleo-temperature file $TEMPSERIES; for option -atmosphere ...,delta_T"

ncks -O -v sealeveltimes,sealevel_time_series $DATANAME $SLSERIES
ncrename -O -d sealeveltimes,t -v sealeveltimes,t -v sealevel_time_series,delta_SL $SLSERIES
ncpdq -O --rdr=-t $SLSERIES $SLSERIES  # reverse time dimension so that
ncap -O -s "t=-t" $SLSERIES $SLSERIES  #   times follow same convention as PISM
ncatted -O -a units,t,a,c,"years since 1-1-1" $SLSERIES
echo "  PISM-readable paleo-sea-level file $SLSERIES; for option -ocean ...,delta_SL"
echo


# produce scaled future forcing data suitable to be used with options 
#-atmosphere pik,anomaly -surface simple -anomaly_temp ar4_ant_artm_anomaly_scalefactor_X.nc -anomaly_precip ar4_ant_precip_anomaly_scalefactor_X.nc
# from ANT_climate_forcing_2004_2098_v3.nc to be downloaded from http://websrv.cs.umt.edu/isis/index.php/Future_Climate_Data
# direct link: http://www.cs.umt.edu/files/ANT_climate_forcing_2004_2098_v3.nc
echo "downloading future forcing data ... "
#wget -nc http://www.cs.umt.edu/files/ANT_climate_forcing_2004_2098_v3.nc
wget -nc http://kluis.cs.umt.edu/isis/ANT_climate_forcing_2004_2098_v3.nc

echo "creating unscaled precip anomaly file ... "
ncks -O -v preciptation ANT_climate_forcing_2004_2098_v3.nc ar4_ant_precip_anomaly_scalefactor_1.0.nc
# change name and convert to ice-equivalent units;
# email 13 July 2011 from Bindshadler says Charles Jackson confirms density 1000.0 is correct
ncap2 -O -s 'precip=float(preciptation*(1000.0/910.0))' ar4_ant_precip_anomaly_scalefactor_1.0.nc ar4_ant_precip_anomaly_scalefactor_1.0.nc
# remove the unneeded 'preciptation' var
ncks -O -v preciptation -x ar4_ant_precip_anomaly_scalefactor_1.0.nc ar4_ant_precip_anomaly_scalefactor_1.0.nc
# remove now incorrect lwe_... standard name
ncatted -O -a standard_name,precip,d,, ar4_ant_precip_anomaly_scalefactor_1.0.nc

echo "creating unscaled temp anomaly file ... "
ncks -O -v annualtemp ANT_climate_forcing_2004_2098_v3.nc ar4_ant_artm_anomaly_scalefactor_1.0.nc
ncrename -O -v annualtemp,artm ar4_ant_artm_anomaly_scalefactor_1.0.nc

ncap2 -O -s 'precip(:,:,:)= (precip(:,:,:)-precip(0,:,:))' ar4_ant_precip_anomaly_scalefactor_1.0.nc ar4_ant_precip_anomaly_scalefactor_1.0.nc
ncap2 -O -s 'artm(:,:,:)= (artm(:,:,:)-artm(0,:,:))' ar4_ant_artm_anomaly_scalefactor_1.0.nc ar4_ant_artm_anomaly_scalefactor_1.0.nc

ncrename -O -v artm,temp_anomaly ar4_ant_artm_anomaly_scalefactor_1.0.nc ar4_ant_artm_anomaly_scalefactor_1.0.nc
ncrename -O -v precip,precip_anomaly ar4_ant_precip_anomaly_scalefactor_1.0.nc ar4_ant_precip_anomaly_scalefactor_1.0.nc

echo -n "creating scaled anomaly files ... factor 1.5 "
echo -n "precip_anomaly.. "
ncap2 -O -s 'precip_anomaly(:,:,:)= 1.5 * precip_anomaly(:,:,:)' ar4_ant_precip_anomaly_scalefactor_1.0.nc ar4_ant_precip_anomaly_scalefactor_1.5.nc
echo "temp_anomaly.. "
ncap2 -O -s 'temp_anomaly(:,:,:)= 1.5 * temp_anomaly(:,:,:)' ar4_ant_artm_anomaly_scalefactor_1.0.nc ar4_ant_artm_anomaly_scalefactor_1.5.nc

echo -n "creating scaled anomaly files ... factor 2.0 "
echo -n "precip_anomaly.. "
ncap2 -s 'precip_anomaly(:,:,:)= 2.0 * precip_anomaly(:,:,:)' ar4_ant_precip_anomaly_scalefactor_1.0.nc ar4_ant_precip_anomaly_scalefactor_2.0.nc
echo "temp_anomaly.. "
ncap2 -O -s 'temp_anomaly(:,:,:)= 2.0 * temp_anomaly(:,:,:)' ar4_ant_artm_anomaly_scalefactor_1.0.nc ar4_ant_artm_anomaly_scalefactor_2.0.nc

# For the temperature anomalies to be interpreted correctly they need to be in Kelvin:
ncatted -O -a units,temp_anomaly,o,c,"K" ar4_ant_artm_anomaly_scalefactor_1.0.nc ar4_ant_artm_anomaly_scalefactor_1.0.nc
ncatted -O -a units,temp_anomaly,o,c,"K" ar4_ant_artm_anomaly_scalefactor_1.5.nc ar4_ant_artm_anomaly_scalefactor_1.5.nc
ncatted -O -a units,temp_anomaly,o,c,"K" ar4_ant_artm_anomaly_scalefactor_2.0.nc ar4_ant_artm_anomaly_scalefactor_2.0.nc

echo "now run spin-up script 'antspin.sh'"
echo

