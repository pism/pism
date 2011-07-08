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


# produce scaled future forcing data suitable to be used with options -surface given -bc_file foo.nc
# from ANT_climate_forcing_2004_2098_v3.nc to be downloaded from http://websrv.cs.umt.edu/isis/index.php/Future_Climate_Data
# direct link: http://www.cs.umt.edu/files/ANT_climate_forcing_2004_2098_v3.nc
echo -n "downloading future forcing data ... "
wget -nc http://www.cs.umt.edu/files/ANT_climate_forcing_2004_2098_v3.nc
echo "creating scaled files ... "
ncks -O ANT_climate_forcing_2004_2098_v3.nc ar4_ant_scalefactor_1.0.nc
ncrename -O -v preciptation,acab ar4_ant_scalefactor_1.0.nc
ncrename -O -v annualtemp,artm ar4_ant_scalefactor_1.0.nc

echo -n "creating scaled files ... times 1.5 "
ncap2 -s 'artm(:,:,:)= artm(0,:,:) + 1.5 * (artm(:,:,:)-artm(0,:,:))' -s 'acab(:,:,:)= acab(0,:,:) + 1.5 * (acab(:,:,:)-acab(0,:,:))' ar4_ant_scalefactor_1.0.nc ar4_ant_scalefactor_1.5.nc

echo "creating scaled files ... times 2.0 "
ncap2 -s 'artm(:,:,:)= artm(0,:,:) + 2.0 * (artm(:,:,:)-artm(0,:,:))' -s 'acab(:,:,:)= acab(0,:,:) + 2.0 * (acab(:,:,:)-acab(0,:,:))' ar4_ant_scalefactor_1.0.nc ar4_ant_scalefactor_2.0.nc


# as an alternative, create anomaly files, suitable (TODO check whether this also works with -surface pik) to be used
# with options -surface given,anomaly --anomaly_artm ar4_ant_artm_anomaly_scalefactor_X.nc -anomaly_acab ar4_ant_acab_anomaly_scalefactor_X.nc

echo "creating anomaly files ... "
ncks -v acab ar4_ant_scalefactor_1.0.nc ar4_ant_acab_anomaly_scalefactor_1.0.nc
ncks -v artm ar4_ant_scalefactor_1.0.nc ar4_ant_artm_anomaly_scalefactor_1.0.nc

ncap2 -O -s 'acab(:,:,:)= (acab(:,:,:)-acab(0,:,:))' ar4_ant_acab_anomaly_scalefactor_1.0.nc ar4_ant_acab_anomaly_scalefactor_1.0.nc
ncap2 -O -s 'artm(:,:,:)= (artm(:,:,:)-artm(0,:,:))' ar4_ant_artm_anomaly_scalefactor_1.0.nc ar4_ant_artm_anomaly_scalefactor_1.0.nc

echo -n "creating scaled files ... times 1.5 "
echo -n "acab.. "
ncap2 -s 'acab(:,:,:)= 1.5 * acab(:,:,:)' ar4_ant_acab_anomaly_scalefactor_1.0.nc ar4_ant_acab_anomaly_scalefactor_1.5.nc
echo "artm.. "
ncap2 -s 'artm(:,:,:)= 1.5 * artm(:,:,:)' ar4_ant_artm_anomaly_scalefactor_1.0.nc ar4_ant_artm_anomaly_scalefactor_1.5.nc

echo -n "creating scaled files ... times 2.0 "
echo -n "acab.. "
ncap2 -s 'acab(:,:,:)= 2.0 * acab(:,:,:)' ar4_ant_acab_anomaly_scalefactor_1.0.nc ar4_ant_acab_anomaly_scalefactor_2.0.nc
echo "artm.. "
ncap2 -s 'artm(:,:,:)= 2.0 * artm(:,:,:)' ar4_ant_artm_anomaly_scalefactor_1.0.nc ar4_ant_artm_anomaly_scalefactor_2.0.nc



echo "now run spin-up script 'antspin.sh'"
echo

