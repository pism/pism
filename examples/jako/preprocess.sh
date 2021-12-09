#!/bin/bash

# Copyright (C) 2011-2014, 2017, 2019, 2021 the PISM authors

# downloads SeaRISE "1km Greenland data set" NetCDF file,
# adjusts metadata, saves under new name with fields needed
# for pism_regional.py

# depends on wget and NCO (ncrename, ncap, ncwa, ncks)

set -e  # exit on error

# get file; see page http://websrv.cs.umt.edu/isis/index.php/1km_Greenland_data_set

DATAURL=http://websrv.cs.umt.edu/isis/images/a/ab/
DATANAME=Greenland1km.nc
DATASIZE=80Mb
echo "fetching $DATASIZE master file $DATANAME ... "
wget -nc ${DATAURL}${DATANAME}

WORKING=gr1km.nc
echo "creating PISM-readable file $WORKING from master ..."

# copies over preserving history and global attrs; drops climate data:
ncks -O -v x,y,mapping,bheatflx,topg,thk $DATANAME $WORKING

# remove time dimension
ncwa -O -a t $WORKING $WORKING

echo "adding lat and lon fields by using nc2cdo.py (which is in pism/util/)"
nc2cdo.py $WORKING

# create usurf, needed by regional-tools not pism
ncap2 -O -s 'usurf=thk+topg' $WORKING $WORKING
ncap2 -O -s 'where(usurf<0.0) usurf=0.0' $WORKING $WORKING
ncatted -a standard_name,usurf,d,, $WORKING # remove it
ncatted -O -a units,usurf,o,c,"m" $WORKING
ncatted -O -a long_name,usurf,o,c,"ice surface elevation" $WORKING

echo "copying geometry fields for boundary conditions in no_model area..."
# create fields thkstore and usrfstore so that pism is able to appropriately
# assign Dirichlet b.c. for surface gradient & driving stress
ncap2 -O -s "usurfstore=1.0*usurf" $WORKING $WORKING
ncatted -a standard_name,usurfstore,d,, $WORKING # remove it
ncatted -O -a units,usurfstore,a,c,"m" $WORKING
ncatted -O -a long_name,usurfstore,a,c,"stored ice surface elevation (for regional b.c.)" $WORKING
ncap2 -O -s "thkstore=1.0*thk" $WORKING $WORKING
ncatted -a standard_name,thkstore,d,, $WORKING # remove it
ncatted -O -a units,thkstore,a,c,"m" $WORKING
ncatted -O -a long_name,thkstore,a,c,"stored ice thickness (for regional b.c.)" $WORKING

echo "... done with cleaning file $WORKING"
echo


# get file containing surface mass balance and other data on 5km grid
# see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATAURL=http://websrv.cs.umt.edu/isis/images/a/a5/
DATANAME=Greenland_5km_v$DATAVERSION.nc
echo "fetching 5km SeaRISE data file which contains surface mass balance ... "
wget -nc ${DATAURL}${DATANAME}
CLIMATEFILE=g5km_climate.nc
echo "creating PISM-readable climate file $CLIMATEFILE from airtemp2m and smb in data file ..."
ncks -O -v mapping,smb,airtemp2m $DATANAME $CLIMATEFILE
ncrename -O -v airtemp2m,ice_surface_temp $CLIMATEFILE
ncatted -O -a units,ice_surface_temp,a,c,"Celsius" $CLIMATEFILE
# convert SMB from liquid water equivalent thickness per year to [kg m-2 year-1];
# assume water density of 1000.0 [kg m-3]
ncap2 -O -s "climatic_mass_balance=1000.0*smb" \
      -s 'climatic_mass_balance@standard_name="land_ice_surface_specific_mass_balance_flux"' \
      -s 'climatic_mass_balance@units="kg m-2 year-1"' \
      $CLIMATEFILE $CLIMATEFILE
ncks -O -x -v smb $CLIMATEFILE $CLIMATEFILE
echo "... done"
echo


# get locally-generated or pre-computed PISM result
echo "checking for locally-generated, or fetching, pre-computed PISM whole ice-sheet result on 5km grid"
URL=https://raw.githubusercontent.com/pism/example-inputs/main/jako/
WHOLE=g5km_gridseq.nc
wget -nc ${URL}/$WHOLE
BCFILE=g5km_bc.nc
echo "creating PISM-readable boundary conditions file $BCFILE from whole ice sheet result ..."
ncks -O -v u_ssa,v_ssa,bmelt,tillwat,enthalpy,litho_temp $WHOLE $BCFILE
# rename bmelt and u_ssa and v_ssa so that they are used as b.c.
ncrename -O -v bmelt,basal_melt_rate_grounded -v u_ssa,u_bc -v v_ssa,v_bc $BCFILE
echo "... done"
echo
