#!/bin/bash

# Copyright (C) 2011-2012 the PISM authors

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

# create usurf, needed by regional-tools not pismo
ncap -O -s 'usurf=thk+topg' $WORKING $WORKING
ncatted -O -a units,usurf,a,c,"m" $WORKING
ncatted -O -a long_name,usurf,a,c,"ice surface elevation" $WORKING

echo "copying geometry fields for boundary conditions in no_model area..."
# create fields thkstore and usrfstore so that pismo is able to appropriately
# assign Dirichlet b.c. for surface gradient & driving stress
ncap -O -s "usurfstore=1.0*usurf" $WORKING $WORKING
ncatted -a standard_name,usurfstore,d,, $WORKING # remove it
ncatted -O -a units,usurfstore,a,c,"m" $WORKING
ncatted -O -a long_name,usurfstore,a,c,"stored ice surface elevation (for regional b.c.)" $WORKING
ncap -O -s "thkstore=1.0*thk" $WORKING $WORKING
ncatted -a standard_name,thkstore,d,, $WORKING # remove it
ncatted -O -a units,thkstore,a,c,"m" $WORKING
ncatted -O -a long_name,thkstore,a,c,"stored ice thickness (for regional b.c.)" $WORKING

echo "... done with cleaning file $WORKING"
echo

