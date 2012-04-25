#!/bin/bash

# Copyright (C) 2011-2012 the PISM authors

# downloads SeaRISE "1km Greenland data set" NetCDF file,
# adjusts metadata, saves under new name with fields needed for pism_regional.py

# depends on wget and NCO (ncrename, ncap, ncwa, ncks)

set -e  # exit on error

# get file; see page http://websrv.cs.umt.edu/isis/index.php/1km_Greenland_data_set

DATAURL=http://websrv.cs.umt.edu/isis/images/a/ab/
DATANAME=Greenland1km.nc
DATASIZE=80Mb

echo "fetching $DATASIZE master file $DATANAME ... "
wget -nc ${DATAURL}${DATANAME}
echo "  ... done."
echo

WORKING=gr1km.nc
echo "creating simplified file $WORKING with vars x,y,usurf,thk from master ..."

# copies over preserving history and global attrs:
ncks -O -v topg,thk,x,y $DATANAME $WORKING

# remove time dimension
ncwa -O -a t $WORKING $WORKING
# create usurf
ncap -O -s 'usurf=thk+topg' $WORKING $WORKING

echo "done with cleaning file $WORKING"
echo "now do 'python pism_regional.py' and open $WORKING"
echo

