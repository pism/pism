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

echo "fetching pre-computed whole ice-sheet result on 5km grid"
URL=http://www.pism-docs.org/download
WHOLE=g5km_0_ftt.nc
wget -nc ${URL}/$WHOLE
echo "... done."
echo

WORKING=gr1km.nc
echo "creating PISM-readable file $WORKING from master ..."

# copies over preserving history and global attrs:
ncks -O $DATANAME $WORKING

# remove time dimension
ncwa -O -a t $WORKING $WORKING

# rename climate fields to work for "-atmosphere given -surface pdd"
ncap -O -s "precipitation=presprcp*(1000.0/910.0)" $WORKING $WORKING
ncatted -O -a units,precipitation,a,c,"m a-1" $WORKING
ncatted -a standard_name,precipitation,d,, $WORKING  # remove standard_name; none available
ncatted -O -a long_name,precipitation,a,c,"ice-equivalent mean annual precipitation rate" $WORKING
ncatted -a standard_name,bheatflx,d,, $WORKING  # remove standard_name; none available
ncks -O -v x,y,bheatflx,topg,thk,precipitation,presartm,mapping $WORKING $WORKING

# create usurf, needed by regional-tools not pismo
ncap -O -s 'usurf=thk+topg' $WORKING $WORKING
ncatted -O -a units,usurf,a,c,"m" $WORKING
ncatted -O -a long_name,usurf,a,c,"ice surface elevation" $WORKING

echo "copying geometry fields as boundary conditions in no_model area..."
# create fields thkstore and usrfstore so that pismo is able to appropriately
# assign Dirichlet b.c. for surface gradient & driving stress
ncap -O -s "usurfstore=1.0*usurf" $WORKING $WORKING
ncatted -a standard_name,usurfstore,d,, $WORKING # remove it
ncatted -O -a units,usurfstore,a,c,"m" $WORKING
ncatted -O -a long_name,usurfstore,a,c,"stored ice surface elevation for regional boundary condition" $WORKING
ncap -O -s "thkstore=1.0*thk" $WORKING $WORKING
ncatted -a standard_name,thkstore,d,, $WORKING # remove it
ncatted -O -a units,thkstore,a,c,"m" $WORKING
ncatted -O -a long_name,thkstore,a,c,"stored ice thickness for regional boundary condition" $WORKING

echo "... done with cleaning file $WORKING"
echo

BCFILE=g5km_bc.nc
echo "creating PISM-readable file $BCFILE from master ..."
ncks -v u_ssa,v_ssa,bmelt,enthalpy $WHOLE $BCFILE
# rename u_ssa and v_ssa so that they are specified as b.c.
ncrename -O -v u_ssa,u_ssa_bc -v v_ssa,v_ssa_bc $BCFILE
echo "... done with creating bc file $BCFILE"
echo


echo "now do 'python pism_regional.py' and open $WORKING"
echo

