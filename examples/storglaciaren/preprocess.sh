#!/bin/bash

./sg_create_flowline.py
DATANAME=storglaciaren_flowline.nc
PISM_DATANAME=pism_$DATANAME
flowline.py -e -o $PISM_DATANAME $DATANAME


./sg_create_3d.py
ncap2 -O -s "where(thk>0) ftt_mask=0;" pism_storglaciaren_3d.nc pism_storglaciaren_mask.nc
ncap2 -O  -s "ice_area_specific_volume=0*ftt_mask;  tillwat=0*ftt_mask;" pism_storglaciaren_3d.nc  pism_storglaciaren_3d.nc
ncatted -a units,ice_area_specific_volume,o,c,"m3/m2" -a units,tillwat,o,c,"m" pism_storglaciaren_3d.nc

./create_warming_climate.py sg_warming_1K.nc

# config file
CDLCONFIG=psg_config.cdl
PCONFIG=psg_config.nc
ncgen -o $PCONFIG $CDLCONFIG
