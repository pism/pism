#!/bin/bash

./sg_create_flowline.py
DATANAME=storglaciaren_flowline.nc
PISM_DATANAME=pism_$DATANAME
flowline.py -e -o $PISM_DATANAME $DATANAME


./sg_create_3d.py
ncap2 -O -s "where(thk<1) ftt_mask=0;" pism_storglaciaren_3d.nc pism_storglaciaren_mask.nc

# config file
CDLCONFIG=psg_config.cdl
PCONFIG=psg_config.nc
ncgen -o $PCONFIG $CDLCONFIG
