#!/bin/bash

./sg_create_flowline.py
DATANAME=sg_flowline.nc
PISM_DATANAME=pism_$DATANAME
flowline.py -e -o $PISM_DATANAME $DATANAME
# config file
CDLCONFIG=sg.cdl
PCONFIG=sg_config.nc
ncgen -o $PCONFIG $CDLCONFIG
