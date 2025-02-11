#!/bin/bash

export PARAM_PPQ=0.5
export REGRIDFILE=g20km_10ka_hy.nc
./spinup.sh 8 paleo 25000 20 hybrid g20km_25ka_paleo.nc &> out.g20km_25ka_paleo &
