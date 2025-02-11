#!/bin/bash

PISM_DO=echo PARAM_PPQ=0.5 ./spinup.sh 8 const 10000 20 hybrid g20km_10ka_hy.nc
