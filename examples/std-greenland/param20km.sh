#!/bin/bash

NN=2
DURA=100000

for PPQ in 0.1 0.25 0.8 ; do
  for TEFO in 0.01 0.02 0.05 ; do
     PARAM_PPQ=$PPQ PARAM_TEFO=$TEFO ./spinup.sh $NN const $DURA 20 hybrid g20km_${PPQ}_${TEFO}_SGL.nc
     PARAM_PPQ=$PPQ PARAM_TEFO=$TEFO PARAM_NOSGL=foo ./spinup.sh $NN const $DURA 20 hybrid g20km_${PPQ}_${TEFO}_NOSGL.nc
  done
done
