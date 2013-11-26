#!/bin/bash

NN=2
DUR=10000

for PPQ in 0.1 0.5 1.0 ; do
  for SIAE in 1 3 9 ; do
     PARAM_PPQ=$PPQ PARAM_SIAE=$SIAE ./spinup.sh $NN const $DUR 20 hybrid par20km_${PPQ}_${SIAE}.nc
  done
done
