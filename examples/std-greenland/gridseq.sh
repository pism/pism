#!/bin/bash

NN=6

export REGRIDFILE=
export EXSTEP=500
./spinup.sh $NN const 50000 20 hybrid g20km.nc

export REGRIDFILE=g20km.nc
export EXSTEP=50
./spinup.sh $NN const 5000  10 hybrid g10km.nc

export REGRIDFILE=g10km.nc
export EXSTEP=1
./spinup.sh $NN const 200    5 hybrid  g5km.nc

