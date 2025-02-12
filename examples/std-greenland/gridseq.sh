#!/bin/bash
NN=8
export PARAM_PPQ=0.5

export REGRIDFILE=g20km_10ka_hy.nc EXSTEP=100
./spinup.sh $NN const 2000 10 hybrid g10km_gridseq.nc

export REGRIDFILE=g10km_gridseq.nc EXSTEP=10
export PISM_EXEC="pism -stress_balance.sia.max_diffusivity 1000"
./spinup.sh $NN const 30days 5 hybrid g5km_smoothing.nc

export PISM_EXEC=pism
export REGRIDFILE=g5km_smoothing.nc EXSTEP=10
./spinup.sh $NN const 200 5 hybrid g5km_gridseq.nc
