#!/bin/bash

./spinup.sh 4 const 30000 20 hybrid g20km.nc &> out.gridseq

REGRIDFILE=g20km.nc EXSTEP=10 ./spinup.sh 4 const 10000 10 hybrid g10km.nc &>> out.gridseq

REGRIDFILE=g10km.nc EXSTEP=1 ./spinup.sh 4 const 100 5 hybrid g5km.nc &>> out.gridseq

