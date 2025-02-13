#!/bin/bash

cd examples/std-greenland/           # g10km_gridseq.nc should be in this directory
mkdir paramstudy
cd paramstudy
ln -s ../g10km_gridseq.nc .            # these four lines make links to ...
ln -s ../pism_Greenland_5km_v1.1.nc .  #
ln -s ../spinup.sh .                   #
ln -s ../param.sh .                    # ... existing files in examples/std-greenland/
./param.sh
