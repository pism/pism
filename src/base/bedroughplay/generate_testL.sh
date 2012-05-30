#!/bin/bash

set -e  # exit on error
set -x  # see commands as they are issued

for M in 21 41 81 161;
do
  pismv -test L -Mx $M -My $M -Mz 3 -Mbz 1 -y 5000.0 -verbose 1 -eo -o testL_$M.nc
  ncks -O -v x,y,climatic_mass_balance,topg,thk,usurf testL_$M.nc testL_$M.nc
done

