#!/bin/bash

set -u
set -x
set -e

pism -eisII A \
   -Lz 5000 \
   -Mx 3 \
   -My 3 \
   -Mz 3 \
   -o eis2-a.nc \
   -y 1e4 \
   -output.snapshot.file ${filename:-snapshot.nc} \
   -output.snapshot.times 2000 \
   -output.snapshot.split ${split:-no} \
   ""
