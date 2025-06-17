#!/bin/bash

python preprocess_sym.py
ncgen -o config.nc config.cdl


N=8
run_length=20000
sb="ssa+sia"
resolution="8km"
out=g${resolution}_${sb}_${run_length}a.nc
mpirun -np $N pism -config_override config.nc \
       -stress_balance.model $sb \
       -grid.dx $resolution \
       -grid.dy $resolution \
       -i mismip+_sym.nc \
       -input.bootstrap yes \
       -o_size medium \
       -surface.given.file climate_sym.nc \
       -geometry.front_retreat.prescribed.file mismip+_sym.nc \
       -output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
       -output.extra.file spatial_$out \
       -output.extra.times 100year \
       -output.extra.vars bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag,ice_mass_transport_across_grounding_line -output.file state_$out \
       -time.run_length $run_length


