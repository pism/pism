#!/bin/bash

NN=4

rm -f nbreen.nc
ncgen -o nbreen_config.nc nbreen_config.cdl

# -topg_to_phi 5.0,20.0,-300.0,700.0

#myMx=264
#myMy=207
myMx=133
myMy=104

mpiexec -n $NN pismo -config_override nbreen_config.nc -no_mass -no_energy -hydrology distributed -report_mass_accounting -no_model_strip 1.0 -boot_file pismnbreen.nc -skip -skip_max 10 -Mx $myMx -My $myMy -Mz 121 -z_spacing equal -Lz 600 -max_dt 0.02 -y 1.0 -o nbreen_y1.nc

