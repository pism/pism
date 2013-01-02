#!/bin/bash

NN=4

#dx=125
#myMx=264
#myMy=207

#dx=250
#myMx=133
#myMy=104

dx=500
myMx=67
myMy=52

grid="-Mx $myMx -My $myMy -Mz 11 -z_spacing equal -Lz 600"

physics="-config_override nbreen_config.nc -no_mass -no_energy"

hydro="-hydrology distributed -report_mass_accounting"

pismexec="pismo -no_model_strip 1.0"

YY=5

#FIXME:   For now, runs generate over-large velocities and (locally) too-large water
# depths in areas of thick ice.  Compared to hydrolakes/matlab/nbreenwater.m, that is.
# This behavior may be because topg is different (subsample versus interpolate).
# But here are possible actions:
#  1. Add outline, outside of which W is reset to zero.
#  (that is the behavior in hydrolakes/matlab/doublediff.m)
#  2. [DONE.  At this point the velocities are NOT too large but the too-large
#  water thickness issue remains.]  Avoid differencing bed topography and pressure
#  across the periodic boundary by using a no_model_mask strip near boundary of
#  computational domain.
#  (such differencing is avoided in hydrolakes/matlab/doublediff.m)
#  3. provide the water velocity diagnostically

mpiexec -n $NN $pismexec -boot_file pismnbreen.nc $physics $hydro \
  $grid -max_dt 0.1 -y $YY -o nbreen_y${YY}_${dx}m.nc

