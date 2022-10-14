#!/bin/bash

#(routingtest.sh) test 'routing' hydrology on Antarctic geometry and constant
#  basal melt rate of 1 cm/year; compare
#  http://www2.gi.alaska.edu/snowice/glaciers/iceflow/bueler-igs-fairbanks-june2012.pdf

# run preprocess.sh before this script

NN=4

DOIT=
#DOIT=echo

PREFIX=
PISMGO="mpiexec -n $NN ${PREFIX}pismr"

VERTGRID="-Lz 5000 -Lbz 2000 -Mz 51 -Mbz 21 -z_spacing equal"

OPTIONS="-sia_e 5.6 -atmosphere given -atmosphere_given_file pism_Antarctica_5km.nc -surface simple -ocean pik -meltfactor_pik 1.5e-2"

HYDRO="-hydrology routing -hydrology_use_const_bmelt -hydrology_const_bmelt 3.1689e-10 -hydrology_hydraulic_conductivity 1.0e-3 -report_mass_accounting"

ENDTIME=20000

dorun () {
  GRID=$1
  LABEL=$2

  #(from antspinCC.sh)  bootstrapping plus short SIA run for 100 years
  cmd="$PISMGO -skip -skip_max 10 -i pism_Antarctica_5km.nc -bootstrap $GRID $VERTGRID $OPTIONS -front_retreat_file pism_Antarctica_5km.nc -y 100 -o pre${LABEL}.nc"
  $DOIT $cmd

  EXTRA="-spatial_file ex_routing${LABEL}.nc -spatial_times 200:100:$ENDTIME -spatial_vars bwat,bwp,bwatvel,subglacial_water_input_rate"

  #hydrology only run for $ENDTIME years
  cmd="$PISMGO -i pre${LABEL}.nc -bootstrap -Lz 5000 $OPTIONS -front_retreat_file pism_Antarctica_5km.nc $HYDRO -no_mass -energy none -stress_balance ssa -max_dt 10.0 -ys 0 -ye $ENDTIME $EXTRA -o routing${LABEL}.nc"
  $DOIT $cmd
}

HUNDREDKMGRID="-Mx 60 -My 60"
FIFTYKMGRID="-Mx 120 -My 120"
TWENTYFIVEKMGRID="-Mx 240 -My 240"
FIFTEENKMGRID="-Mx 400 -My 400"
TENKMGRID="-Mx 600 -My 600"
FIVEKMGRID="-Mx 1200 -My 1200"

# these first three regenerate results from IGS talk:
#dorun "$HUNDREDKMGRID" 100km
#dorun "$FIFTYKMGRID" 50km
dorun "$TWENTYFIVEKMGRID" 25km

# these are more expensive, naturally
#dorun "$FIFTEENKMGRID" 15km
#dorun "$TENKMGRID" 10km
#dorun "$FIVEKMGRID" 5km

