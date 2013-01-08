#!/bin/bash

# format:
#   run.sh PROCS GRID DURATION
# so example usage is
#   $ ./run.sh 4 500 5 >& out.nbreen_y5_500m &

if [ $# -lt 3 ] ; then
  echo "run.sh ERROR: needs three arguments"
  echo "  example usage: 'run.sh 4 500 5' for 4 processors and 500 m grid and 5 year run"
  exit
fi

NN="$1"

if [ "$2" -eq "500" ]; then
  dx=500
  myMx=67
  myMy=52
elif [ "$2" -eq "250" ]; then
  dx=250
  myMx=133
  myMy=104
elif [ "$2" -eq "125" ]; then
  dx=125
  myMx=264
  myMy=207
else
  echo "invalid second argument: must be in {125,250,500}"
  exit
fi

YY="$3"

# the following is what I want but issue #125 arises:
#pismexec="pismo -no_model_strip 1.0"

pismexec="pismr"

data=pismnbreen.nc

grid="-Mx $myMx -My $myMy -Mz 11 -z_spacing equal -Lz 600"

climate="-surface given -surface_given_file $data"

physics="-config_override nbreen_config.nc -no_mass -no_energy"

runpism () {
  mpiexec -n $NN $pismexec -boot_file $data $climate $physics $hydro \
    $grid -max_dt 0.1 -y $YY $diagnostics -o $oname
}

# lakes run: very fast
oname=nbreen_y${YY}_${dx}m_lakes.nc
diagnostics="-extra_file extras_$oname -extra_times 0:0.1:$YY -extra_vars bmelt,bwat,bwp,bwatvel"
hydro="-hydrology lakes -report_mass_accounting"
runpism

# distributed run
oname=nbreen_y${YY}_${dx}m.nc
diagnostics="-extra_file extras_$oname -extra_times 0:0.1:$YY -extra_vars cbase,bmelt,bwat,bwp,bwatvel"
hydro="-hydrology distributed -report_mass_accounting -ssa_sliding -ssa_dirichlet_bc"
runpism

