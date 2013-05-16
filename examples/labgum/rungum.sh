#!/bin/bash

if [ $# -lt 1 ] ; then
  echo "rungum.sh ERROR: needs one argument ... ENDING NOW"
  echo "  format:"
  echo "    rungum.sh PROCS"
  echo "  where"
  echo "    PROCS     = 1,2,3,... is number of MPI processes"
  exit
fi

NN="$1"

FIXME

pismexec="pismr"

data=initgum.nc

grid="-Mx $myMx -My $myMy -Mz 11 -z_spacing equal -Lz 600"

climate="-surface given -surface_given_file $data"

physics="-config_override nbreen_config.nc -no_mass -no_energy"

diagnostics="-extra_file extras_$oname -extra_times $etimes -extra_vars $evarlist"

mpiexec -n $NN $pismexec -boot_file $data $climate $physics $hydro \
    $grid -max_dt $dtmax -ys 0.0 -y $YY $diagnostics -o $oname

