#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "rungum.sh ERROR: needs two arguments ... ENDING NOW"
  echo "  format:"
  echo "    rungum.sh PROCS MX"
  echo "  where"
  echo "    PROCS     = 1,2,3,... is number of MPI processes"
  echo "    MX        = number of grid points in x,y directions;  MX -> cell width:"
  echo "                51->8mm, 101->4mm, 201->2mm, 401->1mm, 801->500micron"
  exit
fi

NN="$1"
myMx="$2"

pismexec="pismr"

data=initgum.nc

oname=lab$myMx.nc

grid="-Mx $myMx -My $myMx -Mz 26 -Lz 0.025 -Mbz 0 -Lbz 0 -z_spacing equal"

physics="-config_override gumparams.nc -no_energy -cold -sia_flow_law isothermal_glen -sia_e 1.0 -gradient mahaffy"

endtime=9.5066e-06   # = 300 / 31556926 = 300 s
ts_dt=3.1689e-09     # = 0.1 / 31556926 = 0.1 s
diagnostics="-ts_file ts_$oname -ts_times 0:$ts_dt:$endtime"

mpiexec -n $NN $pismexec -boot_file $data $physics $diagnostics \
    $grid -ys 0.0 -y $endtime -max_dt 1e-8 -o $oname

