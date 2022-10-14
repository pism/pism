#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "rungum.sh ERROR: needs two arguments ... ENDING NOW"
  echo "  format:"
  echo "    rungum.sh PROCS MX"
  echo "  where"
  echo "    PROCS     = 1,2,3,... is number of MPI processes"
  echo "    MX        = number of grid points in x,y directions;  MX -> cell width:"
  echo "                52 -> 10mm,  104 -> 5mm, 208 -> 2.5mm, 520 -> 1mm"
  exit
fi

NN="$1"
myMx="$2"

pismexec="pismr"

oname=lab$myMx.nc
initfile=init$oname

grid="-Mx $myMx -My $myMx -Mz 101 -Lz 0.1 -Mbz 0 -Lbz 0 -z_spacing equal"

climate="-surface given -surface_given_file $initfile -surface.given.smb_max 1e10"

physics="-config_override gumparams.nc -energy none -sia_flow_law isothermal_glen -sia_e 1.0 -gradient mahaffy"

endtime=746s                    # Sayag personal communication

ts_dt=0.1s
timediag="-ts_file ts_$oname -ts_times $ts_dt"

ex_dt=10s
exvars="diffusivity,flux_mag,velbar_mag,velsurf_mag,mask,thk,wvelsurf"
exdiag="-spatial_file ex_$oname -spatial_vars $exvars -spatial_times $ex_dt"

dt="-time_stepping.resolution 1e-6 -max_dt $ts_dt"

mpiexec -n $NN $pismexec -i $initfile -bootstrap $grid $climate $physics \
    $timediag $exdiag -ys 0.0 -y $endtime $dt -o $oname
