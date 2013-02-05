#!/bin/bash

if [ $# -lt 4 ] ; then
  echo "run.sh ERROR: needs four arguments ... ending"
  echo "  format:"
  echo "    run.sh PROCS GRID DURATION TYPE"
  echo "  where"
  echo "    PROCS    =1,2,3,... is number of MPI processes"
  echo "    GRID     is in {500, 250, 125}, the grid spacing in meters"
  echo "    DURATION is run duration in years"
  echo "    TYPE     is in {dist, event, lakes}"
  echo "  example usage:"
  echo "    run.sh 4 500 5 dist"
  echo "  i.e. 4 processors, 500 m grid, 5 model year run, and '-hydrology distributed'"
  exit
fi

NN="$1"

if [ "$2" -eq "500" ]; then
  dx=500
  dtmax=0.1
  myMx=67
  myMy=52
elif [ "$2" -eq "250" ]; then
  dx=250
  dtmax=0.1
  myMx=133
  myMy=104
elif [ "$2" -eq "125" ]; then
  dx=125
  dtmax=0.01  # more frequent just because so many hydrology substeps occur
  myMx=264
  myMy=207
else
  echo "invalid second argument: must be in {125,250,500}"
  exit
fi

YY="$3"

# these diagnostics apply to "dist" and "event":
evarlist="thk,cbase,bmelt,hydroinput,bwat,bwp,bwatvel,bwprel,effbwp,enwat"

if [ "$4" = "dist" ]; then

  # distributed run
  oname=nbreen_y${YY}_${dx}m_dist.nc
  hydro="-hydrology distributed -hydrology_null_strip 1.0 -report_mass_accounting -ssa_sliding -ssa_dirichlet_bc"
  etimes="0:0.1:$YY"

elif [ "$4" = "event" ]; then

  # distributed run with summer event
  oname=nbreen_y${YY}_${dx}m_event.nc
  hydro="-hydrology distributed -hydrology_null_strip 1.0 -report_mass_accounting -ssa_sliding -ssa_dirichlet_bc -input_to_bed_file fakesummerevent.nc -input_to_bed_period 1.0 -input_to_bed_reference_year 0.0"
  etimes="0.0:0.005:$YY"
#FIXME: this produced a bug: it gave warning about more than 500 frames
# etimes="0.0:daily:$YY"

elif [ "$4" = "lakes" ]; then

  # lakes run: very fast
  oname=nbreen_y${YY}_${dx}m_lakes.nc
  hydro="-hydrology lakes -hydrology_null_strip 1.0 -report_mass_accounting -hydrology_hydraulic_conductivity_at_large_W 1.0e-3"
  evarlist="thk,bmelt,hydroinput,bwat,bwp,bwatvel"  # revised
  etimes="0:0.1:$YY"

else
  echo "invalid fourth argument; must be in {dist,event,lakes}"
  exit
fi

pismexec="pismr"

data=pismnbreen.nc

grid="-Mx $myMx -My $myMy -Mz 11 -z_spacing equal -Lz 600"

climate="-surface given -surface_given_file $data"

physics="-config_override nbreen_config.nc -no_mass -no_energy"

diagnostics="-extra_file extras_$oname -extra_times $etimes -extra_vars $evarlist"

mpiexec -n $NN $pismexec -boot_file $data $climate $physics $hydro \
    $grid -max_dt $dtmax -ys 0.0 -y $YY $diagnostics -o $oname

