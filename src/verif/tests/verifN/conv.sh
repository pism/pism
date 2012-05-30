#!/bin/bash

# run as 
#   $ ./conv.sh >& conv.txt

MPIDO="mpiexec -n 1"

# dx =   30km 10km  5km  3km  2km  1km  500m  300m
for N in   16   46   91  151  226  451   901  #1501
do
  # second-order upwinding (default)
  #./bodvardsson -snes_fd -da_grid_x $N -snes_monitor -snes_rtol 1.0e-8
  # first-order upwinding:
  #  ./bodvardsson -snes_fd -da_grid_x $N -snes_monitor -snes_rtol 1.0e-8 -bod_up_one

  # second-order upwinding (default) with exact soln as initial guess:
  $MPIDO ./bodvardsson -snes_fd -da_grid_x $N -snes_monitor -snes_rtol 1.0e-8 -bod_exact_init
done

