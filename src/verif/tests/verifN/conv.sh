#!/bin/bash

#for N in 46 91 181 361 721  # last case quite slow
for N in 46 91 181 361
do
  ./bodvardsson -snes_fd -da_grid_x $N -snes_monitor -snes_rtol 1.0e-6
done

