#!/bin/bash

# show it runs for low values; Mbz=0,1 will not even give output file with
# a litho_temp var
for MM in 0 1 2 3
do
  pism_btutest -ys 0.0 -ye 1.0 -dt 0.1 -Mbz $MM -Lbz 1000 -o btuout_Mbz=${MM}.nc
done

# refine spatially only; five cases
for MM in 11 21 41 81 161
do
  echo "case [-dt 0.1 AND] -Mbz ${MM}:"
  pism_btutest -ys 0.0 -ye 1.0 -dt 0.1 -Mbz $MM -Lbz 1000 -verbose 1
done

# refine dt only
for DD in 0.8 0.4 0.2 0.1 0.05
do
  echo "case [-Mbz 41 AND] -dt ${DD}:"
  pism_btutest -ys 10.0 -ye 110.0 -dt $DD -Mbz 41 -Lbz 1000 -verbose 1
done


