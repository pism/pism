#!/bin/bash
# script to play with EISMINT II experiment I and plastic till SSA and superposition

NN=8

mpiexec -n $NN pisms -eisII I -Mx 61 -My 61 -Mz 201 -y 100000 -o eis2I100k

mpiexec -n $NN pisms -eisII I -if eis2I100k.nc -y 80000 -o eis2I180k

mpiexec -n $NN pisms -eisII I -if eis2I180k.nc -y 10000 -track_Hmelt -f3d -o eis2I190k

mpiexec -n $NN pisms -eisII I -Mx 121 -My 121 -Mz 301 \
   -regrid eis2I190k.nc -regrid_vars TBHh \
   -y 10000 -track_Hmelt -f3d -o eis2I_fine_wmelt

mpiexec -n $NN pisms -eisII I -if eis2I_fine_wmelt.nc -ssa -plastic -super \
    -track_Hmelt -till_phi 20.0,5.0 -max_low_temps 200 -y 5 -f3d -o eis2I_plastic5

mpiexec -n $NN pisms -eisII I -if eis2I_plastic5.nc -ssa -plastic -super \
    -track_Hmelt -till_phi 20.0,5.0 -max_low_temps 200 -y 5 -f3d -o eis2I_plastic10

mpiexec -n $NN pisms -eisII I -if eis2I_plastic10.nc -ssa -plastic -super \
    -track_Hmelt -till_phi 20.0,5.0 -max_low_temps 200 -y 10 -f3d -o eis2I_plastic20

mpiexec -n $NN pisms -eisII I -if eis2I_plastic20.nc -ssa -plastic -super \
    -track_Hmelt -till_phi 20.0,5.0 -max_low_temps 200 -y 80 -f3d -o eis2I_plastic100 \
    -mato eis2I_plastic100 -matv bcYTHLCQ0345


