#!/bin/bash
# script to play with plastic till SSA and superposition modification of
# EISMINT II experiment I

NN=8

mpiexec -n $NN pisms -eisII I -Mx 61 -My 61 -Mz 201 -y 100000 -o eis2I100k

mpiexec -n $NN pisms -eisII I -if eis2I100k.nc -y 80000 -o eis2I180k

mpiexec -n $NN pisms -eisII I -if eis2I180k.nc -y 10000 -track_Hmelt -f3d -o eis2I190k

mpiexec -n $NN pisms -eisII I -Mx 121 -My 121 -Mz 301 \
   -regrid eis2I190k.nc -regrid_vars TBHh \
   -y 10000 -track_Hmelt -f3d -o eis2I_fine_wmelt

mpiexec -n $NN pisms -eis2Ipl -if eis2I_fine_wmelt.nc -y 10 -f3d -o eis2Ipl10

mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl10.nc -y 90 -f3d -o eis2Ipl100

mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl100.nc -y 900 -f3d -o eis2Ipl1000
    -mato eis2Ipl1000 -matv bcYTHLCQ0345

mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl1000.nc -y 1000 -f3d -o eis2Ipl2000
    -mato eis2Ipl2000 -matv bcYTHLCQ0345

mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl2000.nc -y 3000 -f3d -o eis2Ipl5000
    -mato eis2Ipl5000 -matv bcYTHLCQ0345

mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl5000.nc -y 5000 -f3d -o eis2Ipl10k
    -mato eis2Ipl10k -matv bcYTHLCQ0345

# pisms -eis2Ipl -if eis2I190k.nc -till_phi 0.0,20.0,5.0,0.0 -y 100 -o eis2I190kpl_lakep100


