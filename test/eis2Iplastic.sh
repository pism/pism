# script to play with EISMINT II experiment I and plastic till; for rev 195 and later

NN = 4

mpiexec -n $NN pisms -eisII I -Mx 61 -My 61 -Mz 201 -y 100000 -o eis2I100k

mpiexec -n $NN pisms -eisII I -if eis2I100k.nc -y 80000 -o eis2I180k

mpiexec -n $NN pisms -eisII I -Mx 121 -My 121 -Mz 301 -regrid eis2I180k.nc -regrid_vars TBHh \
   -y 20000 -track_Hmelt -o eis2I_fine_wmelt -full3Dout

mpiexec -n $NN pisms -eisII I -if eis2I_fine_wmelt.nc -ssa -plastic -super -verbose -y 5 \
   -o eis2I_plastic5 -full3Dout

mpiexec -n $NN pisms -eisII I -if eis2I_plastic5.nc -ssa -plastic -super -verbose -y 5 \
   -o eis2I_plastic10 -full3Dout

