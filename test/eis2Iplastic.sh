#!/bin/bash
# script for plastic till SSA and superposition modification of
# EISMINT II experiment I

NN=8

mpiexec -n $NN pisms -eisII I -Mx 61 -My 61 -Mz 201 -y 100000 -o eis2I100k

mpiexec -n $NN pisms -eisII I -if eis2I100k.nc -y 80000 -o eis2I180k

mpiexec -n $NN pisms -eisII I -if eis2I180k.nc -y 10000 -track_Hmelt \
   -f3d -o eis2I190k

mpiexec -n $NN pisms -eisII I -Mx 121 -My 121 -Mz 251 \
   -regrid eis2I190k.nc -regrid_vars LTBHh \
   -y 10000 -track_Hmelt -f3d -o eis2I_final \
   -mato eis2Iend -matv bcYTHLCQ0345

# starting basic experiment P1
mpiexec -n $NN pisms -eis2Ipl -if eis2I_final.nc -y 10 -f3d \
    -o eis2Ipl10 -mato eis2Ipl10 -matv bcYTHLCQ0345

mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl10.nc -y 90 -f3d \
    -o eis2Ipl100 -mato eis2Ipl100 -matv bcYTHLCQ0345

mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl100.nc -y 900 -f3d \
    -o eis2Ipl1000 -mato eis2Ipl1000 -matv bcYTHLCQ0345

mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl1000.nc -y 1000 -f3d \
    -o eis2Ipl2000 -mato eis2Ipl2000 -matv bcYTHLCQ0345

# result eis2Ipl5000.[nc|m] is experiment P1:
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl2000.nc -y 3000 -f3d \
    -o eis2Ipl5000 -mato eis2Ipl5000 -matv bcYTHLCQ0345
# re speed: marmaduke, with 8 cores, the last 3000 model year run took 
# about 14 hours, so about       5 hours/(1000 model years)
# or, optimistically, about      ( 1 hour/(1000 model years) )/core

# possible stopping point 
#exit
# continuing with experiments P2,P3,P4:

# result eis2IplP2.[nc|m] is experiment P2 (narrower stream):
mpiexec -n $NN pisms -eis2Ipl -if eis2I_final.nc -y 5000 -f3d \
    -stream_width 50.0 -o eis2IplP2 \
    -mato eis2IplP2 -matv bcYTHLCQ0345

# result eis2IplP3.[nc|m] is experiment P3 (stronger downstream till):
mpiexec -n $NN pisms -eis2Ipl -if eis2I_final.nc -y 5000 -f3d \
    -till_phi 20.0,20.0,5.0,8.0,0.0 -o eis2IplP3 \
    -mato eis2IplP3 -matv bcYTHLCQ0345

# result eis2IplP4.[nc|m] is experiment P4 (lake):
mpiexec -n $NN pisms -eis2Ipl -if eis2I_final.nc -y 5000 -f3d \
    -till_phi 0.0,20.0,5.0,5.0,0.0 -o eis2IplP4 \
    -mato eis2IplP4 -matv bcYTHLCQ0345

# possible stopping point 
#exit
# continuing with grid variation experiments P6,P7,P8:

# result eis2IplP6.[nc|m] is experiment P6 (coarser horizontal grid):
mpiexec -n $NN pisms -eis2Ipl -Mx 61 -My 61 -Mz 251 -y 5000 -f3d \
    -regrid eis2I_final.nc -regrid_vars HTBL -o eis2IplP6 \
    -mato eis2IplP6 -matv bcYTHLCQ0345

# result eis2IplP7.[nc|m] is experiment P7 (finer horizontal 7.5km grid):
mpiexec -n $NN pisms -eis2Ipl -Mx 201 -My 201 -Mz 251 -y 5000 -f3d \
    -regrid eis2I_final.nc -regrid_vars HTBL -o eis2IplP7 \
    -mato eis2IplP7 -matv bcYTHLCQ0345

# result eis2IplP8.[nc|m] is experiment P8 (finer vertical grid):
mpiexec -n $NN pisms -eis2Ipl -Mx 121 -My 121 -Mz 501 -y 5000 -f3d \
    -regrid eis2I_final.nc -regrid_vars HTBL -o eis2IplP8 \
    -mato eis2IplP8 -matv bcYTHLCQ0345

# possible stopping point 
#exit
# continuing with long run P5 of 100k model years:

mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl5000.nc -y 5000 -f3d \
    -o eis2Ipl10k -mato eis2Ipl10k -matv bcYTHLCQ0345

# lots of runtime!:
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl10k.nc -y 10000 -f3d \
    -o eis2Ipl20k
    
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl20k.nc -y 10000 -f3d \
    -o eis2Ipl30k
    
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl30k.nc -y 10000 -f3d \
    -o eis2Ipl40k
    
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl40k.nc -y 10000 -f3d \
    -o eis2Ipl50k
    
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl50k.nc -y 10000 -f3d \
    -o eis2Ipl60k
    
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl60k.nc -y 10000 -f3d \
    -o eis2Ipl70k
    
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl70k.nc -y 10000 -f3d \
    -o eis2Ipl80k
    
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl80k.nc -y 10000 -f3d \
    -o eis2Ipl90k

# result eis2Ipl100k.[nc|m] is experiment P5:
mpiexec -n $NN pisms -eis2Ipl -if eis2Ipl90k.nc -y 10000 -f3d \
    -o eis2Ipl100k -mato eis2Ipl100k -matv bcYTHLCQ0345

