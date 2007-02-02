#!/bin/sh
# MPIVERIFYNOW is a script to do fairly quick verification of the isothermal and 
# thermocoupled shallow ice components and of the ice stream components of PISM on
# a multi-processor computer.  Uses tests C, I and G.  It is intended to do roughly
# the minimal amount of computation to show convergence to continuum results.   ELB 1/31/07

N=4
echo "mpiverifynow.sh running with $N processors"

echo "++++++++ starting with isothermal shallow ice approx (SIA) test C and Mx = My = 41, 61, 81, 121"
echo "          (so  dx = dy = 50, 33.3, 25, 16.7 km):"

for myMx in 41 61 81 121
do
   time mpiexec -n $N obj/pismv -test C -Mx $myMx -My $myMx -Mz 31 -y 15208.0 -verbose > _temp_result.txt 
   sed '/  history =/,+1!d' _temp_result.txt | sed 1d
   sed '/Actual ERRORS/,+2!d' _temp_result.txt
   date
done

echo "+++++++++ continuing with ice stream test I and My=49 , 193, 769, 3073, 12289"
echo "          (so  dy = 5000, 1250, 312.5, 78.125, 19.53 *meters*):"

for myMy in 49 193 769 3073 12289
do
   time mpiexec -n $N obj/pismv -test I -Mx 5 -My $myMy -mv_rtol 1e-7 -verbose > _temp_result.txt 
   sed '/  history =/,+1!d' _temp_result.txt | sed 1d
   sed '/Actual ERRORS/,+2!d' _temp_result.txt
   date
done

echo "+++++++++ continuing with thermocoupled SIA test G and Mx = My = 61, 91, 121, 181"
echo "          (so  dx = dy = 30, 20, 15, 10  km  and dz = 66.7, 44.4, 33.3, 22.2 m):"

for myMx in 61 91 121 181
do
   time mpiexec -n $N obj/pismv -test G -Mx $myMx -My $myMx -Mz $myMx -y 25000.0 -verbose > _temp_result.txt 
   sed '/  history =/,+1!d' _temp_result.txt | sed 1d
   sed '/Actual ERRORS/,+4!d' _temp_result.txt
   date
done
