# ANTNC.SH  is designed for modelling the current state of Antarctica on a 
# multi-processor machine.  Ice streams are determined by using balance velocities.
# It assumes  init.nc  is present (contains 5km Antarctica data mostly from 
# the British Antarctica Service; write to ffelb@uaf.edu).
# ELB 1/31/07; 3/15/07; 3/25/07; 4/8/07

N=4
echo "  antNC.sh running with $N processors"

#skip past part already done
if [[ 0 -eq 1 ]]
then

echo "  1.  input BAS data, put on 40km grid, and smooth surface  ****************"
mpiexec -n $N obj/pismr -bif init.nc -Mx 141 -My 141 -Mz 101 -Mbz 20 -gk -e 1.2 \
   -verbose -y 10 -o ant10yr_40km -of n
date

echo "  2.  running for 155k year to equilibriate temp (steady geometry)  ********"
echo "      (SIA only)                                                     ********"
mpiexec -n $N obj/pismr -if ant10yr_40km.nc -gk -e 1.2 -no_mass \
   -y 154990 -o ant155k_40km -of n
date

echo "  3.  regrid to 20km grid: interpolate temp and age                      ***"
mpiexec -n $N obj/pismr -bif init.nc -Mx 281 -My 281 -Mz 201 -Mbz 40 -gk -e 1.2 -no_mass\
   -regrid ant155k_40km.nc -regrid_vars TBe -verbose -y 0.01 -o ant155k_20km -of n
date

echo "  4. smoothing surface on fine grid                             ************"
mpiexec -n $N obj/pismr -if ant155k_20km.nc -gk -e 1.2 -verbose -y 9.99 -o ant155k_20kmSM -of n
date

echo "  5. running for another 3000 yrs to equilibriate temp on fine  ************"
echo "     grid (steady geometry and SIA only)                         ************"
mpiexec -n $N obj/pismr -if ant155k_20kmSM.nc -gk -e 1.2 -no_mass\
   -y 2990 -o ant158k_20km -of n
date

echo "  6.  2 yrs to smooth (with MacAyeal)                              ********"
mpiexec -n $N obj/pismr -if ant158k_20km.nc -gk -e 1.2 -mv -ksp_rtol 1e-6\
   -ocean_kill -verbose -y 2 -o ant158k_20kmSMMV -of n
date

echo "  7.  another 500 yrs to equilibriate temp (steady geometry      ********"
echo "      (with MacAyeal)                                             ********"
mpiexec -n $N obj/pismr -if ant158k_20kmSMMV.nc -gk -e 1.2 -no_mass -mv -ksp_rtol 1e-6\
   -ocean_kill -verbose -y 498 -o ant158p5k_20km -of n
date

echo "  8. regrid to 14km grid   ****************"
mpiexec -n $N obj/pismr -bif init.nc -Mx 401 -My 401 -Mz 201 -Mbz 40 -gk -e 1.2 \
   -regrid ant158p5k_20km.nc -regrid_vars TBe -verbose -y 0.001 -o ant158p5k_14km -of n
date

echo "  9. smooth for 2 years on 14km grid             *****************************"
mpiexec -n $N obj/pismr -if ant158p5k_14km.nc -gk -e 1.2 -y 2 -verbose \
   -o ant158p5k_2yr_14km -of n
date

fi

echo "  10. compute everything for 1 year on 14km grid     ****************************"
mpiexec -n $N obj/pismr -if ant158p5k_2yr_14km.nc -gk -e 1.2 -mv -ksp_rtol 1e-6 \
     -ocean_kill -y 1 -verbose -o ant158p5k_final_14km -of n
date
