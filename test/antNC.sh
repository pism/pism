# ANTNC.SH  is designed for modelling the current state of Antarctica on a 
# multi-processor machine.  Ice streams are determined by using balance velocities.
# It assumes  init.nc  is present (contains 5km Antarctica data mostly from 
# the British Antarctica Service; write to ffelb@uaf.edu).  ELB 1/31/07; 3/15/07; 3/25/07

N=4
echo "  antNC.sh running with $N processors"

#skip past part already done
if [[ 0 -eq 1 ]]
then

echo "  1.  input BAS data, put on 40km grid, and smooth surface  ****************"
mpiexec -n $N obj/pismr -bif init.nc -Mx 141 -My 141 -Mz 101 -Mbz 20 -gk -e 1.2 \
   -verbose -y 10 -o ant10yr_40km -of n
date

echo "  2A.  running for 150k year to equilibriate temp (steady geometry)  ********"
echo "      (SIA only)                                                     ********"
mpiexec -n $N obj/pismr -if ant10yr_40km.nc -gk -e 1.2 -no_mass \
   -y 150000 -o ant150k_40km -of n
date

fi

echo "  2B.  running for 10 year to smooth (with MacAyeal)               ********"
mpiexec -n $N obj/pismr -if ant150k_40km.nc -gk -e 1.2 -mv -ksp_rtol 1e-6\
   -ocean_kill -verbose -y 10 -o ant150k_40kmSMMV -of n
date

echo "  2C.  running for 10k year to equilibriate temp (steady geometry)  ********"
echo "      (with MacAyeal)                                               ********"
mpiexec -n $N obj/pismr -if ant150k_40kmSMMV.nc -gk -e 1.2 -no_mass -mv -ksp_rtol 1e-6\
   -ocean_kill -verbose -y 9990 -o ant160k_40km -of n
date

echo "  3.  regrid to 20km grid: interpolate temp and age                      ***"
#mpiexec -n $N obj/pismr -bif init.nc -Mx 281 -My 281 -Mz 201 -Mbz 40 -gk -e 1.2 -no_mass\
#   -regrid ant150k_40km.nc -regrid_vars TBe -verbose -y 0.01 -o ant150k_20km -of n
mpiexec -n $N obj/pismr -if ant160k_40km.nc -gk -e 1.2 -y 0.0001 -o ant160k_40km -of p
mpiexec -n $N obj/pismr -bif init.nc -Mx 281 -My 281 -Mz 201 -Mbz 40 -gk -e 1.2 \
   -regrid ant160k_40km.pb -regrid_vars TBe -verbose -y 0.01 -o ant160k_20km -of n
date

echo "  4. smoothing surface on fine grid                             ************"
mpiexec -n $N obj/pismr -if ant160k_20km.nc -gk -e 1.2 -verbose -y 9.99 -o ant160k_20kmSM -of n
date

echo "  5A. running for another 3000 yrs to equilibriate temp on fine  ************"
echo "     grid (steady geometry and SIA only)                         ************"
mpiexec -n $N obj/pismr -if ant160k_20kmSM.nc -gk -e 1.2 -no_mass\
   -y 2980 -o ant163k_20km -of n
date

echo "  5B.  10 yrs to smooth (with MacAyeal)                              ********"
mpiexec -n $N obj/pismr -if ant163k_20km.nc -gk -e 1.2 -mv -ksp_rtol 1e-6\
   -ocean_kill -verbose -y 10 -o ant163p5k_20kmSMMV -of n
date

echo "  5C.  another 500 yrs to equilibriate temp (steady geometry      ********"
echo "      (with MacAyeal)                                             ********"
mpiexec -n $N obj/pismr -if ant163k_20kmSMMV.nc -gk -e 1.2 -no_mass -mv -ksp_rtol 1e-6\
   -ocean_kill -verbose -y 490 -o ant163p5k_20km -of n
date

#echo "  5C. refining in the vertical and running for 10 yrs smoothing  ************"
#echo "      to equilibriate temp on fine grid (steady geom)            ************"
#mpiexec -n $N obj/pismr -bif init.nc -Mx 281 -My 281 -Mz 401 -Mbz 80 -gk -e 1.2 \
#   -regrid ant153k_20km.nc -regrid_vars TBe -verbose -y 10.0 -o ant154k_20kmVFINE -of n
#date
#echo "  5D. running on fine for 100 yrs with steady geom.  ************"
#mpiexec -n $N obj/pismr -if ant154k_20kmVFINE.nc -gk -e 1.2 -y 100 -no_mass \
#   -o ant154k_20kmVFINE_p100yr -of n

echo "  6. run 100 years w everything               ****************************"
mpiexec -n $N obj/pismr -if ant163p5k_20km.nc -gk -e 1.2 -mv -ksp_rtol 1e-6 -ocean_kill\
   -verbose -y 100 -o ant163p6k_20km -of n
date

echo "  7. another 1400 years w everything             *****************************"
mpiexec -n $N obj/pismr -if ant163p6k_20km.nc -gk -e 1.2 -mv -ksp_rtol 1e-6 -ocean_kill\
   -verbose -y 1400 -o ant165k_20km -of n
date

echo "  8. generate a drag map                                                 *****"
obj/get_drag -if ant165k_20km.nc -balvel init.nc -gk -e 1.2 -verbose -o ant_20km_drag

echo "  9. regrid to 14km grid  (note vertical grid is too coarse!)  ****************"
#mpiexec -n $N obj/pismr -bif init.nc -Mx 401 -My 401 -Mz 51 -Mbz 10 -gk -e 1.2 \
#   -regrid ant155k_20km.nc -regrid_vars TBe -verbose -y 0.001 -o ant155k_14kmTHIN -of n
mpiexec -n $N obj/pismr -if ant165k_20km.nc -gk -e 1.2 -y 0.0001 -o ant165k_20km -of p
mpiexec -n $N obj/pismr -bif init.nc -Mx 401 -My 401 -Mz 51 -Mbz 10 -gk -e 1.2 \
   -regrid ant165k_20km.pb -regrid_vars TBe -verbose -y 0.001 -o ant165k_14kmTHIN -of n
date

echo "  10. smooth for 10 years on 14km grid             *****************************"
mpiexec -n $N obj/pismr -if ant165k_14kmTHIN.nc -gk -e 1.2 -y 10 -verbose \
   -no_temp -o ant165k_p10yr_14kmTHIN -of n
date

echo "  11. generate a drag map                            ****************************"
obj/get_drag -if ant165k_p10yr_14kmTHIN.nc -balvel init.nc -gk -e 1.2 -verbose -o ant_14km_drag

echo "  12. compute everything for 1 year on 14km grid     ****************************"
mpiexec -n $N obj/pismr -if ant165k_p10yr_14kmTHIN.nc -gk -e 1.2 -mv -ksp_rtol 1e-6 \
     -ocean_kill -y 1 -verbose -no_temp -o ant165k_p10yr_mv1yr_14kmTHIN -of n
date
