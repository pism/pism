# ANTNC.SH  is designed for modelling the current state of Antarctica on a 
# multi-processor machine.  Ice streams are determined by using balance velocities.
# It assumes  init.nc  is present (contains 5km Antarctica data mostly from 
# the British Antarctica Service; write to ffelb@uaf.edu).  ELB 1/31/07; 3/15/07

N=4
echo "  antNC.sh running with $N processors"

#skip past part already done
if [[ 0 -eq 1 ]]
then

echo "  1.  input BAS data, put on 40km grid, and smooth surface  ****************"
mpiexec -n $N obj/pismr -bif init.nc -Mx 141 -My 141 -Mz 101 -Mbz 20 -gk -e 1.2 \
   -verbose -y 10 -o ant10yr_40km -of n
date

echo "  2.  running for 150k year to equilibriate temp (steady geometry)  ********"
mpiexec -n $N obj/pismr -if ant10yr_40km.nc -gk -e 1.2 -no_mass \
   -y 150000 -o ant150k_40km -of n
date

fi

echo "  3.  regrid to 20km grid: interpolate temp and age but get other fields ***"
echo "      from BAS data                                                      ***"
#mpiexec -n $N obj/pismr -bif init.nc -Mx 281 -My 281 -Mz 201 -Mbz 40 -gk -e 1.2 -no_mass\
#   -regrid ant150k_40km.nc -regrid_vars TBe -verbose -y 0.01 -o ant150k_20km -of n
mpiexec -n $N obj/pismr -if ant150k_40km.nc -gk -e 1.2 -y 0.0001 -no_mass -o ant150k_40km -of p
mpiexec -n $N obj/pismr -bif init.nc -Mx 281 -My 281 -Mz 201 -Mbz 40 -gk -e 1.2 -no_mass\
   -regrid ant150k_40km.pb -regrid_vars TBe -verbose -y 0.01 -o ant150k_20km -of n
date

echo "  4. smoothing surface on fine grid                             ************"
mpiexec -n $N obj/pismr -if ant150k_20km.nc -gk -e 1.2\
   -verbose -y 9.99 -o ant150k_20km_smooth -of n
date

echo "  5. running for another 3000 yrs to equilibriate temp on fine  ************"
echo "     grid (steady geometry)                                     ************"
mpiexec -n $N obj/pismr -if ant150k_20km_smooth.nc -gk -e 1.2 -no_mass\
   -y 2980 -o ant153k_20km -of n
date

echo "  6. this generates the drag map (which we could use to put the model  *****"
echo "     in steady state                                                   *****"
obj/get_drag -if ant153k_20km.nc -balvel init.nc -gk -d n -verbose -o ant_20km_drag -of n

echo "  6b. refining in the vertical and running for 10 yrs smoothing  ************"
echo "      to equilibriate temp on fine grid (steady geom)            ************"
mpiexec -n $N obj/pismr -bif init.nc -Mx 281 -My 281 -Mz 401 -Mbz 80 -gk -e 1.2 \
   -regrid ant153k_20km.nc -regrid_vars TBe -verbose -y 10.0 -o ant154k_20kmVFINE -of n
date

echo "  6c. running on fine for 100 yrs with no mass cons  ************"
mpiexec -n $N obj/pismr -if ant154k_20kmVFINE.nc -gk -e 1.2 -y 100 -no_mass \
   -o ant154k_20kmVFINE_p100yr -of n

echo "  7. run mass continuity and MacAyeal-Morland equations for   **************"
echo "     streams and shelves for one year w constant temps        **************"
mpiexec -n $N obj/pismr -if ant153k_20km.nc -gk -e 1.2 -mv -mv_rtol 1e-4 -ocean_kill\
   -no_temp -verbose -y 1 -o ant153k_mv1_20km -of n
date

echo "  8. complete 10 years w constant temps          ***************************"
mpiexec -n $N obj/pismr -if ant153k_mv1_20km.nc -gk -e 1.2 -mv -mv_rtol 1e-4 -ocean_kill\
   -no_temp -verbose -y 9 -o ant153k_mv10_20km -of n
date

echo "  9. complete 100 years w constant temps        ****************************"
mpiexec -n $N obj/pismr -if ant153k_mv10_20km.nc -gk -e 1.2 -mv -mv_rtol 1e-4 -ocean_kill\
   -no_temp -verbose -y 90 -o ant153k_mv100_20km -of n
date

echo "  10. complete 1000 years w constant temps     *****************************"
mpiexec -n $N obj/pismr -if ant153k_mv100_20km.nc -gk -e 1.2 -mv -mv_rtol 1e-4 -ocean_kill\
   -no_temp -verbose -y 900 -o ant153k_mv1k_20km -of n
date

echo "  11A. regrid SIA equilibriated from 20km onto 14km grid    *********************"
mpiexec -n $N obj/pismr -bif init.nc -Mx 401 -My 401 -Mz 51 -Mbz 10 -gk -e 1.2 \
   -regrid ant153k_20km.nc -regrid_vars TBe -verbose -y 0.001 -o ant153k_14kmTHIN -of n
date

echo "  11B. smooth for 10 years on 14km grid             *****************************"
mpiexec -n $N obj/pismr -if ant153k_14kmTHIN.nc -gk -e 1.2 -y 10 -no_temp -verbose \
   -o ant153k_p10yr_14kmTHIN -of n
date

echo "  11C. compute MacAyeal for 1 year on 14km grid     *****************************"
mpiexec -n $N obj/pismr -if ant153k_p10yr_14kmTHIN.nc -gk -e 1.2 -mv -mv_rtol 1e-4 \
     -ocean_kill -no_temp -y 1 -verbose -o ant153k_p10yr_mv1yr_14kmTHIN -of n
date
