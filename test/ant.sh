# ANT.SH  is designed for modelling the current state of Antarctica on a 
# multi-processor machine.  Ice streams are determined by using balance velocities.
# It assumes  init.nc  is present (contains 5km Antarctica data mostly from 
# the British Antarctica Service; write to ffelb@uaf.edu).   It has been
# successfully run on marmaduke.gi.alaska.edu with N=4 processors, with a total
# run time of about 5 days.  ELB 1/31/07

N=4
echo "  ant.sh running with $N processors"

echo "  1.  input BAS data, put on 40km grid, and smooth surface  ****************"
mpiexec -n $N obj/pismr -if init.nc -Mx 141 -My 141 -Mz 101 -Mbz 20 -gk -e 1.2 \
   -verbose -y 10 -o ant10yr_40km
date

echo "  2.  running for 150k year to equilibriate temp (steady geometry)  ********"
mpiexec -n $N obj/pismr -if ant10yr_40km.pb -gk -e 1.2 -no_mass_bal \
   -y 150000 -o ant150k_40km
date

echo "  3.  regrid to 20km grid: interpolate temp and age but get other fields ***"
echo "      from BAS data; using ONLY ONE PROCESSOR FOR REGRID                 ***"
obj/pismr -if init.nc -Mx 281 -My 281 -Mz 201 -Mbz 40 -gk -e 1.2 -no_mass_bal\
   -regrid ant150k_40km.pb -regrid_vars TBe -verbose -y 0.01 -o ant150k_20km
date

echo "  4. smoothing surface on fine grid                             ************"
mpiexec -n $N obj/pismr -if ant150k_20km.pb -gk -e 1.2\
   -verbose -y 9.99 -o ant150k_20km_smooth
date

echo "  5. running for another 3000 yrs to equilibriate temp on fine  ************"
echo "     grid (steady geometry)                                     ************"
mpiexec -n $N obj/pismr -if ant150k_20km_smooth.pb -gk -e 1.2 -no_mass_bal\
   -y 2980 -o ant153k_20km
date

echo "  6. this generates the drag map (which we could use to put the model  *****"
echo "     in steady state                                                   *****"
obj/get_drag -if ant153k_20km.pb -balvel init.nc -gk -d n -verbose -o ant_20km_drag

# FOLLOWING WAS *NOT* DONE on marmaduke.gi.alaska.edu; PRODUCED OUT OF MEMORY MESSAGE
#echo "   refining in the vertical and running for another 1000 yrs  ************"
#echo "      to equilibriate temp on fine grid (steady geom)            ************"
#obj/pismr -if init.nc -Mx 281 -My 281 -Mz 401 -Mbz 80 -gk -e 1.2 -no_mass_bal\
#   -regrid ant153k_20km.pb -regrid_vars TBeHh -verbose -y 0.01 -o ant154k_20kmVFINE
#date

echo "  7. run mass continuity and MacAyeal-Morland equations for   **************"
echo "     streams and shelves for one year w constant temps        **************"
mpiexec -n $N obj/pismr -if ant153k_20km.pb -gk -e 1.2 -mv -mv_rtol 1e-4 -ocean_kill\
   -no_temp -verbose -y 1 -o ant153k_mv1_20km
date

echo "  8. complete 10 years w constant temps          ***************************"
mpiexec -n $N obj/pismr -if ant153k_mv1_20km.pb -gk -e 1.2 -mv -mv_rtol 1e-4 -ocean_kill\
   -no_temp -verbose -y 9 -o ant153k_mv10_20km
date

echo "  9. complete 100 years w constant temps        ****************************"
mpiexec -n $N obj/pismr -if ant153k_mv10_20km.pb -gk -e 1.2 -mv -mv_rtol 1e-4 -ocean_kill\
   -no_temp -verbose -y 90 -o ant153k_mv100_20km
date

echo "  10. complete 1000 years w constant temps     *****************************"
mpiexec -n $N obj/pismr -if ant153k_mv100_20km.pb -gk -e 1.2 -mv -mv_rtol 1e-4 -ocean_kill\
   -no_temp -verbose -y 900 -o ant153k_mv1k_20km
date
