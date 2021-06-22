#!/bin/bash

# Copyright (C) 2021 PISM authors
# created by torsten.albrecht@pik-potsdam.de

###############################################################################
# Test simulation of Antarctic ice sheet model to check on PICO adjustments 
# concerning large ice shelves crossing basin boundaries.
###############################################################################

INNAME=bedmap2_schmidtko14_50km.nc

# means in 19 basins
tlist=( 0.0 -1.7552066 -1.644025 -1.6328654 -1.5730474 -1.5132295 -1.7245011 -1.6822174 -0.6628277 -1.6146358 -1.3121665 -1.6178559 -1.6178559 -0.36849985 0.44666043 1.0385389 1.1641185 0.23013027 -1.1589965 -1.7181017 ) 
slist=( 34.5 34.65044 34.529594 34.312176 34.319885 34.32759 34.531147 34.484543 34.567535 34.576805 34.668884 34.648018 34.648018 34.41013 34.5505 34.69003 34.66446 34.53287 34.59139 34.664528)


# create mean salinity and ocean temperature fields
cp $INNAME input_test01.nc
ncap2 -O -s "mask=0*basins" input_test01.nc input_test01.nc
ncatted -O -a units,mask,o,c,"" \
           -a flag_values,mask,o,c,'0, 2, 3, 4' \
           -a long_name,salinity_ocean,d,, \
           -a pism_intent,mask,o,c,'diagnostic' \
           -a long_name,mask,o,c,'ice-type (ice-free/grounded/floating/ocean) integer mask' \
           -a flag_meanings,mask,o,c,"ice_free_bedrock grounded_ice floating_ice ice_free_ocean" \
           input_test01.nc
           ncap2 -O -s "where(thk>topg/(-910.0/1028.0)) mask=2;where(thk<=topg/(-910.0/1028.0)) mask=3;where(thk==0 && topg<0) mask=4" input_test01.nc input_test01.nc

ncatted -O -a history,global,d,, \
           -a history_of_appended_files,global,d,, \
           input_test01.nc


#SIXTEENKMGRID="-Mx 381 -My 381 -Lz 6000 -Lbz 2000 -Mz 81 -Mbz 21"
FIFTYKMGRID="-Mx 120 -My 120 -Lz 6000 -Lbz 2000 -Mz 81 -Mbz 21"

stressbalance="-pik -stress_balance ssa+sia -ssa_method fd"
pico="-ocean pico -ocean_pico_file $INNAME -gamma_T 1.0e-5 -overturning_coeff 0.8e6 -exclude_icerises -continental_shelf_depth -2500"
surface="-atmosphere uniform -surface simple"   #"-surface pik"
regrid_opts="-bootstrap $FIFTYKMGRID" # -regrid_file $INNAME"


# merge Amundsen Sea and Ross ice shelves
ncap2 -O -s "where((basins==12 || basins==14) && mask==2 && topg<-600 && lon<-90) thk=topg/(-910.0/1028.0)-10.0" input_test01.nc input_test02.nc

# advance Ross Ice Shelf into Amundsen Sea basin, but without connection
ncap2 -O -s "where((basins==12 || basins==14) && mask==2 && topg<-600 && lon<-90 && lat<-77.4) thk=topg/(-910.0/1028.0)-10.0" input_test01.nc input_test03.nc

# benchmarks
Tc12a=271.553311083714
Tc12b=271.541340785807
Tc13=272.781500154734
Tc14=273.596660429239

for testcase in "test01" "test02" "test03"; do
  echo
  echo $testcase:
  ../../bin/pismr -i input_${testcase}.nc $regrid_opts $stressbalance $surface $pico -y 0.001 -o ${testcase}.nc >> output_${testcase}.out

  # assert that pico input temperatures in basins 12, 13 and 14 are 271.553311083714, 272.781500154734 and 273.596660429239 K
  T12=`ncks -H -d x,59 -d y,37 -v pico_temperature_box0 ${testcase}.nc | grep -C 1 'pico_temperature_box0 =' | grep -v '^$' | awk 'END{print}' | sed 's/.\;//'`
  T13=`ncks -H -d x,41 -d y,34 -v pico_temperature_box0 ${testcase}.nc | grep -C 1 'pico_temperature_box0 =' | grep -v '^$' | awk 'END{print}' | sed 's/.\;//'`  
  T14=`ncks -H -d x,28 -d y,49 -v pico_temperature_box0 ${testcase}.nc | grep -C 1 'pico_temperature_box0 =' | grep -v '^$' | awk 'END{print}' | sed 's/.\;//'`

  T12f=$(bc -l <<<"${T12}")
  T13f=$(bc -l <<<"${T13}")
  T14f=$(bc -l <<<"${T14}")

  echo $T12f,$T13f,$T14f
  if [ ${testcase} == "test01" ] 
  then
     if [ "$T12f" == "$Tc12a" ] 
     then 
         echo basin 12 equal
     else
       echo basin 12 not-equal
       exit 1
     fi
  else
     if [ "$T12f" == "$Tc12b" ]
     then
         echo basin 12 equal
     else
       echo basin 12 not-equal
       exit 1
     fi
  fi
  if [ "$T13f" == "$Tc13" ]
     then
         echo basin 13 equal
  else
       echo basin 13 not-equal
       exit 1
  fi
  if [ "$T14f" == "$Tc14" ]
     then
         echo basin 14 equal
  else
       echo basin 14 not-equal
       exit 1
  fi
done

echo $?
rm -f output_* test0* input_test0* ; exit 0


