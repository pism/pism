#!/bin/bash

# Copyright (C) 2021 PISM authors
# created by torsten.albrecht@pik-potsdam.de

###############################################################################
# Test simulation of Antarctic ice sheet model to check on PICO adjustments 
# concerning large ice shelves crossing basin boundaries.
###############################################################################


INNAME=bedmap2_schmidtko14_50km.nc

# adding SMB
cp $INNAME input_test01.nc
ncap2 -O -s "climatic_mass_balance=0*salinity_ocean" input_test01.nc input_test01.nc
ncatted -O -a units,climatic_mass_balance,o,c,"kg m-2 year-1" \
           -a long_name,climatic_mass_balance,o,c,"surface mass balance (accumulation/ablation) rate" \
           -a standard_name,climatic_mass_balance,o,c,"land_ice_surface_specific_mass_balance_flux" \
           input_test01.nc


#SIXTEENKMGRID="-Mx 381 -My 381 -Lz 6000 -Lbz 2000 -Mz 81 -Mbz 21"
FIFTYKMGRID="-Mx 120 -My 120 -Lz 6000 -Lbz 2000 -Mz 81 -Mbz 21"

stressbalance="-pik -stress_balance ssa+sia -ssa_method fd"
pico="-ocean pico -ocean_pico_file $INNAME -gamma_T 1.0e-5 -overturning_coeff 0.8e6 -exclude_icerises -continental_shelf_depth -2500"
surface="-surface pik"
regrid_opts="-bootstrap $FIFTYKMGRID" # -regrid_file $INNAME"


# merge Amundsen Sea and Ross ice shelves
ncap2 -O -s "where((basins==12 || basins==14) && mask==2 && topg<-600 && lon<-90) thk=topg/(-910.0/1028.0)-10.0" input_test01.nc input_test02.nc

# advance Ross Ice Shelf into Amundsen Sea basin, but without connection
ncap2 -O -s "where((basins==12 || basins==14) && mask==2 && topg<-600 && lon<-90 && lat<-77.4) thk=topg/(-910.0/1028.0)-10.0" input_test01.nc input_test03.nc

# merge Ross Ie Shelf with neighboring basin to check on weighted average
ncap2 -O -s "where((basins==12 || basins==13) && mask==2 && lon<-90) topg=-600.0" input_test01.nc input_test04.nc
ncap2 -O -s "where((basins==12 || basins==13) && mask==2 && topg<=-600 && lon<-90) thk=topg/(-910.0/1028.0)-10.0" input_test04.nc input_test04.nc


for testcase in "test01" "test02" "test03" "test04"; do
  echo $testcase
  ../../bin/pismr -i input_${testcase}.nc $regrid_opts $stressbalance $surface $pico -y 0.001 -o ${testcase}.nc

  # assert that pico inout temperatures in basins 12, 13 and 14 are 271.553311083714, 272.781500154734 and 273.596660429239 K
  ncks -H -d x,59 -d y,37 -v pico_temperature_box0 ${testcase}.nc | grep -C 1 'pico_temperature_box0 =' | grep -v '^$' | awk 'END{print}'
  ncks -H -d x,41 -d y,34 -v pico_temperature_box0 ${testcase}.nc | grep -C 1 'pico_temperature_box0 =' | grep -v '^$' | awk 'END{print}'  
  ncks -H -d x,28 -d y,49 -v pico_temperature_box0 ${testcase}.nc | grep -C 1 'pico_temperature_box0 =' | grep -v '^$' | awk 'END{print}'
  # in test04 the values in basin 12 and 13 should be averaged to 271.673605540195 K

done



