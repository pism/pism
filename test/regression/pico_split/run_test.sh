#!/bin/bash

# Copyright (C) 2021 PISM authors
# created by torsten.albrecht@pik-potsdam.de

set -e
set -u
set -x

###############################################################################
# Test simulation of Antarctic ice sheet model to check on PICO adjustments
# concerning large ice shelves crossing basin boundaries.
###############################################################################

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

input_file=${PISM_SOURCE_DIR}/test/regression/pico_split/bedmap2_schmidtko14_50km.nc

temp_dir=$(mktemp -d -t pism_pico_split.XXXX) || exit 1

cp $input_file ${temp_dir}/input_01.nc
cp $input_file ${temp_dir}

pushd ${temp_dir}

cat > script.txt <<EOF
mask=basins * 0;
mask@units="";
mask@flag_values="0, 2, 3, 4";
mask@long_name="ice-type (ice-free/grounded/floating/ocean) integer mask";
mask@pism_intent="diagnostic";

where(thk >  topg / (-910.0 / 1028.0)) mask=2;
where(thk <= topg / (-910.0 / 1028.0)) mask=3;
where(thk == 0 && topg < 0) mask=4;
EOF

ncap2 -S script.txt -O input_01.nc input_01.nc

grid="-bootstrap -Mx 120 -My 120 -Lz 6000 -Lbz 2000 -Mz 81 -Mbz 21"

stressbalance="-pik -stress_balance ssa+sia -ssa_method fd"
pico="-ocean pico -ocean_pico_file $input_file -gamma_T 1.0e-5 -overturning_coeff 0.8e6 -exclude_icerises -continental_shelf_depth -2500"
surface="-atmosphere uniform -surface simple"   #"-surface pik"

# Merge Amundsen Sea and Ross ice shelves
ncap2 -O -s "where((basins==12 || basins==14) && mask==2 && topg<-600 && lon<-90) thk=topg/(-910.0/1028.0)-10.0" input_01.nc input_02.nc

# Advance Ross Ice Shelf into Amundsen Sea basin, but without connection
ncap2 -O -s "where((basins==12 || basins==14) && mask==2 && topg<-600 && lon<-90 && lat<-77.4) thk=topg/(-910.0/1028.0)-10.0" input_01.nc input_03.nc

extract() {
# Extracts pico_temperature_box0 at a given location.
filename=$1
X=$2
Y=$3

# - print pico_temperature_box0
# - print the line containing "pico_temperature_box0" and the line after it
# - remove "box0" because it contains 0
# - remove everything *except* for digits and periods
ncks -H -d x,${X} -d y,${Y} -v pico_temperature_box0 ${filename} | \
  grep -A 1 'pico_temperature_box0 =' | \
  sed s/box0// | \
  tr -d -c "[:digit:]."
}

compare() {
# Compares two strings and fails the test if they are different.
basin=$1
A=$2
B=$3

if [ "$A" == "$B" ]
then
  echo Basin ${basin}: OK
else
  echo Basin ${basin}: FAILED
  exit 1
fi
}

for testcase in "01" "02" "03";
do
  echo Test $testcase:

  ${PISM_PATH}/pismr -verbose 2 -i input_${testcase}.nc \
              -config ${PISM_PATH}/pism_config.nc \
              $grid \
              $stressbalance \
              $surface \
              $pico \
              -y 0.001 \
              -o o_${testcase}.nc | tee output_${testcase}.log

  # assert that PICO input temperatures in basins are:
  T12a_ref=271.553311083714
  T12b_ref=271.541340785807
  T13_ref=272.781500154734
  T14_ref=273.596660429239

  T12=$(extract o_${testcase}.nc 59 37)
  T13=$(extract o_${testcase}.nc 41 34)
  T14=$(extract o_${testcase}.nc 28 49)

  # basin 12
  if [ ${testcase} == "01" ]
  then
    compare 12 $T12a_ref $T12
  else
    compare 12 $T12b_ref $T12
  fi

  # basin 13
  compare 13 $T13_ref $T13

  # basin 14
  compare 14 $T13_ref $T13
done

rm -rf ${temp_dir}
