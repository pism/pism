#!/bin/bash

# Copyright (C) 2021 PISM authors
# created by torsten.albrecht@pik-potsdam.de

set -e
set -u
set -x

###############################################################################
# Test simulation of Antarctic ice sheet model to check on pik atmosphere model 
# in combination with elevation_change modifier
###############################################################################

#PISM_PATH=$1
PISM_PATH=/home/albrecht/pism21/pism-dev/bin
#MPIEXEC=$2
MPIEXEC=
#PISM_SOURCE_DIR=$3
PISM_SOURCE_DIR=/home/albrecht/pism21/pism-dev/

input_file=${PISM_SOURCE_DIR}/test/regression/pico_split/bedmap2_schmidtko14_50km.nc

#temp_dir=$(mktemp -d -t pism_pico_split.XXXX) || exit 1
temp_dir=tempdir
mkdir -p $temp_dir

cp $input_file ${temp_dir}/input_01.nc
cp $input_file ${temp_dir}

pushd ${temp_dir}

cat > script.txt <<EOF
usurf=topg+thk;
usurf@standard_name="surface_altitude";
usurf@long_name="ice upper surface elevation";
usurf@units="meters";
precipitation=salinity_ocean * 0.0 + 100.0;
precipitation@units="kg m-2 year-1";
precipitation@long_name="Total Precipitative Flux";
precipitation@standard_name="precipitation_flux";
mask=basins * 0;
mask@units="";
mask@flag_values="0, 2, 3, 4";
mask@long_name="ice-type (ice-free/grounded/floating/ocean) integer mask";
mask@pism_intent="diagnostic";

where(thk >  topg / (-910.0 / 1028.0)) mask=2;
where(thk <= topg / (-910.0 / 1028.0)) mask=3;
where(thk == 0 && topg < 0) mask=4;
where(mask==3) usurf=thk*(1.0-(910.0/1028.0));
usurf=usurf-100.0;
where(usurf < 0) usurf=0.0;

EOF



ncap2 -S script.txt -O input_01.nc input_01.nc

#ncap2 -O -s "usurf=usurf-100.0;where(usurf<0) usurf=0" input_01.nc input_02.nc


initfile=${temp_dir}/input_01.nc
initfile=input_01.nc
#initfile2=${PISM_SOURCE_DIR}/test/regression/pico_split/${temp_dir}/input_02.nc

grid="-bootstrap -Mx 120 -My 120 -Lz 6000 -Lbz 2000 -Mz 81 -Mbz 21"

stressbalance="-pik -stress_balance ssa+sia -ssa_method fd"
pico="-ocean pico -ocean_pico_file $input_file -gamma_T 1.0e-5 -overturning_coeff 0.8e6 -exclude_icerises -continental_shelf_depth -2500"
surface="-atmosphere pik,elevation_change -atmosphere_pik era_interim_lon -atmosphere_pik_file ${initfile} -atmosphere_lapse_rate_file ${initfile} -temp_lapse_rate 8.2 -precip_adjustment scale -surface pdd"
snaps="-save_times 0:100:1000 -save_split -save_size medium"


extract() {
# Extracts ice_surface_temp and climatic_mass_balance at a given location.
filename=$1
var=$2
X=$3
Y=$4

# - print var
# - print the line containing "var" and the line after it
# - remove everything *except* for digits and periods
ncks -H -d x,${X} -d y,${Y} -v ${var} ${filename} | \
  grep -A 1 "${var} =" | \
  tr -d -c "[:digit:]."
}

compare() {
# Compares two strings and fails the test if they are different.
A=$1
B=$2

if [ "$A" == "$B" ]
then
  echo OK
else
  echo FAILED
  exit 1
fi
}

#for testcase in "01" "02" "03";
for testcase in "01";
do
  echo Test $testcase:

  ${PISM_PATH}/pismr -verbose 2 -i input_${testcase}.nc \
              -config ${PISM_SOURCE_DIR}/share/pism/pism_config.nc \
              $grid \
              $stressbalance \
              $surface \
              $pico \
              $snaps -save_file s_${testcase} \
              -y 0.001 \
              -o o_${testcase}.nc | tee output_${testcase}.log

  # assert that temperatures and SMB at South Pole are:
  TSP_ref=223.7456
  PSP_ref=94.38937

  TSP=$(extract o_${testcase}.nc ice_surface_temp 60 60)
  PSP=$(extract o_${testcase}.nc climatic_mass_balance 60 60)
  echo $TSP,$PSP
  if [ ${testcase} == "01" ]
  then
    compare $TSP_ref $TSP
    compare $PSP_ref $PSP
  fi

done

rm -rf ${temp_dir}
