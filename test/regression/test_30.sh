#!/bin/bash

# Test PISM's PDD model. Compare output from pclimate and pypdd.py.

# The test produces:
#
#  * pdd.nc:       pdd model configuration
#  * atm.nc:       artificial climate distribution;
#  * i.nc:         input for pclimate;
#  * pypdd_smb.nc: mass-balance from pypdd.py;
#  * pism_smb.nc:  mass-balance from pclimate.
#  * diff.nc:      difference between pypdd_smb.nc and pism_smb.nc
#
# The fake climate in atm.nc is defined by:
#
#     temp(t, x, y) = -10*y - 5*cos(2*pi*t)
#     prec(t, x, y) = x * (sign(x) - cos(2*pi*t))
#
# With x, y in km and t in yr. In other words:
#
#  * temperature amplitude is the same (5°C) everywhere
#  * temperature always peaks in summer
#  * mean temperature increases from North (-10°C) to South (10°C)
#  * min precipitation is zero everywhere
#  * mean precipitation is zero in the center line of the domain, and
#    increases to 1 m/yr in the western and eastern edges.
#  * precipitation peaks in winter in the western (maritime) domain half
#    and in summer in the eastern (continental) domain half.

PISM_PATH=$1
MPIEXEC=$2
PISM_SOURCE_DIR=$3

# List of files to remove when done:
files="i.nc atm.nc pypdd_smb.nc pism_smb.nc diff.nc"

rm -f $files

set -e

# Generate config overrides
echo "netcdf pdd {
variables:
	byte pism_overrides ;
		pism_overrides:pdd_factor_snow = 0.003 ;
		pism_overrides:pdd_factor_ice = 0.008 ;
		pism_overrides:pdd_refreeze = 0.6 ;
		pism_overrides:pdd_std_dev = 5 ;
		pism_overrides:pdd_balance_year_start_day = 274 ;
		pism_overrides:pdd_max_temperature_evals_per_year = 53 ;
		pism_overrides:air_temp_all_precip_as_snow = 273.15 ;
		pism_overrides:air_temp_all_precip_as_rain = 275.15 ;
}" | ncgen -o pdd.nc

# Run pypdd
./pypdd.py -b -o pypdd_smb.nc --interpolate-n 53

# Prepare an input file for pclimate
$PISM_PATH/pisms -Lx 1 -Ly 1 -Mx 201 -My 201 -y 0 -o i.nc

# Run pclimate
$PISM_PATH/pclimate -i i.nc -o pism_smb.nc -times 0:1:1 -atmosphere given -atmosphere_given_file atm.nc -atmosphere_given_period 1 -surface pdd -config_override pdd.nc -pdd_std_dev 5
ncks -O -d time,0 pism_smb.nc pism_smb.nc

set +e

# Check results:
for var in saccum smelt srunoff climatic_mass_balance; do
  $PISM_PATH/nccmp.py -v $var pypdd_smb.nc pism_smb.nc
done

# If test fails output difference
if [ $? != 0 ];
then
  ncdiff -O -v saccum,smelt,srunoff,climatic_mass_balance pism_smb.nc pypdd_smb.nc diff.nc
  exit 1
fi

rm -f $files; exit 0
