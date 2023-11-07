#!/bin/bash
# Copyright (C) 2009-2018, 2020, 2023 PISM authors
#
# downloads the ALBMAPv1 dataset and prepares it for PISM
#
# depends on wget, unzip, NCO (ncap2, ncks)
#
# Reference:
#
# A. M. L. Brocq, A. J. Payne, and A. Vieli, “An improved Antarctic dataset for high
# resolution numerical ice sheet models (ALBMAP v1),” Earth System Science Data, vol. 2,
# no. 2, pp. 247–260, Oct. 2010, doi: 10.5194/essd-2-247-2010.
#
# Dataset DOI: 10.1594/PANGAEA.734145

# exit on error
set -e
# stop if a variable is not defined
set -u

DATANAME=ALBMAPv1.nc

echo
echo "downloading '$DATANAME'..."
wget -nc https://store.pangaea.de/Publications/LeBrocq_et_al_2010/$DATANAME.zip
unzip -o $DATANAME.zip

echo "making a PISM-readable file by copying parts of $DATANAME"
echo "  and adjusting metadata ..."

PISMVERSION=pism_Antarctica_5km.nc

script='
thk=usrf-lsrf;
thk@units="meter";
thk@long_name="ice thickness";
thk@standard_name="land_ice_thickness";

// choose Van de Berg et al version of accumulation; will treat as
// ice-equivalent snow rate and convert from an ice-equivalent
// thickness rate ("m year-1") to "kg m-2 year-1" by assuming ice
// density of 918 kg m-3 (see the reference)
precipitation[$y1, $x1]=0.0;
where(accr > -9999) precipitation=accr*918.0;
precipitation@units="kg m-2 year-1";
precipitation@long_name = "precipitation";

air_temp=temp;
air_temp@units="Celsius";

where(topg == -9999) topg = -5300;

// Use the Shapiro & Ritzwoller geothermal flux
bheatflx=ghfsr;
bheatflx@units="mW m-2";

land_ice_area_fraction_retreat[$y1, $x1]=0;
where(thk > 0 || topg > 0) land_ice_area_fraction_retreat=1;
land_ice_area_fraction_retreat@long_name="maximum ice extent mask";
land_ice_area_fraction_retreat@units="1";

// The ALBMAPv1 paper states that the dataset uses the EIGEN-GL04C geoid, but unfortunately
// it is not easy to use. Hopefully it is okay to use WGS84 instead.
global@proj = "EPSG:3031";
'

ncap2 -O --script "$script" $DATANAME $PISMVERSION
# remove unneeded variables
ncks -O -v thk,topg,precipitation,air_temp,bheatflx,land_ice_area_fraction_retreat $PISMVERSION $PISMVERSION

echo "  PISM-readable file $PISMVERSION created; only has fields"
echo "    used in bootstrapping."

echo "now run spin-up script 'antspin-coarse.sh'"
echo
