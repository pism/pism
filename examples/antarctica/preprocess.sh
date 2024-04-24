#!/bin/bash
# Copyright (C) 2009-2018, 2020, 2023, 2024 PISM authors
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
*fill = -9999f;

*Nx = $x1.size;
*Ny = $y1.size;

// Define the new grid:
defdim("x", Nx);
defdim("y", Nx);

*dx = x1(1) - x1(0);
x = array(x1(0), dx, $x);
x@units = "meter";

*dy = y1(1) - y1(0);
y = array(y1(0), dy, $y);
y@units = "meter";

// Ice thickness:
thk[$y,$x] = 0.0f;
thk(0:Ny-1,0:Nx-1)=usrf-lsrf;
thk@units="meter";
thk@long_name="ice thickness";
thk@standard_name="land_ice_thickness";

// choose Van de Berg et al version of accumulation; will treat as
// ice-equivalent snow rate and convert from an ice-equivalent
// thickness rate ("m year-1") to "kg m-2 year-1" by assuming ice
// density of 918 kg m-3 (see the reference)
*ice_density = 918.0;
precipitation[$y, $x]=0.0f;
precipitation(0:Ny-1,0:Nx-1) = accr;
where(precipitation > fill)
  precipitation = precipitation * ice_density;
elsewhere
  precipitation = 0;
precipitation@units="kg m-2 year-1";
precipitation@long_name = "precipitation";

air_temp[$y, $x] = 0.0f;
air_temp(0:Ny-1,0:Nx-1) = temp;
air_temp@units="degree_Celsius";
air_temp@long_name = temp@long_name;

bed[$y, $x] = fill;
bed(0:Ny-1,0:Nx-1) = topg;
bed@units = "meter";
bed@standard_name = "bedrock_altitude";
where(bed == fill) bed = -5300;

// Use the Shapiro & Ritzwoller geothermal flux
bheatflx[$y, $x] = 60f;
bheatflx(0:Ny-1,0:Nx-1) = ghfsr;
bheatflx@units = "mW m-2";
bheatflx@long_name = ghfsr@long_name;

land_ice_area_fraction_retreat[$y, $x]=0;
where(thk > 0 || bed > 0) land_ice_area_fraction_retreat=1;
land_ice_area_fraction_retreat@long_name="maximum ice extent mask";
land_ice_area_fraction_retreat@units="1";

// The ALBMAPv1 paper states that the dataset uses the EIGEN-GL04C geoid, but unfortunately
// it is not easy to use. Hopefully it is okay to use WGS84 instead.
global@proj = "EPSG:3031";
'

ncap2 -O --script "$script" $DATANAME $PISMVERSION

# remove unneeded variables
ncks -O -v thk,bed,precipitation,air_temp,bheatflx,land_ice_area_fraction_retreat $PISMVERSION $PISMVERSION

echo "  PISM-readable file $PISMVERSION created; only has fields"
echo "    used in bootstrapping."

echo "now run spin-up script 'antspin-coarse.sh'"
echo
