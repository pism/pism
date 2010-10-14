#!/bin/bash

# Copyright (C) 2010 Ed Bueler

# uses the PISM SeaRISE Greenland example to illustrate the use of regional
# climate model (RCM) output to find PDD parameters which produce closer fit of
# surface mass balance from PISM's PDD model to the RCM model output

# depends on all tools in examples/searise-greenland

set -e  # exit on error

echo
echo "# PDDGREEN example: preprocessing"
echo
# FIXME: could do  echo "# running ../searise-greenland/preprocess.sh and making links to result .."

# get file; see page http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
DATAVERSION=1.1
DATANAME=Greenland_5km_v$DATAVERSION.nc
PISMVERSION=pism_$DATANAME
wget -nc http://websrv.cs.umt.edu/isis/images/a/a5/$DATANAME

ncks -O $DATANAME $PISMVERSION  # just copies over, but preserves history and global attrs

# adjust metadata; uses NCO (http://nco.sourceforge.net/)
ncrename -O -v x1,x -v y1,y -d x1,x -d y1,y $PISMVERSION
ncrename -O -v time,t -d time,t $PISMVERSION
ncrename -O -v usrf,usurf $PISMVERSION

# we will use present surface temps from Fausto et al 2009 (J. Glaciol. vol 55 no 189)
#   so we don't mess with the Ettema-provided air temps

# convert from water equiv to ice thickness change rate; assumes ice density 910.0 kg m-3
ncap -O -s "snowprecip=presprcp*(1000.0/910.0)" $PISMVERSION $PISMVERSION
ncatted -O -a units,snowprecip,a,c,"m a-1" $PISMVERSION
ncatted -O -a long_name,snowprecip,a,c,"ice-equivalent precipitation rate" $PISMVERSION

# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
ncatted -a standard_name,bheatflx,d,, $PISMVERSION

echo "  PISM-readable file $PISMVERSION created from $DATANAME"


