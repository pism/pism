#!/bin/bash

# This script downloads MEaSUREs InSAR-Based Antarctica Ice Velocity Map, Version 2 from NSIDC
#
# It follows directions from
# https://nsidc.org/support/how/how-do-i-programmatically-request-data-services#mac
#
# 1. Create an account at https://urs.earthdata.nasa.gov/
# 2. Run this script as USER=your_username PASSWORD=your_password ./download.sh

set -u
set -e

URL=n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0484.002/1996.01.01/antarctica_ice_velocity_450m_v2.nc

wget \
  --auth-no-challenge=on \
  --execute robots=off \
  --http-password ${PASSWORD} \
  --http-user ${USER} \
  --keep-session-cookies \
  --load-cookies ~/.urs_cookies \
  --no-check-certificate \
  --no-clobber \
  --no-parent \
  --recursive \
  --reject "index.html*" \
  --save-cookies ~/.urs_cookies \
  https://${URL}

ln -s ${URL} .
