#!/bin/bash

set -x

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #12: mass conservation within rounding error (SIA moving margin)."
# The list of files to delete when done.
files="verify-12.nc volume_nonglacierized-12.nc"

rm -f $files
GRID="-Mx 31 -My 31 -Mz 31 -y 5000 -ys 1000"
OPTS="-test B -max_dt 25 -o_size small"
TS_OPTS="-ts_file volume_nonglacierized-12.nc -ts_vars volume_nonglacierized -ts_times 1000:25:1e4"
# run test B 

set -x

$MPIEXEC -n 2 $PISM_PATH/pismv -test B $GRID $OPTS $TS_OPTS

set +x

/usr/bin/env python <<EOF
from numpy import diff, log10, floor
from sys import exit
try:
    from netCDF3 import Dataset
except:
    from netCDF4 import Dataset

nc = Dataset("volume_nonglacierized-12.nc", 'r')
volume = nc.variables['volume_nonglacierized'][:]
volume_max = volume.max()
threshold = 10**(floor(log10(volume_max)) - 14) # 14 digits of accuracy
diff_max = diff(volume).max()

if diff_max < threshold:
    print "diff(volume).max() = %f, threshold = %f" % (diff_max, threshold)
    exit(0)
else:
    print "diff(volume).max() = %f > %f" % (diff_max, threshold)
    exit(1)
EOF

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

