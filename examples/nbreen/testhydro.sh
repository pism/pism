#!/bin/bash

#set -v

MPIDO="mpiexec -n 4"
PATHTO=""

INPUTCHOICE="-hydrology_use_const_bmelt -hydrology_const_bmelt 1.5844e-08"  # 50 cm/year input, a lot

EXTRASNULL="-extra_times 0:0.02:1 -extra_vars thk,tillwat,hydroinput"
EXTRASROUT="-extra_times 0:0.02:1 -extra_vars thk,tillwat,hydroinput,bwat,bwp,bwatvel"

if [ "" ]; then      # zero length string is false
#if [ "pre" ]; then  # <-- do this to create "start.nc" for comparison runs
  $MPIDO ${PATHTO}pisms -Mx 121 -My 121 -Mz 101 -y 10 -o pre.nc
  $MPIDO ${PATHTO}pismr -i pre.nc -ye 8000 -o start.nc
fi

# only goal of this test script is to *LOOK* at the results; verification is another matter

# -hydrology null
$MPIDO ${PATHTO}pismr -i start.nc -hydrology null $INPUTCHOICE -ys 0 -y 1 $EXTRASNULL -extra_file extras_null.nc  -o null.nc

# finer grid -hydrology null
#$MPIDO ${PATHTO}pismr -boot_file start.nc -Mx 241 -My 241 -Lz 4100 -Mz 201 -hydrology null $INPUTCHOICE -ys 0 -y 1 -extra_file extras_finenull.nc $EXTRASNULL -o finenull.nc

# -hydrology routing
cp start.nc start_withbwat.nc
ncrename -v tillwat,bwat start_withbwat.nc
ncks -A -v tillwat start.nc start_withbwat.nc
$MPIDO ${PATHTO}pismr -i start_withbwat.nc -hydrology routing $INPUTCHOICE -report_mass_accounting -ys 0 -y 1 -extra_file extras_routing.nc $EXTRASROUT -o routing.nc

# finer grid -hydrology routing
#$MPIDO ${PATHTO}pismr -boot_file start_withbwat.nc -Mx 241 -My 241 -Lz 4100 -Mz 201 -hydrology routing $INPUTCHOICE -report_mass_accounting -ys 0 -y 1 -extra_file extras_finerouting.nc $EXTRASROUT -o finerouting.nc

# -hydrology distributed
# see for basic tests see
#    test/regression/test_29.py
#    examples/nbreen/run.sh

#set +v

