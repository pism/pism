#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

set -x

# Test name:
echo "Test #26: Cold part of Paterson-Budd flow law regression."
# The list of files to delete when done.
files="flowtable.txt diff.txt"

rm -f $files

$PISM_PATH/flowlaw_test -flow_law arr > flowtable.txt
# compare results
diff flowtable.txt - > diff.txt <<END-OF-OUTPUT
flow law:   "arr"
pressure = 1.785e+07 Pa = (hydrostatic at depth 2000.00 m)
flowtable:
  (dev stress)   (abs temp) (liq frac) =   (flow)
      1.00e+04      241.740      0.000 = 3.917295e-18
      1.00e+04      266.740      0.000 = 6.428034e-17
      1.00e+04      271.740      0.000 = 1.057468e-16
      1.00e+04      271.740      0.005 = 1.057468e-16
      5.00e+04      241.740      0.000 = 9.793238e-17
      5.00e+04      266.740      0.000 = 1.607008e-15
      5.00e+04      271.740      0.000 = 2.643671e-15
      5.00e+04      271.740      0.005 = 2.643671e-15
      1.00e+05      241.740      0.000 = 3.917295e-16
      1.00e+05      266.740      0.000 = 6.428034e-15
      1.00e+05      271.740      0.000 = 1.057468e-14
      1.00e+05      271.740      0.005 = 1.057468e-14
      1.50e+05      241.740      0.000 = 8.813914e-16
      1.50e+05      266.740      0.000 = 1.446308e-14
      1.50e+05      271.740      0.000 = 2.379304e-14
      1.50e+05      271.740      0.005 = 2.379304e-14
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

