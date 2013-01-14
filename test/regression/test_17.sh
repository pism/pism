#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #17: verif test G regression: thermo SIA w. time-dependent SMB."
# The list of files to delete when done.
files="test_17-G-out.txt"

rm -f $files

# run test G
OPTS="-test G -Mbz 1 -Mz 31 -y 1000 -o_size none -verbose 1"
$PISM_PATH/pismv -Mx 31 -My 31 $OPTS   > test_17-G-out.txt
$PISM_PATH/pismv -Mx 41 -My 41 $OPTS  >> test_17-G-out.txt

# compare results
diff test_17-G-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.780389   32.423434    7.148981    0.016282
temp      :        maxT         avT    basemaxT     baseavT
               0.835721    0.250719    0.765509    0.153551
Sigma     :      maxSig       avSig
               8.715311    1.012199
surf vels :     maxUvec      avUvec        maxW         avW
               0.944746    0.200100    0.028355    0.004030
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.782888   30.761282    7.242044    0.015465
temp      :        maxT         avT    basemaxT     baseavT
               0.895267    0.249156    0.757865    0.157201
Sigma     :      maxSig       avSig
               8.294728    0.954785
surf vels :     maxUvec      avUvec        maxW         avW
               0.885667    0.194997    0.027396    0.004221
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

