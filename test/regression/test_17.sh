#!/bin/bash

PISM_PATH=$1
MPIEXEC=$2

# Test name:
echo "Test #17: verif test G regression: thermo SIA w. time-dependent SMB."
# The list of files to delete when done.
files="test_17-G-out.txt verify.nc verify.nc~"

rm -f $files

# run test G
OPTS="-test G -Mbz 1 -Mz 31 -y 1000 -o_size small -verbose 1"
$PISM_PATH/pismv -Mx 31 -My 31 $OPTS   > test_17-G-out.txt
$PISM_PATH/pismv -Mx 41 -My 41 $OPTS  >> test_17-G-out.txt

# compare results
diff test_17-G-out.txt -  <<END-OF-OUTPUT
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.780389   32.430770    7.148950    0.016288
temp      :        maxT         avT    basemaxT     baseavT
               0.835718    0.250497    0.759924    0.153262
Sigma     :      maxSig       avSig
               8.570017    0.996037
surf vels :     maxUvec      avUvec        maxW         avW
               0.944836    0.200079    0.028356    0.004029
NUM ERRORS DONE
NUMERICAL ERRORS evaluated at final time (relative to exact solution):
geometry  :    prcntVOL        maxH         avH   relmaxETA
               0.782888   30.763934    7.242048    0.015469
temp      :        maxT         avT    basemaxT     baseavT
               0.895269    0.249036    0.754972    0.157030
Sigma     :      maxSig       avSig
               8.211503    0.945963
surf vels :     maxUvec      avUvec        maxW         avW
               0.885743    0.194988    0.027397    0.004220
NUM ERRORS DONE
END-OF-OUTPUT

if [ $? != 0 ];
then
    exit 1
fi

rm -f $files; exit 0

