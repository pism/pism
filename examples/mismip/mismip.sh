#!/bin/bash
 
# Script for MISMIP = Marine Ice Sheet Intercomparison Project.  See
#    http://homepages.ulb.ac.be/~fpattyn/mismip/
# and especially
#    http://homepages.ulb.ac.be/~fpattyn/mismip/mismip_4.pdf
# (or latest version of intercomparison description).

# This intercomparison describes a flow line ice stream and ice shelf
# system.  PISM is a 3D model so there is a notable lack of efficiency
# relative to flow line models.  Furthermore there is (for now) no 
# boundary layer or other advanced treatment of dynamics in the vicinity 
# of the grounding line; compare C. Schoof, (2007). "Marine ice-sheet
# dynamics. Part 1. The case of rapid sliding," J. Fluid Mech., vol. 573,
# pp. 27–55.  Also there is no calving front boundary condition; compare 
# EISMINT-Ross on this issue.  See also comments at start of
# src/ismip/iceMISMIPModel.cc regarding grounding line flux issues.

# There are two models: Model one ("ABC1_...") uses only the SSA.
# Model two ("ABC2_...") adds PISM option '-super', so that there is
# an average of SSA and SIA velocities.

# A PISM MISMIP run always saves an ASCII file in the MISMIP specified format:
#   ABC1_1b_M1_A1_t
# and a standard PISM NetCDF file
#   ABC1_1b_M1_A1.nc
# If the run achieves steady state by the MISMIP criterion, then the run also
# saves two more ASCII files in the MISMIP spec format:
#   ABC1_1b_M1_A1_ss
#   ABC1_1b_M1_A1_f
# See intercomparison description (URL above) for format of these files.


NN=2           # set default number of processors here

MPIDO=mpiexec  # change this if "mpirun", etc.

MYINITIALS=EBU # MISMIP says "first character of the first name followed
               #    by the first two characters of the last name"


if [ $# -gt 0 ] ; then  # if user says "mismip.sh 8" then NN = 8
  NN="$1"
fi

SHOWONLY=0
if [ $# -gt 1 ] ; then  # if user says "mismip.sh 8 D" then NN = 8 and only shows; no run 
  SHOWONLY=1
fi

set -e  # exit on error
shopt -s nullglob  # don't return list of files if none found

# function to run
#   "pisms -mismip EXPER -initials MYINITIALS -ksp_rtol 1e-11 OTHEROPTS"
#   on NN processors
mpimismip()
{
    # change this if "bin/pisms", etc.:
    cmd="$MPIDO -np $1 pisms -mismip $2 -initials $MYINITIALS -ksp_rtol 1e-11 -extras $3"
    
    echo 
    echo "date = '`date`' on host '`uname -n`':"
    echo "trying '$cmd'"
    echo
    if [ $SHOWONLY = 0 ] ; then
      $cmd
    fi
}


#uncomment these (and "#fi" at end) and move around to bypass completed stuff:
#if [ ]; then   # always goes to "else"
#else           # put this before restart location

## model 1 is pure SSA, model 2 is SIA+SSA

## exper 1:
#for MODEL in 1 2
for MODEL in 1 2
do

  #for GRID in 1 2 3  # 2 is slow!!
  for GRID in 1
  do

    case $GRID in
      1   ) GRIDMY=151;;
      2   ) GRIDMY=1501;;
      3   ) GRIDMY=601;;
    esac
    case $GRID in
      1   ) SKIP=;;
      2   ) SKIP="-skip 10 ";;
      3   ) SKIP="-skip 10 ";;
    esac

    for ES in 1a 1b
    do
      mpimismip $NN $ES "-model ${MODEL} -step 1 ${SKIP}-Mx 3 -Mz 11 -My ${GRIDMY}"
      for STEP in 2 3 4 5 6 7 8 9
      do
        PREV=$(($STEP-1))
        INNAME=${MYINITIALS}${MODEL}_${ES}_M${GRID}_A${PREV}.nc
        mpimismip $NN $ES "-model ${MODEL} -step ${STEP} ${SKIP}-if ${INNAME}"
      done # for STEP
    done # for ES

  done # for GRID

done # for MODEL

exit


