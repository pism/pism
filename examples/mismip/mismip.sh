#!/bin/bash
 
# Script for MISMIP = Marine Ice Sheet Intercomparison Project.  See
#    http://homepages.ulb.ac.be/~fpattyn/mismip/
# and especially mismip_4.pdf (or latest version of intercomparison
# description) at that site.

# This intercomparison describes a flow line ice stream and ice shelf
# system.  PISM is a 3D model so there is a notable lack of efficiency
# relative to flow line models.  Furthermore there is (for now) no 
# boundary layer or other advanced treatment of dynamics in the vicinity 
# of the grounding line; compare C. Schoof, (2007a). "Marine ice-sheet
# dynamics. Part 1. The case of rapid sliding," J. Fluid Mech., vol. 573,
# pp. 27–55.  Also there is no calving front boundary condition; compare 
# EISMINT-Ross on this issue.  See also comments at start of
# src/ismip/iceMISMIPModel.cc regarding grounding line flux issues.

# There are two models: Model one ("ABC1_...") uses only the SSA.
# Model two ("ABC2_...") is an average of SSA and SIA velocities.

# A PISM MISMIP run always saves an ASCII file in the MISMIP specified format:
#   ABC1_1b_M1_A1_t
# and a standard PISM NetCDF file
#   ABC1_1b_M1_A1.nc
# If the run achieves steady state by the MISMIP criterion, then the run also
# saves two more ASCII files in the MISMIP spec format:
#   ABC1_1b_M1_A1_ss
#   ABC1_1b_M1_A1_f
# See intercomparison description (URL above) for format of these files.  Also
#   ABC1_1b_M1_A1_extras
# may be written.

# The initial runs (ABC?_*_M?_A1 for *=1a,1b,3a,3b) may start from the 
# (Schoof, 2007a)-computed profile (thickness).  These profiles are generated 
# by the solverSM.py script.

# See scripts figMISMIP.py, makefigs.sh, and showflux.py for creating figures
# from outputs.

# Usage:
#     $ ./mismip.sh $N          use $N processors
#     $ ./mismip.sh $N D        debug mode: show the generated pisms calls
#                               but don't do them
#     $ ./mismip.sh $N SM       run solverSM.py to generate initial thicknesses
#                               (SM_??_A?.nc) before running PISM
#     $ ./mismip.sh $N D SM     this *is* allowed


NN=2               # set default number of processors here

MPIDO=mpiexec      # change this if need "mpirun", etc.

MYINITIALS=ABC     # MISMIP says "first character of the first name followed
                   #    by the first two characters of the last name"

#MODELRANGE="1 2"   # model 1 is pure SSA, model 2 is SIA+SSA
MODELRANGE="2"

#GRIDRANGE="1 2 3" # grid 1 is "-My 151", grid 2 is "-My 1501", grid 3 is "-My 601"
                   # 3 is slow, 2 is practically impossible!!
GRIDRANGE="1 3"


if [ $# -gt 3 ] ; then
  echo "error: too many parameters to mismip.sh"
  exit
fi

if [ $# -gt 0 ] ; then  # if user says "mismip.sh 8" then NN = 8
  NN="$1"
fi

SHOWONLY=0
if [ "$2" = D -o "$3" = D ] ; then  # "mismip.sh 8 D" gives NN = 8 and
                                    # only shows calls; no run 
  SHOWONLY=1
fi

SOLVERFIRST=0
if [ "$2" = SM -o "$3" = SM ] ; then  # "mismip.sh 8 SM" gives NN = 8 and  
  SOLVERFIRST=1                       # runs solverSM.py
fi                                    

# "mismip.sh 8 D SM" makes SHOWONLY=1 and SOLVERFIRST=1

set -e  # exit on error


# the function which turns
#    mismip NP NUM EXPER OPTS
# into
#    "mpiexec -np NP pisms -model NUM -mismip EXPER -initials MYINITIALS 
#         -extras -ksp_rtol 1.0e-7 OPTS"
# on NN processors
mpimismip()
{
    # change this if "bin/pisms", etc.:
    cmd="$MPIDO -np $1 pisms -model $2 -mismip $3 -initials $MYINITIALS -extras -ksp_rtol 1.0e-7 $4"
    if [ $SHOWONLY = 1 ] ; then
      echo "would try '$cmd'"
    else
      echo "date = '`date`' on host '`uname -n`':"
      echo "trying '$cmd'"
      echo
      $cmd
    fi
}


### compute initial thicknesses from Schoof 2007 if desired
if [ $SOLVERFIRST = 1 ] ; then
  for GRID in $GRIDRANGE
  do
    case $GRID in
      1   ) GRIDMY=151;;
      2   ) GRIDMY=1501;;
      3   ) GRIDMY=601;;
    esac
      
    for EXPER in 1 3
    do
      for SLIDING in a b
      do
        name=SMthk_${EXPER}${SLIDING}_M${GRID}_A1
        cmd="./solverSM.py -m $GRIDMY -e $EXPER -l $SLIDING -s 1 -n $name.nc"
        echo "trying '$cmd'"
        $cmd
      done
    done

  done   # GRID
fi   # $SOLVERFIRST = 1

### if you want to instead generate .png for each solverSM.py result,
### instead of .nc above, use these commands inside loops above, but be sure
### to uncomment "exit" here so rest of script does not run
# name=SM_${EXPER}${SLIDING}_M${GRID}_A1
# cmd="./solverSM.py -m $GRIDMY -e $EXPER -l $SLIDING -s 1 -o $name.png"

#exit


### now do MISMIP runs
for MODEL in $MODELRANGE
do

  for GRID in $GRIDRANGE
  do
    case $GRID in
      1   ) GRIDMY=151; SKIP=;;
      2   ) GRIDMY=1501; SKIP="-skip 10";;
      3   ) GRIDMY=601; SKIP="-skip 10";;
    esac

    ## experiment 1:
    for ES in 1a 1b
    do
      smncname=SMthk_${ES}_M${GRID}_A1.nc
      regridstart="-regrid_file ${smncname} -regrid_vars H"
      ### replace regridstart by empty string if we want to start w 10 m thickness
      #regridstart=
      options="-step 1 ${SKIP} -Mx 3 -Lz 6000 -Mz 15 Lbz 0 -Mbz 1 -My ${GRIDMY} ${regridstart}"
      mpimismip $NN $MODEL $ES "$options"
      for STEP in 2 3 4 5 6 7 8 9
      do
        PREV=$(($STEP-1))
        INNAME=${MYINITIALS}${MODEL}_${ES}_M${GRID}_A${PREV}.nc
        mpimismip $NN $MODEL $ES "-step ${STEP} ${SKIP} -i ${INNAME}"
      done # for STEP
    done # for ES

    # experiment 2 starts from end of 1, so above loop should be checked
    # for completion:
    for SLIDING in a b
    do
      INNAME=${MYINITIALS}${MODEL}_1${SLIDING}_M${GRID}_A9.nc
      mpimismip $NN $MODEL 2$SLIDING "-step 8 ${SKIP} -i ${INNAME}"
      for STEP in 7 6 5 4 3 2 1
      do
        PREV=$(($STEP+1))
        INNAME=${MYINITIALS}${MODEL}_2${SLIDING}_M${GRID}_A${PREV}.nc
        mpimismip $NN $MODEL 2$SLIDING "-step ${STEP} ${SKIP} -i ${INNAME}"
      done # for STEP
    done # for SLIDING

    # exper 3a
    smncname=SMthk_3a_M${GRID}_A1.nc
    regridstart="-regrid ${smncname} -regrid_vars H"
    ### replace regridstart by empty string if we want to start w 10 m thickness
    #regridstart=
    options="-step 1 ${SKIP}-Mx 3 -Mz 15 -My ${GRIDMY} ${regridstart}"
    mpimismip $NN $MODEL 3a "$options"
    for STEP in 2 3 4 5 6 7 8 9 10 11 12 13
    do
      PREV=$(($STEP-1))
      INNAME=${MYINITIALS}${MODEL}_3a_M${GRID}_A${PREV}.nc
      mpimismip $NN $MODEL 3a "-step ${STEP} ${SKIP} -i ${INNAME}"
    done # for STEP

    # exper 3b
    smncname=SMthk_3b_M${GRID}_A1.nc
    regridstart="-regrid ${smncname} -regrid_vars H"
    ### replace regridstart by empty string if we want to start w 10 m thickness
    #regridstart=
    options="-step 1 ${SKIP} -Mx 3 -Mz 15 -My ${GRIDMY} ${regridstart}"
    mpimismip $NN $MODEL 3b "$options"
    for STEP in 2 3 4 5 6 7 8 9 10 11 12 13 14 15
    do
      PREV=$(($STEP-1))
      INNAME=${MYINITIALS}${MODEL}_3b_M${GRID}_A${PREV}.nc
      mpimismip $NN $MODEL 3b "-step ${STEP} ${SKIP} -i ${INNAME}"
    done # for STEP

  done # for GRID

done # for MODEL

