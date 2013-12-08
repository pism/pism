#!/bin/bash

# Copyright (C) 2009-2013 Ed Bueler and Andy Aschwanden

#  creates 18 scripts, each with NN processors, for a parameter study
#  scripts are suitable for PBS job scheduler
#  (see  http://www.adaptivecomputing.com/products/open-source/torque/)
#
#  usage: to use NN=8 processors, 2 4-core nodes, and duration 4:00:00,
#     $ export PISM_WALLTIME=4:00:00
#     $ export PISM_NODES=2
#     $ ./paramspawn.sh 8
#  then, assuming you like the resulting scripts:
#     $ ./paramsubmit.sh      ### <--- REALLY SUBMITS using qsub


set -e # exit on error
SCRIPTNAME=paramspawn.sh

NN=8  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "paramspawn.sh 8" then NN = 8
  NN="$1"
fi

# set wallclock time
if [ -n "${PISM_WALLTIME:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME                    PISM_WALLTIME = $PISM_WALLTIME  (already set)"
else
  PISM_WALLTIME=4:00:00
  echo "$SCRIPTNAME                     PISM_WALLTIME = $PISM_WALLTIME"
fi
WALLTIME=$PISM_WALLTIME

# set number of nodes
if [ -n "${PISM_NODES:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME                    PISM_NODES = $PISM_NODES  (already set)"
else
  PISM_NODES=2
  echo "$SCRIPTNAME                     PISM_NODES = $PISM_NODES"
fi
NODES=$PISM_NODES

 SHEBANGLINE="#!/bin/bash"
MPIQUEUELINE="#PBS -q standard_4"
 MPITIMELINE="#PBS -l walltime=$WALLTIME"
 MPISIZELINE="#PBS -l nodes=$NODES:ppn=4"
  MPIOUTLINE="#PBS -j oe"

GRID=20
DURA=10000
for PPQ in 0.1 0.25 0.8 ; do
  for TEFO in 0.01 0.02 0.05 ; do

      SCRIPT="do_ppq_${PPQ}_tefo_${TEFO}_sgl_true.sh"
      rm -f $SCRIPT
      EXPERIMENT=ppq_${PPQ}_tefo_${TEFO}_sgl_true
      OUTFILE=g${GRID}km_${PPQ}_${TEFO}.nc

      # insert preamble
      echo $SHEBANGLINE >> $SCRIPT
      echo >> $SCRIPT # add newline
      echo $MPIQUEUELINE >> $SCRIPT
      echo $MPITIMELINE >> $SCRIPT
      echo $MPISIZELINE >> $SCRIPT
      echo $MPIOUTLINE >> $SCRIPT
      echo >> $SCRIPT # add newline
      echo "cd \$PBS_O_WORKDIR" >> $SCRIPT
      echo >> $SCRIPT # add newline

      export PISM_EXPERIMENT=$EXPERIMENT
      export PISM_TITLE="Greenland Parameter Study"

      cmd="PISM_DO="" PARAM_PPQ=$PPQ PARAM_TEFO=$TEFO ./spinup.sh $NN const $DURA $GRID hybrid $OUTFILE"
      echo "$cmd 2>&1 | tee job.\${PBS_JOBID}" >> $SCRIPT

      echo "($SPAWNSCRIPT)  $SCRIPT written"

      SCRIPT="do_ppq_${PPQ}_tefo_${TEFO}_sgl_false.sh"
      rm -f $SCRIPT
      EXPERIMENT=ppq_${PPQ}_tefo_${TEFO}_sgl_false
      OUTFILE=g${GRID}km_${PPQ}_${TEFO}_NOSGL.nc

      # insert preamble
      echo $SHEBANGLINE >> $SCRIPT
      echo >> $SCRIPT # add newline
      echo $MPIQUEUELINE >> $SCRIPT
      echo $MPITIMELINE >> $SCRIPT
      echo $MPISIZELINE >> $SCRIPT
      echo $MPIOUTLINE >> $SCRIPT
      echo >> $SCRIPT # add newline
      echo "cd \$PBS_O_WORKDIR" >> $SCRIPT
      echo >> $SCRIPT # add newline

      export PISM_EXPERIMENT=$EXPERIMENT
      export PISM_TITLE="Greenland Parameter Study"

      cmd="PISM_DO="" PARAM_PPQ=$PPQ PARAM_TEFO=$TEFO PARAM_NOSGL=foo ./spinup.sh $NN const $DURA $GRID hybrid $OUTFILE"
      echo "$cmd 2>&1 | tee job.\${PBS_JOBID}" >> $SCRIPT

      echo "($SPAWNSCRIPT)  $SCRIPT written"

  done
done


echo
echo "($SPAWNSCRIPT)  use paramsubmit.sh to submit the scripts"
echo

