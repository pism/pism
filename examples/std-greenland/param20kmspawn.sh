#!/bin/bash

# Copyright (C) 2009-2014 Ed Bueler and Andy Aschwanden

#  creates 18 scripts each with NN processors and (potentially) submits them
#    on pacman.arsc.edu

#  needs rawparamscript, trimparam.sh

#  usage: to use NN=8 processors, 2 4-core nodes, and duration 4:00:00,
#     $ export PISM_WALLTIME=4:00:00
#     $ export PISM_NODES=2
#     $ ./spawnparam.sh 
#     (assuming you like the resulting scripts)
#     $ ./submitparam.sh      ### <--- REALLY SUBMITS using qsub

#  see param20kmsubmit.sh


set -e # exit on error
SPAWNSCRIPT=param20kmspwan.sh


NN=8  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "spawnparam.sh 8" then NN = 8
  NN="$1"
fi

# set wallclocktime
if [ -n "${PISM_WALLTIME:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME                    PISM_WALLTIME = $PISM_WALLTIME  (already set)"
else
  PISM_WALLTIME=4:00:00
  echo "$SCRIPTNAME                     PISM_WALLTIME = $PISM_WALLTIME"
fi
WALLTIME=$PISM_WALLTIME

# set no of nodes
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


DURA=10000
for PPQ in 0.1 0.25 0.8 ; do
  for TEFO in 0.01 0.02 0.05 ; do

      SCRIPT="do_ppq_${PPQ}_tefo_${TEFO}_sgl_true.sh"
      rm -f $SCRIPT
      EXPERIMENT=ppq_${PPQ}_tefo_${TEFO}_sgl_true
      OUTFILE=g20km_${PPQ}_${TEFO}.nc

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
      
      PISM_DO="echo" PARAM_PPQ=$PPQ PARAM_TEFO=$TEFO  ./spinup.sh $NN const $DURA 20 hybrid $OUTFILE >> $SCRIPT
	  
      echo "($SPAWNSCRIPT)  $SCRIPT written"

      SCRIPT="do_ppq_${PPQ}_tefo_${TEFO}_sgl_false.sh"
      rm -f $SCRIPT
      EXPERIMENT=ppq_${PPQ}_tefo_${TEFO}_sgl_false
      OUTFILE=g20km_${PPQ}_${TEFO}_NOSGL.nc

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
      
      cmd="PISM_DO="echo" PARAM_PPQ=$PPQ PARAM_TEFO=$TEFO PARAM_NOSGL=foo ./spinup.sh $NN const $DURA 20 hybrid $OUTFILE"
      echo "$cmd 2>&1 | tee job.${PBS_JOBID}" >> $SCRIPT
	  
      echo "($SPAWNSCRIPT)  $SCRIPT written"

  done
done


echo
echo "($SPAWNSCRIPT)  use param20kmsubmit.sh to submit the scripts"
echo

