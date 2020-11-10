#!/bin/bash

RUNDIR=/Users/flo/Recherche/MODELS/simulations/PISM/Ant_paleo_spinup_30km
OUTDIR=/Users/flo/Recherche/MODELS/WORK_PISM/ANT_paleo_spinup
DATAFILE=./



cd $RUNDIR
cp pism_config_ant_paleo.cdl $DATAFILE

# -------------------------------------------------------------------------
# Run the model
# -------------------------------------------------------------------------
./run_paleo.sh 4

#&> out.paleo
