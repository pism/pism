#!/bin/bash

# Copyright (C) 2009-2011 Ed Bueler and Andy Aschwanden

#  submits scripts produced by spawnparam.sh; uses QSUB environment variable if set

# usage for real, using qsub:
#   $ ./submitparam.sh

# usage for test:
#   $ export PISM_QSUB=cat
#   $ ./submitparam.sh

set -e -x # exit on error

SCRIPTNAME=submitparam.sh

# submission command
if [ -n "${PISM_QSUB:+1}" ] ; then  # check if env var PREFIX is already set
    QSUB=$PISM_QSUB
    echo "($SCRIPTNAME) QSUB = $PISM_QSUB"
else
    QSUB="qsub"
    echo "($SCRIPTNAME) QSUB = $QSUB"
fi

for SCRIPT in do_*_*_*.sh
do
  echo "($SCRIPTNAME) doing '$QSUB $SCRIPT' ..."
  $QSUB $SCRIPT
done
