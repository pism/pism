#!/bin/bash

# Copyright (C) 2009-2013 Ed Bueler and Andy Aschwanden

# submits scripts produced by paramspawn.sh; uses QSUB environment variable if set
# "qsub" is from PBS job scheduler
# (see  http://www.adaptivecomputing.com/products/open-source/torque/)
#
# usage for real, using qsub:
#   $ ./paramsubmit.sh
#
# usage for test:
#   $ PISM_QSUB=cat ./paramsubmit.sh

set -e -x # exit on error

SCRIPTNAME=paramsubmit.sh

# submission command
if [ -n "${PISM_QSUB:+1}" ] ; then  # check if env var PREFIX is already set
    QSUB=$PISM_QSUB
    echo "($SCRIPTNAME) QSUB = $PISM_QSUB"
else
    QSUB="qsub"
    echo "($SCRIPTNAME) QSUB = $QSUB"
fi

for SCRIPT in do_*.sh
do
  echo "($SCRIPTNAME) doing '$QSUB $SCRIPT' ..."
  $QSUB $SCRIPT
done
