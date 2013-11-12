#!/bin/bash

# Copyright (C) 2009-2013 Ed Bueler and Andy Aschwanden

#  submits scripts produced by param20kmspawn.sh; uses QSUB environment variable if set

# usage for real, using qsub:
#   $ ./submitparam.sh

# usage for test:
#   $ export PISM_QSUB=cat
#   $ ./param20kmsubmit.sh

set -e -x # exit on error

SCRIPTNAME=param20kmsubmit.sh

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
