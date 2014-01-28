#!/bin/bash

if [ $# -lt 1 ] ; then
  echo "preprocess.sh ERROR: needs an Mx argument ... ENDING NOW"
  echo "  format:"
  echo "    preprocess.sh MX"
  echo "  where"
  echo "    MX        = number of grid points in x,y directions;  MX -> cell width:"
  echo "                53 -> 10mm,  105 -> 5mm, 209 -> 2.5mm, 521 -> 1mm"
  echo "  See also rungum.sh."
  exit
fi

myMx="$1"

# preprocess stage 1: create config file
echo "creating PISM-readable config override file gumparams.nc ..."
CMD="rm -f gumparams.nc"
echo $CMD
$CMD
CMD="ncgen -o gumparams.nc gumparams.cdl"
echo $CMD
$CMD

# preprocess stage 2: create bootstrap file specific to this Mx value
initfile=initlab$myMx.nc
echo "creating PISM-readable bootstrap file $initfile ..."
CMD="python buildgum.py $myMx $initfile"
echo $CMD
$CMD

