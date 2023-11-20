#!/bin/bash

# Copyright (C) 2013  Ed Bueler

SCRIPTNAME="#(topgreplace.sh)"

NAME=$1
DO=$PISM_DO

# replace topg with topgsmooth (does nothing if $NAME is absent)
echo
echo "$SCRIPTNAME  replacing topg in $NAME with topgsmooth"
cmd="ncatted -O -a standard_name,topg,d,, $NAME"
$DO $cmd
cmd="ncrename -O -v topg,topgoriginal $NAME"
$DO $cmd
cmd="ncrename -O -v topgsmooth,topg $NAME"
$DO $cmd
cmd="ncatted -O -a standard_name,topg,o,c,bedrock_altitude $NAME"
$DO $cmd
cmd="ncatted -O -a long_name,topg,d,, $NAME"
$DO $cmd
cmd="ncatted -O -a pism_intent,topg,d,, $NAME"
$DO $cmd
echo
