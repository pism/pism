#!/bin/bash

#set -x

NN=4  # number of processors
pismdir=../../builddev
#pismdir=../../buildtill

for MM in 26 51 101 201 401;
#for MM in 11 21 41 81 161 321 641;
#for MM in 11 21 41;
do
  rm -f foo.txt
  #./runTestP.py $pismdir "mpiexec -n ${NN}" $MM &> foo.txt
  time ./runTestP.py --pism_path=$pismdir --mpiexec="mpiexec -n ${NN}" --Mx=$MM --keep &> runP$MM.txt
  echo "results for Mx=My=${MM}:"
  cat runP$MM.txt |grep -A 1 "NUMERICAL ERRORS"
done
