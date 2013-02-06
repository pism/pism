#!/bin/bash

NN=4  # number of processors
pismdir=../../builddev

for MM in 11 21 41 81 161 321 641;
#for MM in 11 21 41 81;
do
  rm -f foo.txt
  #./runTestP.py $pismdir "mpiexec -n ${NN}" $MM &> foo.txt
  time ./runTestP.py --pism_path=$pismdir --mpiexec="mpiexec -n ${NN}" --Mx=$MM --keep &> runP$MM.txt
  echo "results for Mx=My=${MM}:"
  cat runP$MM.txt |grep "Drift in"
done

