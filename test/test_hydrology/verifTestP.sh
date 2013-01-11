#!/bin/bash

NN=4  # number of processors
pismdir=../../buildbwp

#for MM in 11 21 41 81 161 321;
for MM in 11 21 41 81;
do
  rm -f foo.txt
  #./runTestP.py $pismdir "mpiexec -n ${NN}" $MM &> foo.txt
  ./runTestP.py --pism_path=$pismdir --mpiexec="mpiexec -n ${NN}" --Mx=$MM &> foo.txt
  echo "results for Mx=My=${MM}:"
  cat foo.txt |grep "Drift in"
done

