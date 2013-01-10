#!/bin/bash

NN=4  # number of processors

for MM in 11 21 41 81 161 321;
do
  rm -f foo.txt
  ./test_29.py ../../buildbwp "mpiexec -n ${NN}" $MM &> foo.txt
  echo "results for Mx=My=${MM}:"
  cat foo.txt |grep errbw
done

