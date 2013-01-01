#!/bin/bash

#for M in 26 51 101 201;
for M in 26 51 101 201;
do
  rm -f foo.txt
  ./test_29.py ../../buildbwp "mpiexec -n 4" $M &> foo.txt
  cat foo.txt |grep errbwat
  cat foo.txt |grep errbwp
done

