#!/bin/bash

#set -v

MPIDO="mpiexec -n 4"
PATHTO="./"

if [ "" ]; then      # zero length string is false
#if [ "pre" ]; then  # <-- do this to create "start.n" for comparison runs
  $MPIDO ${PATHTO}pisms -Mx 121 -My 121 -Mz 101 -y 10 -o stage1.nc
  $MPIDO ${PATHTO}pismr -i stage1.nc -ye 8000 -o stage2-tillcan.nc
  cp stage2-tillcan.nc start.nc
  # assert:  dbwat == 0
  #ncap2 -o dbwat_start.nc -v -s 'dbwat=bwattillcan-bwat' start.nc
fi

# only goal of this test script is to *LOOK* at the results; verification is another matter
$MPIDO ${PATHTO}pismr -i start.nc                        -y 100 -o tillcan.nc
  # assert:  dbwat == 0
  #ncap2 -o dbwat_tillcan.nc -v -s 'dbwat=bwattillcan-bwat' tillcan.nc
$MPIDO ${PATHTO}pismr -i start.nc -diffuse_bwat          -y 100 -o diffuse.nc
  # assert:  dbwat == 0
  #ncap2 -o dbwat_diffuse.nc -v -s 'dbwat=bwattillcan-bwat' diffuse.nc
#$MPIDO ${PATHTO}pismr -i start.nc -lakes_hydrology       -y 100 -o lakes.nc
#$MPIDO ${PATHTO}pismr -i start.nc -distributed_hydrology -y 100 -o distributed.nc

#set +v

