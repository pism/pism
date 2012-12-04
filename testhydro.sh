#!/bin/bash

#set -v

MPIDO="mpiexec -n 4"
PATHTO="./"

$MPIDO ${PATHTO}pisms -o stage1.nc -y 10
$MPIDO ${PATHTO}pismr -i stage1.nc -ye 6000 -o stage2-tillcan.nc
$MPIDO ${PATHTO}pismr -i stage1.nc -diffuse_bwat -ye 6000 -o stage2-diffusebwat.nc

#set +v

