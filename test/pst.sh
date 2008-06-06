#!/bin/bash
# 
# Script for Plastic till ice Stream with Thermocoupling experiment.   Note 
# this includes bedrock thermal and tracking of Hmelt unlike EIS II, just
# to get started.  Then plastic SSA is turned in 4 streams.
# 
# See preprint Bueler and Brown 2008, "The shallow shelf approximation as a
# `sliding law' in an ice sheet model with streaming flow".

# FIXME: need new speed estimates
# re speed: with exit after completion of P5 (as below; P8, P10, and P0cont not included),
#           the whole script should take less than 500 processor-hours [???]
# re speed: on experiment P0, marmaduke.gi.alaska.edu (8 cores) took about
#           5 hours/(1000 model years) or, optimistically, about
#           ( 1 hour/(1000 model years) )/core

NN=2  # set default number of processors here

if [ $# -gt 0 ] ; then  # if user says "pst.sh 8" then NN = 8
  NN="$1"
fi

set -e  # exit on error

# function to run "pisms -pst" on NN processors
mpst()
{
    cmd="mpiexec -n $1 pisms -pst $2"  # change if "mpirun" or "bin/pisms", etc.
    
    echo 
    echo "date = '`date`' on host '`uname -n`':"
    echo "trying '$cmd'"
    echo
    $cmd
}

# function to run pisms on NN processors and set vertical grid to standard choices
mpst_vg()
{
    if [ "$2" = "e" ] ; then  # equal spaced default vertical grid
        vg="-Mz 251 -Mbz 51 $3"
    elif [ "$2" = "u" ] ; then # un-equal spaced default vertical grid
        vg="-Mz 101 -Mbz 41 -quadZ $3"
    else
        return 1
    fi
    
    mpst $1 "$vg"
}


#uncomment these (and "#fi" at end) and move around to bypass completed stuff:
#if [ ]; then   # always goes to "else"
#else           # put this before restart location


# P0A: run without trough on refining grid for total of 200k years:
mpst_vg $NN u "-P0A -Mx 61 -My 61 -y 100000 -tempskip 10 -o P0A_100k"

mpst $NN "-P0A -if P0A_100k.nc -y 50000 -tempskip 10 -o P0A_150k"

mpst_vg $NN u "-P0A -Mx 121 -My 121 -y 40000 \
 -tempskip 10 -regrid P0A_150k.nc -regrid_vars LTBHhe -o P0A_190k"

mpst_vg $NN u "-P0A -Mx 151 -My 151 -y 10000 \
 -tempskip 10 -regrid P0A_190k.nc -regrid_vars LTBHhe -f3d -o P0A"


# P1 (also save 100 year and 1000 year states): flat with variable
#   width but grid-aligned ice streams, on 10km grid
mpst $NN "-P1 -if P0A.nc -ys 0 -y 100 -o P1_100"

mpst $NN "-P1 -if P1_100.nc -y 900 -o P1_1000"

mpst $NN "-P1 -if P1_1000.nc -y 4000 -f3d -o P1"

# P2: flat with THREE same width and NOT grid-aligned ice streams
mpst $NN "-P2 -if P0A.nc -ys 0 -y 5000 -o P2"

# P4: flat with variable width but grid-aligned ice streams 
#   and different down-stream till phi
mpst $NN "-P4 -if P0A.nc -ys 0 -y 5000 -o P4"


# P0I: run with trough on refining grid for total of 200k years:
mpst_vg $NN u "-P0I -Mx 61 -My 61 -y 100000 -tempskip 10 -o P0I_100k"

mpst $NN "-P0I -if P0I_100k.nc -y 50000 -tempskip 10 -o P0I_150k"

mpst_vg $NN u "-P0I -Mx 121 -My 121 -y 40000 \
 -tempskip 10 -regrid P0I_150k.nc -regrid_vars LTBHhe -o P0I_190k"

mpst_vg $NN u "-P0I -Mx 151 -My 151 -y 10000 \
 -tempskip 10 -regrid P0I_190k.nc -regrid_vars LTBHhe -f3d -o P0I"


# P3: troughs, with variable width but grid-aligned ice streams
mpst $NN "-P3 -if P0I.nc -ys 0 -y 5000 -o P3"


exit


# P1coarse: as P1, but on 15km grid
mpst_vg $NN u "-P1 -Mx 101 -My 101 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P1coarse"


# P1fine: as P1, but on 7.5km grid
mpst_vg $NN u "-P1 -Mx 201 -My 201 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P1fine"


# P1vertfine: as P1, but with finer vertical grid (x2 points)
mpst $NN u "-P1 -Mx 151 -My 151 -Mz 201 -Mbz 81 -quadZ -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P1vertfine"

 
# P1finest: as P1, but on 5km grid; slow so save intermediate
mpst_vg $NN u "-P1 -Mx 301 -My 301 -y 1000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P1finest_1k"

mpst $NN "-P1 -if P1finest_1k -y 2000 -o P1finest_3k"

mpst $NN "-P1 -if P1finest_3k -y 2000 -o P1finest"


# P1cont: as P1, but continue from saved state to 100k model years
mpst $NN "-P1 -if P1.nc -y 15000 -o P1_20k"

mpst $NN "-P1 -if P1_20k.nc -y 20000 -o P1_40k"

mpst $NN "-P1 -if P1_40k.nc -y 20000 -o P1_60k"

mpst $NN "-P1 -if P1_60k.nc -y 20000 -o P1_80k"

mpst $NN "-P1 -if P1_80k.nc -y 20000 -o P1_100k"



