#!/bin/bash
 
# Script for Plastic till ice Stream with Thermocoupling experiment.  See
# preprint Bueler and Brown 2008, "The shallow shelf approximation as a
# `sliding law' in an ice sheet model with streaming flow".
#
# Starts with a run like EISMINT II exper A (and exper I), but 
# this includes bedrock thermal model, inclusion of basal melt rate in mass
# continuity, and tracking of Hmelt.  Then plastic SSA is turned on 
# for 4 (or 3 in P2) streams in several parameter sensitivity studies.
# See src/eismint/icePSTexModel.cc for meaning of P1,P2,P3,P4, but roughly:
#   * P1 studies stream width parameter
#   * P2 studies stream orientation (relative to grid) parameter
#   * P3 studies bed slope
#   * P4 studies downstream till friction angle changes
#
# Each of P1,P2,P3,P4 is run on each of 15km,10km,7.5km,5km grids.
# The 10km grid cases are done first.  The 5km grid cases are saved for last.
# Each of P1,P2,P3,P4 is run on 10km horizontal grid with doubled vertical resolution.
# Just before the 5km cases, the 10km grid P1 is extended to 100k years.
#
# re SPEED: On experiment P1 with a 10km grid, bueler-pogo with 8 cores 
# (two quad core Xeon processors at 2.33GHz) took about 1.5 hours for 
# 5000 model years.  If sustained and cleanly parallelizable this means
# about 400 model years per processor-hour on the 10km grid.  (We can
# determine whether it is sustained-parallel, which is likely, from timing of
# experiment P1cont below.)

NN=2  # set default number of processors here

if [ $# -gt 0 ] ; then  # if user says "pst.sh 8" then NN = 8
  NN="$1"
fi

SHOWONLY=0
if [ $# -gt 1 ] ; then  # if user says "pst.sh 8 D" then NN = 8 and only shows; no run 
  SHOWONLY=1
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
    if [ $SHOWONLY = 0 ] ; then
      $cmd
    fi
}

# function to run pisms on NN processors and set vertical grid to standard choices
mpst_vg()
{
    if [ "$2" = "e" ] ; then  # equal spaced default vertical grid
        vg="-Mz 251 -Mbz 51 $3"
    elif [ "$2" = "u" ] ; then # un-equal spaced default vertical grid
        vg="-Mz 101 -Mbz 41 -quadZ $3"
    elif [ "$2" = "t" ] ; then # un-equal spaced default vertical grid, *t*wice finer
        vg="-Mz 201 -Mbz 81 -quadZ $3"
    else
        return 1
    fi
    
    mpst $1 "$vg"
}


#uncomment these (and "#fi" at end) and move around to bypass completed stuff:
#if [ ]; then   # always goes to "else"
#else           # put this before restart location


# P0A: run without trough on refining grid for total of 200k years:
mpst_vg $NN u "-P0A -Mx 61 -My 61 -y 1e5 -skip 20 -o P0A_100k"
mpst $NN "-P0A -if P0A_100k.nc -y 50000 -skip 20 -o P0A_150k"
mpst_vg $NN u "-P0A -Mx 121 -My 121 -y 40000 \
 -skip 20 -regrid P0A_150k.nc -regrid_vars LTBHhe -o P0A_190k"
mpst_vg $NN u "-P0A -Mx 151 -My 151 -y 10000 \
 -skip 20 -regrid P0A_190k.nc -regrid_vars LTBHhe -o P0A"


# P1: flat with variable width but grid-aligned ice streams
mpst $NN "-P1 -if P0A.nc -ys 0 -y 5000 -o P1"

# P2: flat with THREE same width and NOT grid-aligned ice streams
mpst $NN "-P2 -if P0A.nc -ys 0 -y 5000 -o P2"

# P4: flat with variable width but grid-aligned ice streams 
#   and different down-stream till phi
mpst $NN "-P4 -if P0A.nc -ys 0 -y 5000 -o P4"


# P0I: run with troughs on refining grid for total of 200k years:
mpst $NN "-P0I -if P0A_100k.nc -y 50000 -skip 20 -o P0I_150k"
mpst_vg $NN u "-P0I -Mx 121 -My 121 -y 40000 \
 -skip 20 -regrid P0I_150k.nc -regrid_vars LTBHhe -o P0I_190k"
mpst_vg $NN u "-P0I -Mx 151 -My 151 -y 10000 \
 -skip 20 -regrid P0I_190k.nc -regrid_vars LTBHhe -o P0I"


# P3: troughs, with variable width but grid-aligned ice streams
mpst $NN "-P3 -if P0I.nc -ys 0 -y 5000 -o P3"


# possible stopping point before grid coarsening & refinement experiments
#exit  


# COARSE: as above but on 15km grid
mpst_vg $NN u "-P1 -Mx 101 -My 101 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P1coarse"

mpst_vg $NN u "-P2 -Mx 101 -My 101 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P2coarse"

mpst_vg $NN u "-P3 -Mx 101 -My 101 -y 5000 \
 -regrid P0I.nc -regrid_vars LTBHhe -o P3coarse"  # initial state has troughs

mpst_vg $NN u "-P4 -Mx 101 -My 101 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P4coarse"

#exit

# FINE: as above but on 7.5km grid
mpst_vg $NN u "-P1 -Mx 201 -My 201 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P1fine"

mpst_vg $NN u "-P2 -Mx 201 -My 201 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P2fine"

mpst_vg $NN u "-P3 -Mx 201 -My 201 -y 5000 \
 -regrid P0I.nc -regrid_vars LTBHhe -o P3fine"

mpst_vg $NN u "-P4 -Mx 201 -My 201 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P4fine"


# VERTFINE: as above but with finer vertical grid (x2 points)
mpst_vg $NN t "-P1 -Mx 151 -My 151 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P1vertfine"

mpst_vg $NN t "-P2 -Mx 151 -My 151 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P2vertfine"

mpst_vg $NN t "-P3 -Mx 151 -My 151 -y 5000 \
 -regrid P0I.nc -regrid_vars LTBHhe -o P3vertfine"

mpst_vg $NN t "-P4 -Mx 151 -My 151 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P4vertfine"


# possible stopping point before P1 long run
exit

# P1cont: as P1, but continue from saved state to 100k model years
mpst $NN "-P1 -if P1.nc -y 15000 -o P1_20k"
mpst $NN "-P1 -if P1_20k.nc -y 20000 -o P1_40k"
mpst $NN "-P1 -if P1_40k.nc -y 20000 -o P1_60k"
mpst $NN "-P1 -if P1_60k.nc -y 20000 -o P1_80k"
mpst $NN "-P1 -if P1_80k.nc -y 20000 -o P1_100k"


# possible stopping point before finest grid
exit

# FINEST: as P1, but on 5km grid; slow!
mpst_vg $NN u "-P1 -Mx 301 -My 301 -tempskip 10 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P1finest"

mpst_vg $NN u "-P2 -Mx 301 -My 301 -tempskip 10 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P2finest"

mpst_vg $NN u "-P3 -Mx 301 -My 301 -tempskip 10 -y 5000 \
 -regrid P0I.nc -regrid_vars LTBHhe -o P3finest"

mpst_vg $NN u "-P4 -Mx 301 -My 301 -tempskip 10 -y 5000 \
 -regrid P0A.nc -regrid_vars LTBHhe -o P4finest"


#fi

