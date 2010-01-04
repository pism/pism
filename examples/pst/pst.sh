#!/bin/bash
 
# Script for Plastic till ice Stream with Thermocoupling experiment.
#
# See Bueler and Brown 2009, "Shallow shelf approximation as a
# "sliding law" in a thermomechanically coupled ice sheet model",
# J. Geophys. Res. 114, F03008, doi:10.1029/2008JF001179.
#
# Starts with a run like EISMINT II exper A (and exper I), but 
# including bedrock thermal model, inclusion of basal melt rate in mass
# continuity, and tracking of Hmelt.  Then plastic-base SSA is turned on 
# for 4 (or 3 in P2) streams in several parameter sensitivity studies.
# See src/eismint/icePSTexModel.cc for meaning of P1,P2,P3,P4, but roughly:
#   * P1 studies stream width parameter
#   * P2 studies stream orientation (relative to grid) parameter
#   * P3 studies bed slope
#   * P4 studies downstream till friction angle changes
#
# Each of P1,P2,P3,P4 is run on each of 15km,10km,7.5km,5km grids.
# The 10km grid cases are done first.  The 5km grid cases are 
# saved for last.  Each of P1,P2,P3,P4 is run on 10km horizontal
# grid with doubled vertical resolution.  Just before the 5km cases,
# the 10km grid P1 is extended to 100k years.
#
# re SPEED: On experiment P1 with a 10km grid, bueler-pogo with 8 cores 
# (two quad core Xeon processors at 2.33GHz) took about 1.5 hours for 
# 5000 model years.  If sustained and cleanly parallelizable this means
# about 400 model years per processor-hour on the 10km grid.  (We can
# determine whether it is sustained-parallel, which is likely, from timing of
# experiment P1cont below.)
#
# This script exits early so that it can be used as a standardized test of
# the ssa-as-sliding mechanism.  (There are no coupled ssa-as-sliding verification
# tests.)  Comment out the first "exit" to do full experiment.

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
        vg="-Mz 251 -Lbz 515 -Mbz 51 -z_spacing equal $3"
    elif [ "$2" = "u" ] ; then # un-equal spaced default vertical grid
        vg="-Mz 101 -Lbz 515 -Mbz 21 $3"
    elif [ "$2" = "t" ] ; then # un-equal spaced default vertical grid, *t*wice finer
        vg="-Mz 201 -Lbz 515 -Mbz 41 $3"
    else
        return 1
    fi
    
    mpst $1 "$vg"
}


#uncomment these (and "#fi" at end) and move around to bypass completed stuff:
#if [ ]; then   # always goes to "else"
#else           # put this before restart location

regridv="-regrid_vars bwat,temp,litho_temp,thk"

# P0A: run without trough on refining grid for total of 200k years:
mpst_vg $NN u "-P0A -Mx 61 -My 61 -y 1e5 -skip 20 -o P0A_100k.nc"
mpst $NN "-P0A -i P0A_100k.nc -y 50000 -skip 20 -o P0A_150k.nc"
mpst_vg $NN u "-P0A -Mx 121 -My 121 -y 40000 -ys 150000 \
 -skip 20 -regrid_from P0A_150k.nc $regridv -o P0A_190k.nc"
mpst_vg $NN u "-P0A -Mx 151 -My 151 -y 10000 -ys 190000 \
 -skip 20 -regrid_from P0A_190k.nc $regridv -o P0A.nc"


#  FIXME:  do these with -ts_file ... (and modify icePSTexModel.cc to produce)

# P1: flat with variable width but grid-aligned ice streams
mpst $NN "-P1 -i P0A.nc -ys 0 -y 5000 -skip 5 -o P1.nc"

exit # early exit so script is usable as test

# P2: flat with THREE same width and NOT grid-aligned ice streams
mpst $NN "-P2 -i P0A.nc -ys 0 -y 5000 -skip 5 -o P2.nc"

# P4: flat with variable width but grid-aligned ice streams 
#   and different down-stream till phi
mpst $NN "-P4 -i P0A.nc -ys 0 -y 5000 -skip 5 -o P4.nc"

# P0I: run with troughs on refining grid for total of 200k years:
mpst $NN "-P0I -i P0A_100k.nc -y 50000 -skip 20 -o P0I_150k.nc"
mpst_vg $NN u "-P0I -Mx 121 -My 121 -y 40000 -ys 150000 \
 -skip 20 -regrid_from P0I_150k.nc $regridv -o P0I_190k.nc"
mpst_vg $NN u "-P0I -Mx 151 -My 151 -y 10000 -ys 190000 \
 -skip 20 -regrid_from P0I_190k.nc $regridv -o P0I.nc"

# P3: troughs, with variable width but grid-aligned ice streams
mpst $NN "-P3 -i P0I.nc -ys 0 -y 5000 -skip 5 -o P3.nc"


#  grid coarsening & refinement experiments:

# COARSE: as above but on 15km grid
mpst_vg $NN u "-P1 -Mx 101 -My 101 -y 5000 -skip 5 \
 -regrid_from P0A.nc $regridv -o P1coarse.nc"

mpst_vg $NN u "-P2 -Mx 101 -My 101 -y 5000 -skip 5 \
 -regrid_from P0A.nc $regridv -o P2coarse.nc"

mpst_vg $NN u "-P3 -Mx 101 -My 101 -y 5000 -skip 5 \
 -regrid_from P0I.nc $regridv -o P3coarse.nc"  # initial state has troughs

mpst_vg $NN u "-P4 -Mx 101 -My 101 -y 5000 -skip 5 \
 -regrid_from P0A.nc $regridv -o P4coarse.nc"

# FINE: as above but on 7.5km grid
mpst_vg $NN u "-P1 -Mx 201 -My 201 -y 5000 -skip 5 \
 -regrid_from P0A.nc $regridv -o P1fine.nc"

mpst_vg $NN u "-P2 -Mx 201 -My 201 -y 5000 -skip 5 \
 -regrid_from P0A.nc $regridv -o P2fine.nc"

mpst_vg $NN u "-P3 -Mx 201 -My 201 -y 5000 -skip 5 \
 -regrid_from P0I.nc $regridv -o P3fine.nc"  # initial state has troughs

mpst_vg $NN u "-P4 -Mx 201 -My 201 -y 5000 -skip 5 \
 -regrid_from P0A.nc $regridv -o P4fine.nc"

# VERTFINE: as above but with finer vertical grid (x2 points)
mpst_vg $NN t "-P1 -Mx 151 -My 151 -y 5000 -skip 5 \
 -regrid_from P0A.nc $regridv -o P1vertfine.nc"

mpst_vg $NN t "-P2 -Mx 151 -My 151 -y 5000 -skip 5 \
 -regrid_from P0A.nc $regridv -o P2vertfine.nc"

mpst_vg $NN t "-P3 -Mx 151 -My 151 -y 5000 -skip 5 \
 -regrid_from P0I.nc $regridv -o P3vertfine.nc"  # initial state has troughs

mpst_vg $NN t "-P4 -Mx 151 -My 151 -y 5000 -skip 5 \
 -regrid_from P0A.nc $regridv -o P4vertfine.nc"


# expensive stuff follows:
#    *  100k year run of P1 on 10km grid
#    *  5km (finest) grid

#  FIXME:  do this with -save_to ...

# P1cont: as P1, but continue from saved state to 100k model years
mpst $NN "-P1 -i P1.nc -y 15000 -skip 5 -o P1_20k.nc"
mpst $NN "-P1 -i P1_20k.nc -y 20000 -skip 5 -o P1_40k.nc"
mpst $NN "-P1 -i P1_40k.nc -y 20000 -skip 5 -o P1_60k.nc"
mpst $NN "-P1 -i P1_60k.nc -y 20000 -skip 5 -o P1_80k.nc"
mpst $NN "-P1 -i P1_80k.nc -y 20000 -skip 5 -o P1_100k.nc"

# FINEST: as P1, but on 5km grid; slow!
mpst_vg $NN u "-P1 -Mx 301 -My 301 -skip 10 -y 5000 \
 -regrid_from P0A.nc $regridv -o P1finest.nc"

mpst_vg $NN u "-P2 -Mx 301 -My 301 -skip 10 -y 5000 \
 -regrid_from P0A.nc $regridv -o P2finest.nc"

mpst_vg $NN u "-P3 -Mx 301 -My 301 -skip 10 -y 5000 \
 -regrid_from P0I.nc $regridv -o P3finest.nc"  # initial state has troughs

mpst_vg $NN u "-P4 -Mx 301 -My 301 -skip 10 -y 5000 \
 -regrid_from P0A.nc $regridv -o P4finest.nc"

#fi

