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


#THE EXPERIMENTS:

# P0A: run without trough on refining grid for total of 200k years:
mpst_vg $NN u "-P0A -Mx 61 -My 61 -y 100000 -track_Hmelt \
 -tempskip 10 -o P0A_100k"

mpst $NN "-P0A -if P0A_100k.nc -y 50000 -track_Hmelt \
 -tempskip 10 -o P0A_150k"

mpst_vg $NN u "-P0A -Mx 121 -My 121 -y 40000 -track_Hmelt \
 -tempskip 10 -regrid P0A_150k.nc -regrid_vars LTBHhe -o P0A_190k"

mpst_vg $NN u "-P0A -Mx 151 -My 151 -y 10000 -track_Hmelt \
 -tempskip 10 -regrid P0A_190k.nc -regrid_vars LTBHhe -f3d -o P0A"


# P0I: run with trough on refining grid for total of 200k years:
mpst_vg $NN u "-P0I -Mx 61 -My 61 -y 100000 -track_Hmelt \
 -tempskip 10 -o P0I_100k"

mpst $NN "-P0I -if P0I_100k.nc -y 50000 -track_Hmelt \
 -tempskip 10 -o P0I_150k"

mpst_vg $NN u "-P0I -Mx 121 -My 121 -y 40000 -track_Hmelt \
 -tempskip 10 -regrid P0I_150k.nc -regrid_vars LTBHhe -o P0I_190k"

mpst_vg $NN u "-P0I -Mx 151 -My 151 -y 10000 -track_Hmelt \
 -tempskip 10 -regrid P0I_190k.nc -regrid_vars LTBHhe -f3d -o P0I"


#exit    # possible stopping point; uncomment to stop


# P1 (also save 100 year and 1000 year states): flat with variable
#   width but grid-aligned ice streams
mpst $NN "-P1 -if P0A.nc -ys 0 -y 100 -o P1_100"

mpst $NN "-P1 -if P1_100.nc -y 900 -o P1_1000"

mpst $NN "-P1 -if P1_1000.nc -y 4000 -f3d -o P1"

# P2: flat with THREE same width and NOT grid-aligned ice streams
mpst $NN "-P2 -if P0A.nc -ys 0 -y 5000 -o P2"

# P3: troughs, with variable width but grid-aligned ice streams
mpst $NN "-P3 -if P0I.nc -ys 0 -y 5000 -o P3"

# P4: flat with variable width but grid-aligned ice streams 
#   and different down-stream till phi
mpst $NN "-P4 -if P0A.nc -ys 0 -y 5000 -o P4"

exit


## experiment P6 (coarser horizontal 25km grid):
mpisms_vg $NN u "-eis2pl -Mx 61 -My 61 -ys 0 -y 5000 \
 -regrid eis2I_final.nc -regrid_vars LTBHh -f3d -o eis2plP6"

# experiment P7 (finer horizontal *7.5km* grid):
mpisms_vg $NN u "-eis2pl -Mx 201 -My 201 -ys 0 -y 5000 \
 -regrid eis2I_final.nc -regrid_vars LTBHh -f3d -o eis2plP7"

# [finest horizontal grid (5km) run P8 put at end; slowest]

# experiment P9 (finer vertical [less than near base] 10m grid):
mpisms $NN "-eis2pl -Mx 121 -My 121 -Mz 201 -Mbz 81 -quadZ -ys 0 -y 5000 \
 -regrid eis2I_final.nc -regrid_vars LTBHh -f3d -o eis2plP9"

# [finest vertical grid run P10 put at end; slowest]



# experiment P1 (no trough)
mpisms $NN "-eis2pl -if eis2A_final.nc -ys 0 -y 5000 \
 -no_trough -f3d -o eis2plP1"

# experiment P2 (narrower stream):
mpisms $NN "-eis2pl -if eis2I_final.nc -ys 0 -y 5000 \
 -stream_width 50.0 -f3d -o eis2plP2"

# experiment P3 (stronger downstream till):
mpisms $NN "-eis2pl -if eis2I_final.nc -ys 0 -y 5000 \
 -till_phi 20.0,20.0,5.0,8.0,8.0 -f3d -o eis2plP3"

# experiment P4 (lake):
mpisms $NN "-eis2pl -if eis2I_final.nc -ys 0 -y 5000 -f3d \
 -till_phi 0.0,20.0,5.0,5.0,5.0 -o eis2plP4"

# experiment P5 (fjord):
mpisms $NN "-eis2pl -if eis2I_final.nc -ys 0 -y 5000 -f3d \
 -till_phi 20.0,20.0,5.0,5.0,0.0 -o eis2plP5"


exit   # possible stopping point


# experiment P0cont (run another 95k years;  lots of runtime!):
mpisms $NN "-eis2pl -if eis2plP0.nc -y 5000"

mpisms $NN "-eis2pl -if eis2pl10k.nc -y 10000 -o eis2pl20k"
    
mpisms $NN "-eis2pl -if eis2pl20k.nc -y 10000 -o eis2pl30k"
    
mpisms $NN "-eis2pl -if eis2pl30k.nc -y 10000 -o eis2pl40k"
    
mpisms $NN "-eis2pl -if eis2pl40k.nc -y 10000 -o eis2pl50k"
    
mpisms $NN "-eis2pl -if eis2pl50k.nc -y 10000 -o eis2pl60k"
    
mpisms $NN "-eis2pl -if eis2pl60k.nc -y 10000 -o eis2pl70k"
    
mpisms $NN "-eis2pl -if eis2pl70k.nc -y 10000 -o eis2pl80k"
    
mpisms $NN "-eis2pl -if eis2pl80k.nc -y 10000 -o eis2pl90k"

mpisms $NN "-eis2pl -if eis2pl90k.nc -y 10000 -f3d -o eis2plP0cont"


# experiment P8 (finest horizontal **5km** grid); save intermediate as it is long:
mpisms_vg $NN u "-eis2pl -Mx 301 -My 301 -ys 0 -y 1000 \
 -regrid eis2I_final.nc -regrid_vars LTBHh -o eis2plP8_1k"

mpisms $NN "-eis2pl -if eis2plP8_1k.nc -y 2000 -o eis2plP8_3k"

mpisms $NN "-eis2pl -if eis2plP8_3k.nc -y 2000 -f3d -o eis2plP8"
 
 
# experiment P10 (finest vertical [less than near base] *5m* grid):
#MPISMS -eis2pl -Mx 121 -My 121 -Mz 1001 -Mbz 201 ...
mpisms $NN "-eis2pl -Mx 121 -My 121 -Mz 401 -Mbz 161 -quadZ -ys 0 -y 5000 \
 -regrid eis2I_final.nc -regrid_vars LTBHh -f3d -o eis2plP10"
 
#fi

