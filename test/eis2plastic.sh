#!/bin/bash
# script for plastic till SSA and superposition modification of
# EISMINT II experiment I (and A); note this includes bedrock thermal, unlike EIS II
# see preprint Bueler and Brown 2008, "A model of an ice sheet with one ice stream"

# re speed: with exit after completion of P5 (as below; P8, P10, and P0cont not included),
#           the whole script should take less than 500 processor-hours [???]

# re speed: on experiment P0, marmaduke.gi.alaska.edu (8 cores) took about
#           5 hours/(1000 model years) or, optimistically, about
#           ( 1 hour/(1000 model years) )/core

NN=2  # set default number of processors here

if [ $# -gt 0 ] ; then  # if user says "eis2plastic.sh 8" then NN = 8
  NN="$1"
fi

set -e  # exit on error

# function to run pisms on NN processors
mpisms()
{
    cmd="mpiexec -n $1 pisms $2"  # change if "mpirun" or "bin/pisms", etc.
    
    echo 
    echo "date = '`date`' on host '`uname -n`':"
    echo "trying '$cmd'"
    echo
    $cmd
}

# function to run pisms on NN processors and set vertical grid to standard choices
mpisms_vg()
{
    if [ "$2" = "e" ] ; then  # equal spaced default vertical grid
        vg="-Mz 251 -Mbz 51 $3"
    elif [ "$2" = "u" ] ; then # un-equal spaced default vertical grid
        vg="-Mz 101 -Mbz 41 -quadZ $3"
    else
        return 1
    fi
    
    mpisms $1 "$vg"
}


#uncomment these (and "#fi" at end) and move around to bypass completed stuff:
#if [ ]; then   # always goes to "else"
#else           # put this before restart location


#THE EXPERIMENTS:

# run without trough on coarse 25km grid for 100k years:
mpisms_vg $NN u "-eisII A -Mx 61 -My 61 -y 100000 -track_Hmelt \
 -tempskip 10 -o eis2A100k"

#continue WITHOUT trough
mpisms $NN "-eisII A -if eis2A100k.nc -y 90000 -track_Hmelt \
 -tempskip 10 -o eis2A190k"

   # refine to 12.5km grid, run 10k, save lots:
mpisms_vg $NN u "-eisII A -Mx 121 -My 121 -ys 190000 -y 10000 -track_Hmelt \
 -tempskip 10 -regrid eis2A190k.nc -regrid_vars LTBHh \
 -f3d -o eis2A_final -mato eis2A_final -matv bcYTHLCQ0345"

#continue WITH trough (starting from 100k sans trough):
mpisms $NN "-eisII I -if eis2A100k.nc -y 90000 -track_Hmelt \
 -tempskip 10 -o eis2I190k"

   # refine to 12.5km grid, run 10k, save lots:
mpisms_vg $NN u "-eisII I -Mx 121 -My 121 -ys 190000 -y 10000 -track_Hmelt \
 -tempskip 10 -regrid eis2I190k.nc -regrid_vars LTBHh \
 -f3d -o eis2I_final -mato eis2I_final -matv bcYTHLCQ0345"


#exit    # possible stopping point; uncomment to stop


# basic experiment P0 (also save 100 year and 1000 year states)
mpisms $NN "-eis2pl -if eis2I_final.nc -ys 0 -y 100 -o eis2plP0_100"

mpisms $NN "-eis2pl -if eis2plP0_100.nc -y 900 -o eis2plP0_1000"

mpisms $NN "-eis2pl -if eis2plP0_1000.nc -y 4000 -f3d \
 -o eis2plP0 -mato eis2plP0 -matv bcYTHLCQ0345"

#exit


# experiment P6 (coarser horizontal 25km grid):
mpisms_vg $NN u "-eis2pl -Mx 61 -My 61 -ys 0 -y 5000 \
 -regrid eis2I_final.nc -regrid_vars HTBL \
 -f3d -o eis2plP6 -mato eis2plP6 -matv bcYTHLCQ0345"

# experiment P7 (finer horizontal *7.5km* grid):
mpisms_vg $NN u "-eis2pl -Mx 201 -My 201 -ys 0 -y 5000 \
 -regrid eis2I_final.nc -regrid_vars HTBL \
 -f3d -o eis2plP7 -mato eis2plP7 -matv bcYTHLCQ0345"

# [finest horizontal grid run P8 put at end; slowest]

# experiment P9 (finer vertical [less than near base] 10m grid):
#MPISMS -eis2pl -Mx 121 -My 121 -Mz 501 -Mbz 101 ...
mpisms $NN "-eis2pl -Mx 121 -My 121 -Mz 201 -Mbz 81 -quadZ -ys 0 -y 5000 \
 -regrid eis2I_final.nc -regrid_vars HTBL \
 -f3d -o eis2plP9 -mato eis2plP9 -matv bcYTHLCQ0345"

# [finest vertical grid run P10 put at end; slowest]


# experiment P1 (no trough)
mpisms $NN "-eis2pl -if eis2A_final.nc -ys 0 -y 5000 \
 -no_trough -f3d -o eis2plP1 -mato eis2plP1 -matv bcYTHLCQ0345"

# experiment P2 (narrower stream):
mpisms $NN "-eis2pl -if eis2I_final.nc -ys 0 -y 5000 \
 -stream_width 50.0 -f3d -o eis2plP2 -mato eis2plP2 -matv bcYTHLCQ0345"

# experiment P3 (stronger downstream till):
mpisms $NN "-eis2pl -if eis2I_final.nc -ys 0 -y 5000 \
 -till_phi 20.0,20.0,5.0,8.0,20.0 -f3d -o eis2plP3 \
 -mato eis2plP3 -matv bcYTHLCQ0345"

# experiment P4 (lake):
mpisms $NN "-eis2pl -if eis2I_final.nc -ys 0 -y 5000 -f3d \
 -till_phi 0.0,20.0,5.0,5.0,20.0 -o eis2plP4 \
 -mato eis2plP4 -matv bcYTHLCQ0345"

# experiment P5 (fjord):
mpisms $NN "-eis2pl -if eis2I_final.nc -ys 0 -y 5000 -f3d \
 -till_phi 20.0,20.0,5.0,5.0,0.0 -o eis2plP5 \
 -mato eis2plP5 -matv bcYTHLCQ0345"


exit   # possible stopping point


# experiment P0cont (run another 95k years;  lots of runtime!):
mpisms $NN "-eis2pl -if eis2plP0.nc -y 5000 \
 -o eis2pl10k -mato eis2pl10k -matv bcYTHLCQ0345"

mpisms $NN "-eis2pl -if eis2pl10k.nc -y 10000 -o eis2pl20k"
    
mpisms $NN "-eis2pl -if eis2pl20k.nc -y 10000 -o eis2pl30k"
    
mpisms $NN "-eis2pl -if eis2pl30k.nc -y 10000 -o eis2pl40k"
    
mpisms $NN "-eis2pl -if eis2pl40k.nc -y 10000 -o eis2pl50k"
    
mpisms $NN "-eis2pl -if eis2pl50k.nc -y 10000 -o eis2pl60k"
    
mpisms $NN "-eis2pl -if eis2pl60k.nc -y 10000 -o eis2pl70k"
    
mpisms $NN "-eis2pl -if eis2pl70k.nc -y 10000 -o eis2pl80k"
    
mpisms $NN "-eis2pl -if eis2pl80k.nc -y 10000 -o eis2pl90k"

mpisms $NN "-eis2pl -if eis2pl90k.nc -y 10000 -f3d \
 -o eis2plP0cont -mato eis2plP0cont -matv bcYTHLCQ0345"


# experiment P8 (finest horizontal **5km** grid); save intermediate as it is long:
mpisms_vg $NN u "-eis2pl -Mx 301 -My 301 -ys 0 -y 1000 \
 -regrid eis2I_final.nc -regrid_vars HTBL -o eis2plP8_1k"

mpisms $NN "-eis2pl -if eis2plP8_1k.nc -y 2000 -o eis2plP8_3k"

mpisms $NN "-eis2pl -if eis2plP8_3k.nc -y 2000 -f3d -o eis2plP8 \
 -mato eis2plP8 -matv bcYTHLCQ0345"
 
 
# experiment P10 (finer vertical [less than near base] *5m* grid):
#MPISMS -eis2pl -Mx 121 -My 121 -Mz 1001 -Mbz 201 ...
mpisms $NN "-eis2pl -Mx 121 -My 121 -Mz 401 -Mbz 161 -quadZ -ys 0 -y 5000 \
 -regrid eis2I_final.nc -regrid_vars HTBL \
 -f3d -o eis2plP10 -mato eis2plP10 -matv bcYTHLCQ0345"
 
#fi

