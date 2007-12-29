#!/bin/bash
# script for plastic till SSA and superposition modification of
# EISMINT II experiment I (and A); note this includes bedrock thermal, unlike EIS II
# see preprint Bueler and Brown 2007, "A model of an ice sheet with one ice stream"

# re speed: on experiment P0, marmaduke.gi.alaska.edu (8 cores) took about
# 5 hours/(1000 model years) or, optimistically, about
# ( 1 hour/(1000 model years) )/core

# this alias makes it easier to move to other machines w. "mpirun" or unset path
shopt -s expand_aliases
alias MPISMS='mpiexec -n 2 pisms '

#uncomment these (and "#fi" at end) and move around to bypass completed stuff:
#if [ ]; then   # always goes to "else"
#else           # put this before restart location

  
# run without trough on coarse 25km grid for 100k years:
MPISMS -eisII A -Mx 61 -My 61 -Mz 251 -Mbz 51 -y 100000 -track_Hmelt \
   -o eis2A100k

#continue WITHOUT trough
MPISMS -eisII A -if eis2A100k.nc -y 90000 -track_Hmelt -o eis2A190k

   # refine to 12.5km grid, run 10k, save lots:
MPISMS -eisII A -Mx 121 -My 121 -Mz 251 -Mbz 51 \
   -regrid eis2A190k.nc -regrid_vars LTBHh \
   -ys 190000 -y 10000 -track_Hmelt -f3d -o eis2A_final \
   -mato eis2A_final -matv bcYTHLCQ0345

#continue WITH trough (starting from 100k sans trough):
MPISMS -eisII I -if eis2A100k.nc -y 90000 -track_Hmelt -o eis2I190k

   # refine to 12.5km grid, run 10k, save lots:
MPISMS -eisII I -Mx 121 -My 121 -Mz 251 -Mbz 51 \
   -regrid eis2I190k.nc -regrid_vars LTBHh \
   -ys 190000 -y 10000 -track_Hmelt -f3d -o eis2I_final \
   -mato eis2I_final -matv bcYTHLCQ0345


#exit    # possible stopping point; uncomment to stop


# basic experiment P0 (also save 100 year and 1000 year states)
MPISMS -eis2pl -if eis2I_final.nc -ys 0 -y 100 -o eis2plP0_100

MPISMS -eis2pl -if eis2plP0_100.nc -y 900 -o eis2plP0_1000

MPISMS -eis2pl -if eis2plP0_1000.nc -y 4000 -f3d \
    -o eis2plP0 -mato eis2plP0 -matv bcYTHLCQ0345


#exit   # possible stopping point


# experiment P6 (coarser horizontal 25km grid):
MPISMS -eis2pl -Mx 61 -My 61 -Mz 251 -Mbz 51 -ys 0 -y 5000 -f3d \
    -regrid eis2I_final.nc -regrid_vars HTBL -o eis2plP6 \
    -mato eis2plP6 -matv bcYTHLCQ0345

# experiment P7 (finer horizontal *7.5km* grid):
MPISMS -eis2pl -Mx 201 -My 201 -Mz 251 -Mbz 51 -ys 0 -y 5000 -f3d \
    -regrid eis2I_final.nc -regrid_vars HTBL -o eis2plP7 \
    -mato eis2plP7 -matv bcYTHLCQ0345

# experiment P8 (finer horizontal **5km** grid); save intermediate as it is long:
MPISMS -eis2pl -Mx 301 -My 301 -Mz 251 -Mbz 51 -ys 0 -y 1000 \
    -regrid eis2I_final.nc -regrid_vars HTBL -o eis2plP8_1k

MPISMS -eis2pl -if eis2plP8_1k.nc -y 2000 -o eis2plP8_3k

MPISMS -eis2pl -if eis2plP8_3k.nc -y 2000 -f3d -o eis2plP8 \
    -mato eis2plP8 -matv bcYTHLCQ0345

# experiment P9 (finer vertical 10m grid):
MPISMS -eis2pl -Mx 121 -My 121 -Mz 501 -Mbz 101 -ys 0 -y 5000 -f3d \
    -regrid eis2I_final.nc -regrid_vars HTBL -o eis2plP9 \
    -mato eis2plP9 -matv bcYTHLCQ0345

# experiment P10 (finer vertical *5m* grid):
MPISMS -eis2pl -Mx 121 -My 121 -Mz 1001 -Mbz 201 -ys 0 -y 5000 -f3d \
    -regrid eis2I_final.nc -regrid_vars HTBL -o eis2plP10 \
    -mato eis2plP10 -matv bcYTHLCQ0345


#exit    # possible stopping point


# experiment P1 (no trough)
MPISMS -eis2pl -if eis2A_final.nc -ys 0 -y 5000 \
    -no_trough -f3d -o eis2plP1 -mato eis2plP1 -matv bcYTHLCQ0345

# experiment P2 (narrower stream):
MPISMS -eis2pl -if eis2I_final.nc -ys 0 -y 5000 \
    -stream_width 50.0 -f3d -o eis2plP2 -mato eis2plP2 -matv bcYTHLCQ0345

# experiment P3 (stronger downstream till):
MPISMS -eis2pl -if eis2I_final.nc -ys 0 -y 5000 \
    -till_phi 20.0,20.0,5.0,8.0,20.0 -f3d -o eis2plP3 \
    -mato eis2plP3 -matv bcYTHLCQ0345

# experiment P4 (lake):
MPISMS -eis2pl -if eis2I_final.nc -ys 0 -y 5000 -f3d \
    -till_phi 0.0,20.0,5.0,5.0,20.0 -o eis2plP4 \
    -mato eis2plP4 -matv bcYTHLCQ0345

# experiment P5 (fjord):
MPISMS -eis2pl -if eis2I_final.nc -ys 0 -y 5000 -f3d \
    -till_phi 20.0,20.0,5.0,5.0,0.0 -o eis2plP5 \
    -mato eis2plP5 -matv bcYTHLCQ0345


#exit   # possible stopping point


# experiment P0cont (run another 95k years;  lots of runtime!):
MPISMS -eis2pl -if eis2plP0.nc -y 5000 \
    -o eis2pl10k -mato eis2pl10k -matv bcYTHLCQ0345

MPISMS -eis2pl -if eis2pl10k.nc -y 10000 -o eis2pl20k
    
MPISMS -eis2pl -if eis2pl20k.nc -y 10000 -o eis2pl30k
    
MPISMS -eis2pl -if eis2pl30k.nc -y 10000 -o eis2pl40k
    
MPISMS -eis2pl -if eis2pl40k.nc -y 10000 -o eis2pl50k
    
MPISMS -eis2pl -if eis2pl50k.nc -y 10000 -o eis2pl60k
    
MPISMS -eis2pl -if eis2pl60k.nc -y 10000 -o eis2pl70k
    
MPISMS -eis2pl -if eis2pl70k.nc -y 10000 -o eis2pl80k
    
MPISMS -eis2pl -if eis2pl80k.nc -y 10000 -o eis2pl90k

MPISMS -eis2pl -if eis2pl90k.nc -y 10000 -f3d \
    -o eis2plP0cont -mato eis2plP0cont -matv bcYTHLCQ0345

#fi

