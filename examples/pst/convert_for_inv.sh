#!/bin/bash

# this script is related to the development of the -surf_vel_to_phi option
#   and the associated inverse model

# generate INFILE below by running pst.sh to get P0A.nc, and then:
#   $  pisms -pst -P1 -if P0A.nc -y 1 -pseudo_plastic_q 0.25 -o P1earlypseudo.nc # writes a tillphi
#   $  pismr -ssa -super -plastic -pseudo_plastic_q 0.25 -eta -if P1earlypseudo.nc -y 199 -f3d -o basestate.nc
# then do:
#   $  ./convert_for_inv.sh
# and try inverse method versus no inverse:
#   $  pismr -ssa -super -plastic -pseudo_plastic_q 0.25 -eta -if inv_me.nc -y 100 \
#   >   -surf_vel_to_phi inv_me.nc -inv_reg_eps 0.0 -inv_write_fields ifields_noreg.nc -o result_noreg.nc
#   $  pismr -ssa -super -plastic -pseudo_plastic_q 0.25 -eta -if inv_me.nc -y 100 \
#   >   -o result_noinv.nc

# also see comments in src/base/iMinverse.cc

INFILE=basestate.nc
OUTFILE=inv_me.nc

cp $INFILE $OUTFILE 

# make points with no ice invalid by making velocities huge
#   (see comments for IceModel::invertSurfaceVelocities())
# ncap2 mis-handles attributes?
ncap2 -O -s 'obsmask=(thk>1200.0)' $OUTFILE $OUTFILE
ncap2 -O -s 'uvelsurf=obsmask*uvelsurf+!obsmask*(2.0*31556926.0)' $OUTFILE $OUTFILE
ncap2 -O -s 'vvelsurf=obsmask*vvelsurf+!obsmask*(2.0*31556926.0)' $OUTFILE $OUTFILE

# convert to m s-1; ncap2 fails?
ncap -O -s 'uvelsurf=uvelsurf/31556926.0' $OUTFILE $OUTFILE
ncatted -O -a units,uvelsurf,c,c,'m s-1' $OUTFILE
ncap -O -s 'vvelsurf=vvelsurf/31556926.0' $OUTFILE $OUTFILE
ncatted -O -a units,vvelsurf,c,c,'m s-1' $OUTFILE

ncatted -O -a _FillValue,uvelsurf,c,d,2.0 $OUTFILE # 2.0 m/s is huge and out of valid range
ncatted -O -a valid_min,uvelsurf,c,d,-0.1 $OUTFILE
ncatted -O -a valid_max,uvelsurf,c,d,0.1  $OUTFILE
ncatted -O -a _FillValue,vvelsurf,c,d,2.0 $OUTFILE
ncatted -O -a valid_min,vvelsurf,c,d,-0.1 $OUTFILE
ncatted -O -a valid_max,vvelsurf,c,d,0.1  $OUTFILE

ncap -O -s 'obsmask=1.0*obsmask' $OUTFILE $OUTFILE # clears hosed attributes

