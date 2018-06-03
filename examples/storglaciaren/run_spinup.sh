#!/bin/bash
# (c) Andy Aschwanden 2018

myyear=$(date +'%Y')
mymonth=$(date +'%m')

if [ -z "$1" ]; then
    N=4
else
    N=$1
fi

./preprocess.sh

odir=${myyear}_${mymonth}_spinup

./psg_3d_ch.py -n $N -b --stress_balance sia --start_year 0 --duration 250 --step 250 --o_dir $odir --exstep 10

./psg_3d_ch.py -n $N --restart_file $odir/state/storglaciaren_g40m_sb_ssa+sia_climate_init_id_CTRL_0_250.nc --start_year 0 --duration 1000 --step 1000 --o_dir $odir --exstep 25

./psg_3d_ch.py -n $N --restart_file $odir/state/storglaciaren_g40m_sb_ssa+sia_climate_init_id_CTRL_0_1000.nc --climate elev+ftt --start_year 0 --duration 1000 --step 1000 --o_dir $odir --exstep 25
