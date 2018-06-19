#!/bin/bash
# (c) Andy Aschwanden 2018

myyear=$(date +'%Y')
mymonth=$(date +'%m')

if [ -z "$1" ]; then
    N=4
else
    N=$1
fi

# A simulation with 1K warming over a hundred years, control simulation
odir=${myyear}_${mymonth}_ctrl
./psg_3d_ch.py  -n $N --spatial_ts default --restart_file 2018_05_spinup/state/storglaciaren_g40m_sb_ssa+sia_climate_init_id_CTRL_0_1000.nc --climate warming --start_year 0 --duration 100 --step 100 --o_dir $odir --exstep monthly

# A simulation with 1K warming over a hundred years, with cryo-hydrologic warming
odir=${myyear}_${mymonth}_ch
./psg_3d_ch.py  -n $N --spatial_ts ch --restart_file 2018_05_spinup/state/storglaciaren_g40m_sb_ssa+sia_climate_init_id_CTRL_0_1000.nc --climate warming --cryohydrologic_warming --start_year 0 --duration 100 --step 100 --o_dir $odir --exstep monthly
