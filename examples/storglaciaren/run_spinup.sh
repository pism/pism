#!/bin/bash
# (c) Andy Aschwanden 2018


odir=2018_05_spinup
./preprocess

./psg_3d_ch.py -b --stress_balance sia --start_year 0 --duration 250 --step 250 --o_dir $odir --exstep 10

./psg_3d_ch.py --restart_file $odir/state/storglaciaren_g40m_sb_ssa+sia_climate_init_id_CTRL_0_250.nc --start_year 0 --duration 1000 --step 1000 --o_dir $odir --exstep 25

./psg_3d_ch.py --restart_file $odir/state/storglaciaren_g40m_sb_ssa+sia_climate_init_id_CTRL_0_1000.nc --climate elev+ftt --start_year 0 --duration 1000 --step 1000 --o_dir $odir --exstep 25
