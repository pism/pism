#!/bin/bash
# (c) Andy Aschwanden 2018


odir=2018_05_ctrl
./psg_3d_ch.py --spatial_ts ch --restart_file 2018_05_spinup/state/storglaciaren_g40m_sb_ssa+sia_climate_init_id_CTRL_0_1000.nc --climate warming --start_year 0 --duration 100 --step 100 --o_dir $odir --exstep monthly

odir=2018_05_ch
./psg_3d_ch.py --spatial_ts ch --restart_file 2018_05_spinup/state/storglaciaren_g40m_sb_ssa+sia_climate_init_id_CTRL_0_1000.nc --climate warming --cryohydrologic_warming --start_year 0 --duration 100 --step 100 --o_dir $odir --exstep monthly
