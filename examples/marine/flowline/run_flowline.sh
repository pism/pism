#!/bin/bash
# Copyright (C) 2011, 2013, 2014, 2015, 2016, 2021 Torsten Albrecht (torsten.albrecht@pik-potsdam.de)
# regression/verification test of SSA ice shelf solution against Van-der-Veen analytical solution

count=0
boxes=("11" "25" "51" "101" "251" "501") #cells

for res in 50 20 10 5 2 1  #km
do
    pismr -verbose 2 -i flowline_setup_${res}km.nc -bootstrap -Mx ${boxes[$count]} -My 3 -Mz 5 -Lz 1000 -energy cold -ssa_flow_law isothermal_glen -energy none -stress_balance ssa -ssa_dirichlet_bc -config_override flowline_config.nc -ssa_method fd -cfbc -part_grid -ssafd_ksp_rtol 1e-7 -ys 0 -ye 3e3 -options_left -extra_file flowline_extra_${res}km.nc -extra_times 0:25:3e3 -extra_vars thk,topg,velbar_mag,flux_mag,mask,dHdt,usurf,hardav,usurf,velbase,velsurf,velbar,vel_bc_mask -ts_file flowline_time_${res}km.nc -ts_times 0:25:3e3 -o flowline_result_${res}km.nc -o_order zyx -o_size big -calving thickness_calving -thickness_calving_threshold 250 | tee flowline_out_${res}km.out

    let count+=1
done
