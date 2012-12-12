#!/bin/bash
# Copyright (C) 2011 Torsten Albrecht (torsten.albrecht@pik-potsdam.de)
# regression/verification test of SSA ice shelf solution against Van-der-Veen analytical solution 

count=0
boxes=("11" "25" "51" "101" "251" "501") #cells

for res in 50 20 10 5 2 1  #km
do

../../../bin/pismr -verbose 3 -boot_file flowline_setup_${res}km.nc -Mx ${boxes[$count]} -My 3 -Mz 51 -Lz 1000 -cold -ssa_flow_law isothermal_glen -no_energy -ssa_sliding -ssa_dirichlet_bc -config_override flowline_config.nc -ssa_method fd -cfbc -part_grid -part_redist -ksp_rtol 1e-7 -ys 0 -ye 3e3 -options_left -no_sia -extra_file flowline_extra_${res}km.nc -extra_times 0:25:3e3 -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,usurf,velbase,velsurf,velbar,bcflag -ts_file flowline_time_${res}km.nc -ts_times 0:25:3e3 -o flowline_result_${res}km.nc -o_order zyx -o_size big -thickness_calving -calving_at_thickness 250 > flowline_out_${res}km.out

#-pik: -cfbc -part_grid -part_redist -kill_icebergs

let count+=1

done

