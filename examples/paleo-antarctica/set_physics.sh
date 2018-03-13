# Copyright (C) 2009-2018  PISM authors

# grids
export THIRTYKMGRID="-Mx 200 -My 200 -Lz 7000 -Lbz 2000 -Mz 41 -Mbz 16"
export TWENTYKMGRID="-Mx 300 -My 300 -Lz 7000 -Lbz 2000 -Mz 81 -Mbz 21"
export FIFTEENKMGRID="-Mx 400 -My 400 -Lz 7000 -Lbz 2000 -Mz 141 -Mbz 21"
export TWELVEKMGRID="-Mx 500 -My 500 -Lz 7000 -Lbz 2000 -Mz 141 -Mbz 31"
export TENKMGRID="-Mx 600 -My 600 -Lz 7000 -Lbz 2000 -Mz 151 -Mbz 31"
export SEVENKMGRID="-Mx 800 -My 800 -Lz 7000 -Lbz 2000 -Mz 151 -Mbz 31"
export FIVEKMGRID="-Mx 1200 -My 1200 -Lz 7000 -Lbz 2000 -Mz 201 -Mbz 51"

# skips
export SKIPTHIRTYKM=10
export SKIPTWENTYKM=10
export SKIPFIFTEENKM=10
export SKIPTWELVEKM=50
export SKIPTENKM=100
export SKIPSEVENKM=100
export SKIPFIVEKM=200


#pscale=`echo "8.2*(1.07-1.0)" | bc -l` #motivated by 7degree temperature change over 1000m height
pscale=0.574
tlapse=0.0 #already paramterized in pik_temp
export pre_atm_opts="-atmosphere pik_temp -temp_era_interim \
          -atmosphere_pik_temp_file $atmfile -surface pdd"
export fit_atm_opts="-atmosphere pik_temp,lapse_rate -temp_era_interim \
          -atmosphere_pik_temp_file $atmfile -atmosphere_lapse_rate_file $origfile \
          -temp_lapse_rate $tlapse -precip_lapse_scaling -precip_scale_factor $pscale \
          -surface pdd,forcing -force_to_thickness_file $FTTMASK -force_to_thickness_alpha 2e-4 \
          -prescribe_gl -iterative_phi $origfile \
          -tphi_inverse 500.0 -hphi_inverse 250.0 -phimax_inverse 70.0 -phimin_inverse 2.0 -phimod_inverse 2e-3"

#export atm_opts="-atmosphere pik_temp,delta_T,paleo_precip,lapse_rate -temp_era_intierim_lon \
#          -atmosphere_pik_temp_file $infile -atmosphere_delta_T_file $tforcefile \
#          -atmosphere_paleo_precip_file $tforcefile -atmosphere_lapse_rate_file $origfile \
#          -temp_lapse_rate $tlapse -precip_lapse_scaling -precip_scale_factor $pscale \
#          -surface pdd"

export atm_opts="-atmosphere pik_temp,delta_T,paleo_precip -temp_era_interim \
          -atmosphere_pik_temp_file $RESNAMEFITMOD -atmosphere_delta_T_file $tforcefile \
          -atmosphere_paleo_precip_file $tforcefile \
          -surface pdd,lapse_rate -surface_lapse_rate_file $origfile \
          -temp_lapse_rate $tlapse -smb_lapse_rate 0.0 -precip_scale_factor $pscale"


pico_opts="-gamma_T 1.0e-5 -overturning_coeff 0.8e6 -exclude_icerises \
           -number_of_basins 20 -continental_shelf_depth -2000"
export pre_ocean_opts="-ocean cavity -ocean_cavity_file $oceanfile $pico_opts"
export ocean_opts="-ocean cavity,delta_SL -ocean_cavity_file $toforcefile \
           -ocean_delta_SL_file $slforcefile $pico_opts"


export pre_calv_opts="-calving ocean_kill -ocean_kill_file $origfile"

export calv_opt="-calving eigen_calving,thickness_calving,ocean_kill \
           -ocean_kill_file $OKMASK -ocean_kill_mask \
           -eigen_calving_K 1.0e17 -thickness_calving_threshold 75.0"


export pre_bed_opts="-bed_def none -hydrology null -bed_smoother_range 5.0e3"
export bed_opts="-bed_def lc -hydrology null" # -include_ocean_load"

#already in -pik option included
subgl_opts="-subgl " #-no_subgl_basal_melt" #-tauc_slippery_grounding_lines

###### ice physics
export tillphi_param="-topg_to_phi 5.0,45.0,-300.0,700.0"
basal_opts="-yield_stress mohr_coulomb -pseudo_plastic \
            -pseudo_plastic_q 0.75 -pseudo_plastic_uthreshold 100.0" # -till_effective_fraction_overburden 0.02"

#'-pik' = '-cfbc -part_grid -kill_icebergs -subgl'
export stress_sia_opts="-pik -sia_e 2.0 -sia_flow_law gpbld"

export stress_opts="$stress_sia_opts -stress_balance ssa+sia -ssa_method fd \
             -ssa_flow_law gpbld -ssa_e 0.8 $basal_opts"

###### technical
export config_opts="-config_override ${paramfile}.nc"

export tech_opts="-options_left -verbose 2 -o_order zyx -o_size big -backup_interval 3.0"

# perhaps use these if KSP "diverged" errors occur
#export STRONGKSP="-ssafd_ksp_type gmres -ssafd_ksp_norm_type unpreconditioned -ssafd_ksp_pc_side right -ssafd_pc_type asm -ssafd_sub_pc_type lu"
#-ssafd_ksp_rtol 1e-7 -ssa_rtol 1.0e-4 -ssa_maxi 300 -ssafd_nuH_iter_failure_underrelaxation 0.8 -ssa_eps 1.0e13

