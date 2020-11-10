# Copyright (C) 2009-2018  PISM authors

# grids
export THIRTYKMGRID="-Mx 200 -My 200 -Lz 7000 -Lbz 2000 -Mz 61 -Mbz 16"

# skips
export SKIPTHIRTYKM=10


#MODEL PHYSICS OPTIONS
pscale=0.574
tlapse=8.0

# ATMOSPHERE
export pre_atm_opts="-atmosphere given,elevation_change -atmosphere_given_period 1 \
          -atmosphere_given_file $atmfile -atmosphere_elevation_change_file $origfile \
          -temp_lapse_rate $tlapse -precip_adjustment scale -surface pdd"

# ISMIP6 Ocean model
#export pre_ocean_opts="-ocean ismip6 -ocean_ismip6_file $oceanfile"
export pre_ocean_opts="-ocean ismip6nl -ocean_ismip6nl_file $oceanfile"

# CALVING
export pre_calv_opts="-calving eigen_calving,thickness_calving \
           -front_retreat_file $OKMASK \
           -eigen_calving_K 1.0e17 -thickness_calving_threshold 100.0"

# BED
export pre_bed_opts="-bed_def none -hydrology null -bed_smoother_range 5.0e3"


###### ice physics
export basal_opts="-yield_stress mohr_coulomb -pseudo_plastic \
            -pseudo_plastic_q 0.75 -pseudo_plastic_uthreshold 100.0 \
            -till_effective_fraction_overburden 0.02"

export stress_sia_opts="-pik -sia_e 1 -sia_flow_law gpbld -age"

export stress_opts="$stress_sia_opts -stress_balance ssa+sia -ssa_method fd \
             -ssa_flow_law gpbld -ssa_e 1 $basal_opts"

###### technical
export config_opts="-config_override ${paramfile}.nc"

export tech_opts="-options_left -verbose 2 -o_order zyx -o_size big -backup_interval 3.0"
