# Copyright (C) 2009-2016  PISM authors

# grids
export THIRTYKMGRID="-Mx 200 -My 200 -Lz 5000 -Lbz 2000 -Mz 41 -Mbz 16"
export TWENTYKMGRID="-Mx 300 -My 300 -Lz 5000 -Lbz 2000 -Mz 81 -Mbz 21"
export FIFTEENKMGRID="-Mx 400 -My 400 -Lz 5000 -Lbz 2000 -Mz 81 -Mbz 21"
export TWELVEKMGRID="-Mx 500 -My 500 -Lz 5000 -Lbz 2000 -Mz 101 -Mbz 31"
export TENKMGRID="-Mx 600 -My 600 -Lz 5000 -Lbz 2000 -Mz 101 -Mbz 31"
export SEVENKMGRID="-Mx 800 -My 800 -Lz 5000 -Lbz 2000 -Mz 151 -Mbz 31"
export FIVEKMGRID="-Mx 1200 -My 1200 -Lz 5000 -Lbz 2000 -Mz 201 -Mbz 51"

# skips
export SKIPTHIRTYKM=10
export SKIPTWENTYKM=10
export SKIPFIFTEENKM=10
export SKIPTWELVEKM=50
export SKIPTENKM=100
export SKIPSEVENKM=100
export SKIPFIVEKM=200

#PIK-stuff; notes:
# 1)   '-pik' = '-cfbc -part_grid -kill_icebergs -subgl'
# 2)   -meltfactor_pik 5e-3 is default when using -ocean pik
export PIKPHYS="-ssa_method fd -ssa_e 0.6 -pik -calving eigen_calving,thickness_calving -eigen_calving_K 2.0e18 -thickness_calving_threshold 200.0"
export PIKPHYS_COUPLING="-atmosphere given -atmosphere_given_file $PISM_INDATANAME -surface simple -ocean pik -meltfactor_pik 5e-3"

# dynamics related options
export SIA_ENHANCEMENT="-sia_e 3.0"
export PARAMS="-pseudo_plastic -pseudo_plastic_q 0.25 -till_effective_fraction_overburden 0.02 -tauc_slippery_grounding_lines"
export TILLPHI="-topg_to_phi 15.0,40.0,-300.0,700.0"
export FULLPHYS="-stress_balance ssa+sia -hydrology null $PARAMS $TILLPHI"

# perhaps use these if KSP "diverged" errors occur
export STRONGKSP="-ssafd_ksp_type gmres -ssafd_ksp_norm_type unpreconditioned -ssafd_ksp_pc_side right -ssafd_pc_type asm -ssafd_sub_pc_type lu"
