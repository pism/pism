#/bin/bash
start="1978-01-01"
end="2023-01-01"
res=500

blatter_options="""
  -bp_ksp_rtol 0.001  \
  -bp_mg_coarse_ksp_type preonly  \
  -bp_mg_coarse_pc_type lu  \
  -bp_mg_levels_ksp_type richardson  \
  -bp_mg_levels_pc_type sor  \
  -bp_pc_mg_levels 3  \
  -bp_pc_type mg  \
  -bp_snes_ksp_ew 1  \
  -bp_snes_ksp_ew_version 3  \
  -bp_snes_monitor_ratio   \
  -bp_snes_rtol 0.001  \
  -stress_balance.blatter.Mz 10  \
  -stress_balance.blatter.coarsening_factor 3  \
  -stress_balance.blatter.enhancement_factor 2.0  \
  -stress_balance.blatter.flow_law gpbld  \
  -stress_balance.blatter.use_eta_transform yes  \
  -stress_balance.calving_front_stress_bc yes  \
  -stress_balance.model blatter  \
  -time_stepping.adaptive_ratio 500  \
"""

hybrid_options="""
  -stress_balance.model ssa+sia  \
  -stress_balance.ssa.method fd \
  -stress_balance.sia.enhancement_factor 2.0  \
  -stress_balance.sia.flow_law gpbld  \
  -stress_balance.sia.max_diffusivity 100000.0  \
  -stress_balance.sia.surface_gradient_method eta  \
  -stress_balance.ssa.flow_law gpbld  \
"""

ssa_options="""
  -stress_balance.model ssa  \
  -stress_balance.ssa.method fd \
  -stress_balance.sia.enhancement_factor 2.0  \
  -stress_balance.sia.flow_law gpbld  \
  -stress_balance.sia.max_diffusivity 100000.0  \
  -stress_balance.sia.surface_gradient_method eta  \
  -stress_balance.ssa.flow_law gpbld  \
"""

for sb in hybrid ssa; do
    if [ "$sb" = "hybrid" ]; then
        sb_options="$hybrid_options"
        inv="make_synth_ssa.py -generate_observed"
    else
        sb_options="$ssa_options"
        inv="make_synth_ssa.py -generate_observed"
    fi
    infile=state_${sb}_g${res}m_RGI2000-v7.0-C-01-04374_id_0_${start}_${end}_0.nc
    ofile=synth_obs_${sb}_g${res}m_RGI2000-v7.0-C-01-04374_id_0_${start}_${end}.nc

    mpirun -np 8 python $inv \
       -i $infile \
       $sb_options \
       -o $ofile
done

for sb in blatter; do
    infile=state_${sb}_g${res}m_RGI2000-v7.0-C-01-04374_id_0_${start}_${end}_0.nc
    ofile=synth_obs_${sb}_g${res}m_RGI2000-v7.0-C-01-04374_id_0_${start}_${end}.nc
    cp $infile $ofile
done


