#/bin/bash

for f in grid_g50m_RGI2000-v7.0-C-01-04374.nc bootfile_g50m_RGI2000-v7.0-C-01-04374.nc obs_RGI2000-v7.0-C-01-04374_0.nc era5_wgs84_RGI2000-v7.0-C-01-04374.nc; do
    wget -nc https://pism-cloud-data.s3.amazonaws.com/inverse/$f
done

start="1978-01-01"
end="2023-01-01"
res=500

blatter_options="""
  -bp_ksp_monitor   \
  -bp_ksp_rtol 0.001  \
  -bp_ksp_view_singularvalues   \
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

for sb in blatter hybrid ; do
    if [ "$sb" = "hybrid" ]; then
        sb_options="$hybrid_options"
    else
        sb_options="$blatter_options"
    fi
    postfix=${sb}_g${res}m_RGI2000-v7.0-C-01-04374_id_0_${start}_${end}.nc
    ofile=state_$postfix
    sfile=spatial_$postfix
    ofile_0=state_${sb}_g${res}m_RGI2000-v7.0-C-01-04374_id_0_${start}_${end}_0.nc
    mpirun -np  8 pism \
      $sb_options \
  -atmosphere.given.file era5_wgs84_RGI2000-v7.0-C-01-04374.nc  \
  -atmosphere.models given  \
  -basal_resistance.pseudo_plastic.enabled yes  \
  -basal_resistance.pseudo_plastic.q 0.75  \
  -basal_resistance.pseudo_plastic.u_threshold 100m/yr  \
  -basal_yield_stress.model mohr_coulomb  \
  -basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden 0.025  \
  -basal_yield_stress.mohr_coulomb.till_phi_default 30  \
  -calving.methods float_kill \
  -energy.model enthalpy  \
  -geometry.front_retreat.use_cfl yes  \
  -geometry.part_grid.enabled yes  \
  -geometry.remove_icebergs yes  \
  -grid.Lbz 0  \
  -grid.Lz 2000  \
  -grid.Mbz 1  \
  -grid.Mz 101  \
  -grid.dx ${res}m  \
  -grid.dy ${res}m  \
  -grid.file grid_g50m_RGI2000-v7.0-C-01-04374.nc  \
  -stress_balance.sia.bed_smoother.range  $res \
  -grid.registration center  \
  -hydrology.model null \
  -hydrology.null_diffuse_till_water \
  -input.bootstrap yes  \
  -input.file bootfile_g50m_RGI2000-v7.0-C-01-04374.nc  \
  -input.forcing.buffer_size 390  \
  -input.forcing.time_extrapolation yes  \
  -output.file $ofile \
  -output.spatial.file $sfile \
  -output.spatial.vars velsurf_mag,velbase_mag,usurf,thk \
  -output.spatial.times yearly \
  -output.size medium  \
  -output.sizes.medium sftgif,velsurf_mag,mask,usurf,velbase_mag,velsurf  \
  -surface.debm_simple.c1 30  \
  -surface.debm_simple.c2 -120  \
  -surface.debm_simple.interpret_precip_as_snow no  \
  -surface.force_to_thickness.file bootfile_g50m_RGI2000-v7.0-C-01-04374.nc  \
  -surface.models debm_simple  \
  -time.calendar standard  \
  -time.end $end  \
  -time.reference_date $start  \
  -time.start $start  \
  -time_stepping.skip.enabled yes  \
  -time_stepping.skip.max 100

    cdo setmisstoc,0 $ofile $ofile_0
    ncatted -a _FillValue,uvelsurf,d,, -a _FillValue,vvelsurf,d,, \
            -a _FillValue,u_ssa,d,, -a _FillValue,v_ssa,d,, $ofile_0
    ncks -A -v pism_config $ofile $ofile_0
    
done
