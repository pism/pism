# Add the -config to all PISM test runs (necessary because we want to use the
# config file from *this* build, not the one that might be present in the
# install location).
file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/.petscrc "-config ${PROJECT_BINARY_DIR}/pism_config.nc")

# Copy nccmp.py into the build directory:
configure_file (${PROJECT_SOURCE_DIR}/util/nccmp.py ${PROJECT_BINARY_DIR}/nccmp.py COPYONLY)

# Set the location of PISM regression test scripts
set(PISM_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# A macro simplifying test creation:
macro(pism_test name script)
  add_test(NAME "${name}"
    COMMAND ${PISM_TEST_DIR}/${script} ${PROJECT_BINARY_DIR} ${MPIEXEC} ${PROJECT_SOURCE_DIR})
endmacro()

# Tests themselves:
pism_test (pism_testing_tools test_00.sh)

pism_test (pismr_exact_restartability_[SIA_only] test_01.sh)

pism_test (pismv_processor_independence test_02.sh)

pism_test (pismr_zero_length_run test_03.sh)

pism_test (pismr_regridding_during_bootstrapping test_04.sh)

pism_test (pismr_bootstrap_variable_order test_05.sh)

pism_test (pismr_bootstrap_variable_ranges test_06.sh)

pism_test (regridding_via_-boot_file_and_-regrid_file test_07.sh)

pism_test (regridding_in_vertical_direction test_08.sh)

pism_test (pismr_regridding_variable_order test_09.sh)

pism_test (pismr_regridding_processor_independence test_10.sh)

pism_test (automatic_vertical_grid_extension test_11.sh)

pism_test (SIA_mass_conservation test_12.sh)

pism_test (temperature_continuity_at_ice-bed_interface_[cold_case] test_13.sh)

pism_test (enthalpy_symmetry_near_the_base_[pisms_-no_cold] test_14.sh)

pism_test (verif_test_E_regression:_isothermal_SIA_with_sliding test_15.sh)

pism_test (verif_test_C_regression:_isothermal_SIA_w._time-dependent_SMB test_16.sh)

pism_test (verif_test_L_regression:_isothermal_SIA_with_non-flat_bed test_17.sh)

pism_test (verif_test_G_regression:_thermo_SIA_w._time-dependent_SMB test_18.sh)

pism_test (verif_test_I_SSAFD_regression:_nonlinear_SSA_with_plastic_bed ssa/ssa_testi_fd.sh)

pism_test (verif_test_I_SSAFD_PIK_regression:_nonlinear_SSA_with_plastic_bed ssa/ssa_testi_fd_pik.sh)

pism_test (verif_test_I_SSAFEM_regression:_nonlinear_SSA_with_plastic_bed ssa/ssa_testi_fem.sh)

pism_test (verif_test_J_SSAFD_regression:_linearized_SSA,_floating ssa/ssa_testj_fd.sh)

pism_test (verif_test_J_SSAFD_PIK_regression:_linearized_SSA,_floating ssa/ssa_testj_fd_pik.sh)

pism_test (verif_test_J_SSAFEM_regression:_linearized_SSA,_floating ssa/ssa_testj_fem.sh)

pism_test (verif_SSAFEM_linear_flow_regression ssa/ssafem_test_linear.sh)

pism_test (verif_SSAFEM_plug_flow_regression ssa/ssafem_test_plug.sh)

pism_test (verif_test_K_regression:_cold_ice_method,_bedrock_thermal_layer test_21.sh)

pism_test (EISMINT-Ross_experiment_regression test_22.sh)

pism_test (GPBLD_flow_law_regression test_24.sh)

pism_test (Paterson-Budd_flow_law_regression test_25.sh)

pism_test (Schoof_[2003]_bed_roughness_SIA_parameterization_regression test_26.sh)

pism_test (comparing_restart:_-i_vs_-boot_file_and_-regrid_file test_27.sh)

pism_test (Hooke_flow_law_regression test_28.sh)

pism_test (hybrid_[GK-PB]_flow_law_regression test_29.sh)

pism_test (Cold_part_of_Paterson-Budd_flow_law_regression test_30.sh)

pism_test (runtime_viewers test_31.sh)

add_custom_target (test
  COMMAND ${CMAKE_CTEST_COMMAND}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})