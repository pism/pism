# Add the -config to all PISM test runs (necessary because we want to use the
# config file from *this* build, not the one that might be present in the
# install location).
file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/.petscrc "-config ${PROJECT_BINARY_DIR}/pism_config.nc")

# Copy nccmp.py into the build directory:
configure_file (${PROJECT_SOURCE_DIR}/util/nccmp.py ${PROJECT_BINARY_DIR}/nccmp.py COPYONLY)

# Set the location of PISM regression test scripts
set(PISM_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set (MEMORYCHECK_COMMAND_OPTIONS "--dsymutil=yes --trace-children=yes --trace-children-skip=*python,*nccmp.py,*rm,*which,*env,*diff,*cp,*hydra_*")

function(add_memcheck_test name binary)
  set(memcheck_command "${MEMORYCHECK_COMMAND} ${MEMORYCHECK_COMMAND_OPTIONS}")
  separate_arguments(memcheck_command)
  add_test(memcheck_${name} ${memcheck_command} ${binary} ${ARGN})
endfunction(add_memcheck_test)

function(set_memcheck_test_properties name)
  set_tests_properties(memcheck_${name} ${ARGN})
endfunction(set_memcheck_test_properties)

# A macro simplifying test creation:
macro(pism_test name script)
  add_test("${name}:${script}" ${PISM_TEST_DIR}/${script} ${PROJECT_BINARY_DIR} ${MPIEXEC} ${PROJECT_SOURCE_DIR})
  if (${Pism_TEST_USING_VALGRIND})
    add_memcheck_test("${name}:${script}" ${PISM_TEST_DIR}/${script} ${PROJECT_BINARY_DIR} ${MPIEXEC} ${PROJECT_SOURCE_DIR})
    set_memcheck_test_properties("${name}:${script}"
      PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR SUMMARY: [^0]")
  endif()

  set_property(TEST "${name}:${script}" PROPERTY
    SKIP_REGULAR_EXPRESSION
    "ModuleNotFoundError: No module named"
    "ncpdq: command not found"
    "ncap2: command not found"
    "ncks: command not found"
    "ncrcat: command not found"
  )
endmacro()

macro(pism_python_test name script)
  add_test("${name}:${script}" ${PISM_TEST_DIR}/${script} ${PROJECT_BINARY_DIR} ${MPIEXEC} ${PROJECT_SOURCE_DIR} -python ${PYTHON_EXECUTABLE})
endmacro()

# Tests themselves:
pism_test (pism_testing_tools test_regression_tools.py)

pism_test (pismr_exact_restartability_SIA_only test_01.sh)

pism_test (pismr_processor_independence test_02.sh)

pism_test (pismr_zero_length_run test_03.sh)

pism_test (pismr_regridding_during_bootstrapping test_04.sh)

pism_test (pismr_bootstrap_variable_order test_05.sh)

pism_test (pismr_bootstrap_variable_ranges test_06.sh)

pism_test (regridding:bootstrap_vs_regrid_file test_07.sh)

pism_test (regridding:vertical_direction test_08.sh)

pism_test (regridding:variable_order test_09.sh)

pism_test (regridding:processor_independence test_10.sh)

pism_test (regridding:no_enthalpy test_32.sh)

pism_test (SIA_mass_conservation test_12.sh)

pism_test (temperature_continuity_base_polythermal temp_continuity.py)

pism_test (enthalpy_symmetry_near_base test_13.sh)

pism_test (Verification:test_C test_15.sh)

pism_test (Verification:test_L test_16.sh)

pism_test (Verification:test_G test_17.sh)

pism_test (Verification:test_K test_18.sh)

pism_test (restart:i_vs_bootstrap_and_regrid_file test_23.sh)

pism_test (runtime_viewers test_27.sh)

pism_test (bootstrapping_incomplete_input test_28.sh)

pism_test (distributed_hydrology test_29.py)

pism_test (initialization_without_enthalpy test_31.sh)

pism_test (vertical_grid_expansion vertical_grid_expansion.sh)

pism_test (bed_deformation:LC:exact_restartability beddef_lc_restart.sh)

pism_test (PICO:Split-and-merge pico_split/run_test.sh)

if (Pism_USE_PROJ)
  pism_test (epsg_code_processing test_epsg_processing.py)
endif()

if(Pism_BUILD_EXTRA_EXECS)
  # These tests require special executables. They are disabled unless
  # these executables are built. This way we don't need to explain why
  # "make && make test" shows some test failures when PISM is built
  # with default settings.
  pism_test (Verification:PISMBedThermalUnit_test_K btu_regression.sh)

  pism_test (Verification:test_V_SSAFD_CFBC ssa/ssa_test_cfbc_fd.sh)

  pism_test (Verification:test_V_SSAFEM_CFBC ssa/ssa_test_cfbc_fem.sh)

  pism_test (Verification:test_I_SSAFD ssa/ssa_testi_fd.sh)

  pism_test (Verification:test_I_SSAFEM ssa/ssa_testi_fem.sh)

  pism_test (Verification:test_J_SSAFD ssa/ssa_testj_fd.sh)

  pism_test (Verification:test_J_SSAFEM ssa/ssa_testj_fem.sh)

  pism_test (Verification:SSAFEM_linear_flow ssa/ssafem_test_linear.sh)

  pism_test (Verification:SSAFEM_plug_flow ssa/ssafem_test_plug.sh)
endif()

if(Pism_BUILD_PYTHON_BINDINGS)

	pism_python_test (Python:Verification:test_I_SSAFD ssa/ssa_testi_fd.sh)

	pism_python_test (Python:Verification:test_I_SSAFEM ssa/ssa_testi_fem.sh)

	pism_python_test (Python:Verification:test_J_SSAFD ssa/ssa_testj_fd.sh)

	pism_python_test (Python:Verification:test_J_SSAFEM ssa/ssa_testj_fem.sh)

	pism_python_test (Python:Verification:SSAFEM_linear_flow ssa/ssafem_test_linear.sh)

        pism_python_test (Python:Verification:SSAFEM_plug_flow ssa/ssafem_test_plug.sh)

        pism_python_test (Python:Verification:test_V_SSAFD_CFBC ssa/ssa_test_cfbc_fd.sh)

        pism_python_test (Python:Verification:test_V_SSAFEM_CFBC ssa/ssa_test_cfbc_fem.sh)

        pism_python_test (Python:sia_forward.py test_33.sh)

# Inversion regression tests.

        execute_process (COMMAND ${PYTHON_EXECUTABLE} -c "import siple"
          RESULT_VARIABLE IMPORT_SIPLE_ERRCODE
          OUTPUT_QUIET ERROR_QUIET)

        set (Pism_HAS_SIPLE ${IMPORT_SIPLE_ERRCODE} EQUAL 0)

        # Copy over commonly used files
        foreach (FILE
            ${PROJECT_SOURCE_DIR}/examples/inverse/make_synth_ssa.py
            ${PROJECT_SOURCE_DIR}/examples/inverse/pismi.py
            ${PROJECT_SOURCE_DIR}/examples/inverse/verify_ssa_inv.py
            inverse/build_tiny.py)
          get_filename_component(OUTPUT ${FILE} NAME)
          configure_file (${FILE} ${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT} COPYONLY)
        endforeach()


        if (Pism_HAS_SIPLE) # siple is installed

          message (STATUS "Enabling regression tests that use siple")

          pism_python_test (Python:inversion:nlcg inverse/tiny_nlcg.sh)

        endif()

        pism_python_test (Python:inversion:tikhonov  inverse/tiny_tikhonov_lmvm.sh)

endif()
