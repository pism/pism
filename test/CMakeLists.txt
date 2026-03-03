find_package(Python3 COMPONENTS Interpreter)

execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import nose"
  RESULT_VARIABLE PISM_NOSE
  OUTPUT_QUIET ERROR_QUIET
)

if (NOT ${PISM_NOSE})
  set(Pism_USE_NOSE_TESTS TRUE
    CACHE BOOL "Enable regression tests that use nose" FORCE)
else()
  set(Pism_USE_NOSE_TESTS FALSE
    CACHE BOOL "Enable regression tests that use nose" FORCE)
endif()
mark_as_advanced(Pism_USE_NOSE_TESTS)

# add a "test" that generates the config file to ensure that it is available later
add_test(NAME "Config:generate"
  COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target pism_config
)
set_tests_properties("Config:generate" PROPERTIES FIXTURES_SETUP config)

add_test(NAME "Config:metadata_structure"
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/config_test.py pism_config.nc
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
set_tests_properties("Config:metadata_structure" PROPERTIES
  FIXTURES_REQUIRED config)

add_test(NAME "Config:parameters_are_sorted"
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regression/sorted_parameters.sh ${PROJECT_SOURCE_DIR}/src/pism_config.cdl
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
set_tests_properties("Config:parameters_are_sorted" PROPERTIES
  FIXTURES_REQUIRED config)

if (Pism_BUILD_PYTHON_BINDINGS)

  macro(pism_unittest_test name executable)
    string(CONCAT test_name "Python:unittest:" ${name})
    file(REAL_PATH ${executable} exec)
    cmake_path(GET exec PARENT_PATH workdir)
    cmake_path(GET exec FILENAME filename)
    # message(STATUS "Adding test ${filename} in ${workdir}")
    add_test(NAME ${test_name}
      COMMAND ${Python3_EXECUTABLE} "-m" "unittest" "-v" ${filename}
      WORKING_DIRECTORY ${workdir})
    set_tests_properties(${test_name} PROPERTIES
      ENVIRONMENT "PETSC_OPTIONS=-config \"${PROJECT_BINARY_DIR}/pism_config.nc\""
      FIXTURES_REQUIRED config)
  endmacro()

  pism_unittest_test("misc" miscellaneous.py)
  pism_unittest_test("array:Forcing" icemodelvec2t.py)
  pism_unittest_test("Verification:btu" bedrock_column.py)
  pism_unittest_test("atmosphere" regression/atmosphere_models.py)
  pism_unittest_test("surface" regression/surface_models.py)
  pism_unittest_test("ocean" regression/ocean_models.py)
  pism_unittest_test("atmosphere:LTOP" regression/orographic_precipitation.py)
  pism_unittest_test("frontal_melt" regression/frontal_melt_models.py)
  pism_unittest_test("iceberg_remover" regression/iceberg_remover.py)
  pism_unittest_test("Verification:blatter" blatter_verification.py)
  pism_unittest_test("grounded_cell_fraction" grounded_cell_fraction.py)
  pism_unittest_test("sia:bed_smoother" bed_smoother.py)
  pism_unittest_test("enthalpy:converter" enthalpy/converter.py)
  pism_unittest_test("enthalpy:column" enthalpy/column.py)
  pism_unittest_test("file-io" regression/file.py)
  pism_unittest_test("hydrology:steady" regression/hydrology_steady_test.py)
  pism_unittest_test("bed_deformation:LC:restart" regression/beddef_lc_restart.py)
  pism_unittest_test("Verification:mass_transport" mass_transport.py)
  pism_unittest_test("Verification:bed_deformation:LC:viscous" regression/beddef_lc_viscous.py)
endif()

if (Pism_BUILD_PYTHON_BINDINGS AND Pism_USE_NOSE_TESTS)
  message(STATUS "Enabling regression tests that use 'nose'")

  macro(pism_nose_test name executable)
    string(CONCAT test_name "Python:nose:" ${name})
    add_test(NAME ${test_name}
      COMMAND ${Python3_EXECUTABLE} "-m" "nose" "-v" "-s" ${CMAKE_CURRENT_SOURCE_DIR}/${executable}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  endmacro()

  pism_nose_test("Verification:bed_deformation:LC:elastic" regression/beddef_lc_elastic.py)
  pism_nose_test("Verification:bed_deformation:iso" regression/beddef_iso.py)
  pism_nose_test("Verification:bed_deformation:given" regression/beddef_given.py)
else()
  message(STATUS "Python module 'nose' was not found; some regression tests will be disabled")
endif()
