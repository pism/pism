execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c "import nose"
  RESULT_VARIABLE PISM_NOSE)

if (NOT ${PISM_NOSE})
  set(Pism_USE_NOSE_TESTS TRUE
    CACHE BOOL "Enable regression tests that use nose" FORCE)
else()
  set(Pism_USE_NOSE_TESTS FALSE
    CACHE BOOL "Enable regression tests that use nose" FORCE)
endif()
mark_as_advanced(Pism_USE_NOSE_TESTS)

add_test(NAME "Config:metadata_structure"
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/config_test.py pism_config.nc
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

add_test(NAME "Config:parameters_are_sorted"
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regression/sorted_parameters.sh ${PROJECT_SOURCE_DIR}/src/pism_config.cdl
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

if (Pism_BUILD_PYTHON_BINDINGS AND Pism_USE_NOSE_TESTS)
  message(STATUS "Enabling regression tests that use 'nose'")

  macro(pism_nose_test name executable)
    string(CONCAT test_name "Python:nose:" ${name})
    add_test(NAME ${test_name}
      COMMAND ${PYTHON_EXECUTABLE} "-m" "nose" "-v" "-s" ${CMAKE_CURRENT_SOURCE_DIR}/${executable}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  endmacro()

  pism_nose_test("misc" miscellaneous.py)
  pism_nose_test("IceModelVec2T" icemodelvec2t.py)
  pism_nose_test("enthalpy:converter" enthalpy/converter.py)
  pism_nose_test("enthalpy:column" enthalpy/column.py)
  pism_nose_test("sia:bed_smoother" bed_smoother.py)
  pism_nose_test("bed_deformation:LC:restart" regression/beddef_lc_restart.py)
  pism_nose_test("ocean" regression/ocean_models.py)
  pism_nose_test("surface" regression/surface_models.py)
  pism_nose_test("atmosphere" regression/atmosphere_models.py)
  pism_nose_test("atmosphere:LTOP" regression/orographic_precipitation.py)
  pism_nose_test("Verification:bed_deformation:LC:viscous" regression/beddef_lc_viscous.py)
  pism_nose_test("Verification:bed_deformation:LC:elastic" regression/beddef_lc_elastic.py)
  pism_nose_test("Verification:bed_deformation:iso" regression/beddef_iso.py)
  pism_nose_test("Verification:bed_deformation:given" regression/beddef_given.py)
  pism_nose_test("Verification:mass_transport" mass_transport.py)
  pism_nose_test("Verification:btu" bedrock_column.py)
  pism_nose_test("Verification:blatter" blatter_verification.py)
  pism_nose_test("frontal_melt" regression/frontal_melt_models.py)
  pism_nose_test("hydrology:steady" regression/hydrology_steady_test.py)
  pism_nose_test("file-io" regression/file.py)
  pism_nose_test("grounded_cell_fraction" grounded_cell_fraction.py)
  pism_nose_test("iceberg_remover" regression/iceberg_remover.py)
else()
  message(STATUS "Python module 'nose' was not found; some regression tests will be disabled")
endif()
