find_program(NOSE_EXECUTABLE
  NAMES nosetests nosetests3)
mark_as_advanced(NOSE_EXECUTABLE)

if (Pism_DEBUG)
  add_test(NAME "Config:metadata_structure"
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/config_test.py pism_config.nc
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})

  add_test(NAME "Config:parameters_are_sorted"
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regression/sorted_parameters.sh ${PROJECT_SOURCE_DIR}/src/pism_config.cdl
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endif ()

if (Pism_BUILD_PYTHON_BINDINGS AND NOSE_EXECUTABLE)
  message(STATUS "Enabling PISM Python tests that use nose")

  macro(pism_nose_test name executable)
    add_test(NAME ${name}
      COMMAND ${NOSE_EXECUTABLE} "-v" "-s" ${CMAKE_CURRENT_SOURCE_DIR}/${executable}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  endmacro()

  pism_nose_test("Python:nose:misc" miscellaneous.py)
  pism_nose_test("Python:nose:IceModelVec2T" icemodelvec2t.py)
  pism_nose_test("Python:nose:enthalpy:converter" enthalpy/converter.py)
  pism_nose_test("Python:nose:enthalpy:column" enthalpy/column.py)
  pism_nose_test("Python:nose:sia:bed_smoother" bed_smoother.py)
  pism_nose_test("Python:nose:bed_deformation:LC:restart" regression/beddef_lc_restart.py)
  pism_nose_test("Python:nose:ocean" regression/ocean_models.py)
  pism_nose_test("Python:nose:surface" regression/surface_models.py)
  pism_nose_test("Python:nose:atmosphere" regression/atmosphere_models.py)
  pism_nose_test("Python:nose:atmosphere:LTOP" regression/orographic_precipitation.py)
  pism_nose_test("Python:Verification:nose:bed_deformation:LC:viscous" regression/beddef_lc_viscous.py)
  pism_nose_test("Python:Verification:nose:bed_deformation:LC:elastic" regression/beddef_lc_elastic.py)
  pism_nose_test("Python:Verification:nose:bed_deformation:iso" regression/beddef_iso.py)
  pism_nose_test("Python:Verification:nose:bed_deformation:given" regression/beddef_given.py)
  pism_nose_test("Python:Verification:nose:mass_transport" mass_transport.py)
  pism_nose_test("Python:Verification:nose:btu" bedrock_column.py)
  pism_nose_test("Python:nose:frontal_melt" regression/frontal_melt_models.py)
  pism_nose_test("Python:nose:hydrology:steady" regression/hydrology_steady_test.py)
  pism_nose_test("Python:nose:file-io" regression/file.py)
  pism_nose_test("Python:nose:grounded_cell_fraction" grounded_cell_fraction.py)
else()
  message(STATUS "nose was not found; some regression tests will be disabled")
endif()
