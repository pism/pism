find_program(NOSE_EXECUTABLE "nosetests")
mark_as_advanced(NOSE_EXECUTABLE)

if (Pism_DEBUG)
  add_test(NAME "Config:metadata_structure"
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/config_test.py pism_config.nc
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
  pism_nose_test("Python:nose:bed_deformation:LC:restart" beddef_lc_restart.py)
  pism_nose_test("Python:nose:ocean" regression/ocean_models.py)
  pism_nose_test("Python:nose:surface" regression/surface_models.py)
  pism_nose_test("Python:nose:atmosphere" regression/atmosphere_models.py)
  pism_nose_test("Python:nose:atmosphere:LTOP" regression/orographic_precipitation.py)
  pism_nose_test("Python:Verification:nose:bed_deformation:LC" beddef_lc.py)
  pism_nose_test("Python:Verification:nose:bed_deformation:iso" beddef_iso.py)
  pism_nose_test("Python:Verification:nose:mass_transport" mass_transport.py)
  pism_nose_test("Python:Verification:nose:btu" bedrock_column.py)
endif()
