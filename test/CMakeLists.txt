find_package(Python3 COMPONENTS Interpreter)

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
  macro(pism_add_unittest name executable)
    string(CONCAT test_name "Python:unittest:" ${name})
    file(REAL_PATH ${executable} exec)
    cmake_path(GET exec PARENT_PATH workdir)
    cmake_path(GET exec FILENAME filename)
    # message(STATUS "Adding test ${filename} in ${workdir}")
    add_test(NAME ${test_name}
      COMMAND ${Python3_EXECUTABLE} "-m" "unittest" "-v" ${filename}
      WORKING_DIRECTORY ${workdir})
    set_tests_properties(${test_name} PROPERTIES
      ENVIRONMENT "PETSC_OPTIONS=-config \"${PROJECT_BINARY_DIR}/pism_config.nc\";PYTHONPATH=$ENV{PYTHONPATH}:${PROJECT_BINARY_DIR}/site-packages"
      FIXTURES_REQUIRED config)
  endmacro()

  pism_add_unittest("misc" miscellaneous.py)
  pism_add_unittest("array:Forcing" array_forcing.py)
  pism_add_unittest("Verification:btu" bedrock_column.py)
  pism_add_unittest("atmosphere" regression/atmosphere_models.py)
  pism_add_unittest("surface" regression/surface_models.py)
  pism_add_unittest("ocean" regression/ocean_models.py)
  pism_add_unittest("atmosphere:LTOP" regression/orographic_precipitation.py)
  pism_add_unittest("frontal_melt" regression/frontal_melt_models.py)
  pism_add_unittest("iceberg_remover" regression/iceberg_remover.py)
  pism_add_unittest("Verification:blatter" blatter_verification.py)
  pism_add_unittest("grounded_cell_fraction" grounded_cell_fraction.py)
  pism_add_unittest("sia:bed_smoother" bed_smoother.py)
  pism_add_unittest("enthalpy:converter" enthalpy/converter.py)
  pism_add_unittest("enthalpy:column" enthalpy/column.py)
  pism_add_unittest("file-io" regression/file.py)
  pism_add_unittest("hydrology:steady" regression/hydrology_steady_test.py)
  pism_add_unittest("bed_deformation:LC:restart" regression/beddef_lc_restart.py)
  pism_add_unittest("Verification:mass_transport" mass_transport.py)
  pism_add_unittest("Verification:bed_deformation:LC:viscous" regression/beddef_lc_viscous.py)
  pism_add_unittest("Verification:bed_deformation:LC:elastic" regression/beddef_lc_elastic.py)
  pism_add_unittest("Verification:bed_deformation:iso" regression/beddef_iso.py)
  pism_add_unittest("Verification:bed_deformation:given" regression/beddef_given.py)
endif()
