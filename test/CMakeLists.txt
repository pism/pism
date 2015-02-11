find_program(NOSE_EXECUTABLE "nosetests")
mark_as_advanced(NOSE_EXECUTABLE)

if (Pism_BUILD_PYTHON_BINDINGS AND NOSE_EXECUTABLE)
  message(STATUS "Enabling PISM Python tests that use nose")
  add_test(NAME "Python:nose:misc"
    COMMAND ${NOSE_EXECUTABLE} "-v" "-s" ${CMAKE_CURRENT_SOURCE_DIR}/nosetests.py)
  add_test(NAME "Python:nose:EnthalpyConverter"
    COMMAND ${NOSE_EXECUTABLE} "-v" "-s" ${CMAKE_CURRENT_SOURCE_DIR}/enthalpy_converter.py)
endif()