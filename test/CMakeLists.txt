find_program(NOSE_EXECUTABLE "nosetests")
mark_as_advanced(NOSE_EXECUTABLE)

if (Pism_BUILD_PYTHON_BINDINGS AND NOSE_EXECUTABLE)
  message(STATUS "Enabling PISM Python tests that use nose")
  add_test(NAME "nosetests_py"
    COMMAND ${NOSE_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/nosetests.py)
endif()