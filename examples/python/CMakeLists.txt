ENABLE_TESTING()

# Set the location of PISM regression test scripts
set(PISM_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

macro(pism_python_test name script arg)
  set(test_name "Python:${name}:${script}_${arg}")
  add_test(${test_name} ${PISM_TEST_DIR}/${script} ${PISM_INSTALL_PYTHON_MODULE_DIR}
    ${PISM_TEST_DIR} ${arg})
  set_property(TEST ${test_name}
    PROPERTY FAIL_REGULAR_EXPRESSION "Traceback"
    )
endmacro()

pism_python_test("SSA_test_I" "test_wrapper.py" 1)

pism_python_test("SSA_test_J" "test_wrapper.py" 2)

pism_python_test("SSA_test_linear" "test_wrapper.py" 3)

pism_python_test("SSA_test_plug" "test_wrapper.py" 4)