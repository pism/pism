if (Pism_BUILD_PYTHON_BINDINGS)

  if (NOT (${CMAKE_VERSION} VERSION_LESS "3.13"))
    # Policy CMP0078 was introduced in 3.13. See "cmake --help-policies" for details.
    cmake_policy(SET CMP0078 NEW)       #controls names of SWIG targets
  endif()

  include(${SWIG_USE_FILE})

  include_directories (. ${PETSC4PY_INCLUDES} ${PYTHON_INCLUDE_DIRS} ${Pism_SOURCE_DIR}/src)

  set(PISM_PYTHON_BINDING_LIBRARIES "pism")

  list(APPEND PISM_Python_deps
    AgeModel.i
    geometry.i
    pism_BedDef.i
    pism_ColumnSystem.i
    pism_DM.i
    pism_FlowLaw.i
    pism_Hydrology.i
    pism_IceGrid.i
    pism_IceModelVec.i
    pism_File.i
    pism_SIA.i
    pism_SSA.i
    pism_blatter.i
    pism_VariableMetadata.i
    pism_Vars.i
    pism_Vec.i
    pism_atmosphere.i
    pism_energy.i
    pism_exception.i
    pism_inverse.i
    pism_ocean.i
    pism_options.i
    pism_surface.i
    pism_verification.i
    pism_frontalmelt.i
    pism_calving.i
    )

  set_property(SOURCE PISM.i PROPERTY SWIG_FLAGS "-I${CMAKE_BINARY_DIR}" APPEND)

  if(Pism_USE_JANSSON)
    set_property(SOURCE PISM.i PROPERTY SWIG_FLAGS "-DPism_USE_JANSSON=1" APPEND)
  else()
    set_property(SOURCE PISM.i PROPERTY SWIG_FLAGS "-DPism_USE_JANSSON=0" APPEND)
  endif()

  set_source_files_properties(PISM.i PROPERTIES
    CPLUSPLUS ON
    SWIG_MODULE_NAME cpp)

  # has to be set *before* calling SWIG_ADD_MODULE
  set(SWIG_MODULE_cpp_EXTRA_DEPS ${PISM_Python_deps})

  if (${CMAKE_VERSION} VERSION_LESS "3.8.0")
    SWIG_ADD_MODULE(cpp python PISM.i pism_python.cc)
  else()
    SWIG_ADD_LIBRARY(cpp
      LANGUAGE python
      SOURCES PISM.i pism_python.cc)
  endif()
  SWIG_LINK_LIBRARIES(cpp ${PYTHON_LIBRARIES} ${PISM_PYTHON_BINDING_LIBRARIES})

  # Disable compiler warnings:
  SET_SOURCE_FILES_PROPERTIES(${swig_generated_file_fullname}
    PROPERTIES COMPILE_FLAGS "-w")

  add_custom_target (clean_swig
    COMMAND ${CMAKE_COMMAND} -E touch PISM.i
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM)

  # Disable compiler warnings:
  SET_SOURCE_FILES_PROPERTIES(${swig_generated_file_fullname}
    PROPERTIES COMPILE_FLAGS "-w")

  # Generate version_info.py to make PETSc's version information available to PISM's
  # Python bindings. This is necessary to initialize petsc4py using the correct PETSC_ARCH
  # value *before* importing PISM.cpp.
  configure_file(version_info.py.in version_info.py ESCAPE_QUOTES)

  if(${CMAKE_VERSION} VERSION_LESS "3.13")
    set(TARGET_FILE $<TARGET_FILE:_cpp>)
    set(CPP_REAL_NAME ${SWIG_MODULE_cpp_REAL_NAME})
  else()
    set(TARGET_FILE $<TARGET_FILE:cpp>)
    set(CPP_REAL_NAME cpp)
  endif()

  # Install Python .py files
  install(FILES
    ${TARGET_FILE}
    ${CMAKE_CURRENT_BINARY_DIR}/cpp.py
    ${CMAKE_CURRENT_BINARY_DIR}/version_info.py
    DESTINATION ${PISM_INSTALL_PYTHON_MODULE_DIR}/PISM
    COMPONENT RuntimeLibraries
    )

  add_custom_target (install_python_bindings_in_the_build_directory ALL
    COMMAND ${CMAKE_COMMAND} -E copy ${TARGET_FILE} ${CMAKE_BINARY_DIR}/site-packages/PISM
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/cpp.py ${CMAKE_BINARY_DIR}/site-packages/PISM/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/version_info.py ${CMAKE_BINARY_DIR}/site-packages/PISM/
    VERBATIM
    DEPENDS ${CPP_REAL_NAME} install_site_packages_in_the_build_directory )

  # The Python code snippet below returns 1 if Python can import a PISM
  # and if the directory of said PISM matches the place we are about to
  # install PISM's Python wrappers. If the return value is not 1, we
  # print a message stating how to set your PYTHONPATH.
  set (pism_test_py "${PROJECT_BINARY_DIR}/pism_test.py")
  file (WRITE "${pism_test_py}" "
from os.path import realpath, join
try:
    import PISM
except:
    exit(0)

if realpath(PISM.__file__) == realpath(join('${PISM_INSTALL_PYTHON_MODULE_DIR}', 'PISM', '__init__.pyc')):
    exit(1)
")

  install(CODE "EXEC_PROGRAM( ${PYTHON_EXECUTABLE} ARGS ${pism_test_py} RETURN_VALUE PISM_PYTHON_MATCHES OUTPUT_VARIABLE PYTHON_OUT )
if(NOT PISM_PYTHON_MATCHES)
    MESSAGE(\"\n   ATTENTION! To complete installation of Python bindings you need to add\n\n   export PYTHONPATH=${PISM_INSTALL_PYTHON_MODULE_DIR}:\$PYTHONPATH\n\n   to your .bashrc or equivalent. \")
endif()
" )

endif(Pism_BUILD_PYTHON_BINDINGS)
