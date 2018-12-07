if (Pism_BUILD_PYTHON_BINDINGS)
    include(${SWIG_USE_FILE})

    include_directories (. ${PETSC4PY_INCLUDES} ${PYTHON_INCLUDES} ${Pism_SOURCE_DIR}/src)

    set(PISM_PYTHON_BINDING_LIBRARIES "pism")

    list(APPEND PISM_Python_deps
      geometry.i
      AgeModel.i
      EnergyModel.i
      petsc_version.i
      pism_BedDef.i
      pism_ColumnSystem.i
      pism_DM.i
      pism_FlowLaw.i
      pism_Hydrology.i
      pism_IceGrid.i
      pism_IceModelVec.i
      pism_PIO.i
      pism_SIA.i
      pism_SSA.i
      pism_Timeseries.i
      pism_VariableMetadata.i
      pism_Vars.i
      pism_Vec.i
      pism_exception.i
      pism_inverse.i
      pism_ocean.i
      pism_surface.i
      pism_atmosphere.i
      pism_options.i
      pism_verification.i
      )

    if(Pism_USE_JANSSON)
      set_property(SOURCE PISM.i PROPERTY SWIG_FLAGS "-DPISM_USE_JANSSON" APPEND)
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
    SWIG_LINK_LIBRARIES(cpp ${PYTHON_LIBRARY} ${PISM_PYTHON_BINDING_LIBRARIES})

    # Disable compiler warnings:
    SET_SOURCE_FILES_PROPERTIES(${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-w")

    add_custom_target (clean_swig
      COMMAND ${CMAKE_COMMAND} -E touch PISM.i
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      VERBATIM)

    if (${CMAKE_VERSION} VERSION_LESS "3.8.0")
      SWIG_ADD_MODULE(petsc_version python petsc_version.i petsc_version.cc)
    else()
      SWIG_ADD_LIBRARY(petsc_version
        LANGUAGE python
        SOURCES petsc_version.i petsc_version.cc)
    endif()
    SWIG_LINK_LIBRARIES(petsc_version ${PYTHON_LIBRARY})

    # Disable compiler warnings:
    SET_SOURCE_FILES_PROPERTIES(${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-w")

    # Install Python .py files
    install(FILES
      $<TARGET_FILE:_cpp>
      $<TARGET_FILE:_petsc_version>
      ${CMAKE_CURRENT_BINARY_DIR}/cpp.py
      ${CMAKE_CURRENT_BINARY_DIR}/petsc_version.py
      DESTINATION ${PISM_INSTALL_PYTHON_MODULE_DIR}/PISM
      COMPONENT RuntimeLibraries
    )

add_custom_target (install_python_bindings_in_the_build_directory ALL
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_cpp> ${CMAKE_BINARY_DIR}/site-packages/PISM
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_petsc_version> ${CMAKE_BINARY_DIR}/site-packages/PISM/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/cpp.py ${CMAKE_BINARY_DIR}/site-packages/PISM/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/petsc_version.py ${CMAKE_BINARY_DIR}/site-packages/PISM/
  VERBATIM
  DEPENDS ${SWIG_MODULE_cpp_REAL_NAME}  _petsc_version install_site_packages_in_the_build_directory )

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
