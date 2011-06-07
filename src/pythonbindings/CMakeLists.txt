# A patch by Tristan Carel that will be included in the next release of CMake
# Downloaded from http://www.itk.org/Bug/view.php?id=4147&nbn=4

#
# Get dependencies of the generated wrapper.
#
MACRO(SWIG_GET_WRAPPER_DEPENDENCIES swigFile genWrapper language DEST_VARIABLE)
  GET_FILENAME_COMPONENT(swig_getdeps_basename ${swigFile} NAME_WE)
  GET_FILENAME_COMPONENT(swig_getdeps_outdir ${genWrapper} PATH)
  GET_SOURCE_FILE_PROPERTY(swig_getdeps_extra_flags "${swigFile}" SWIG_FLAGS)
  IF("${swig_getdeps_extra_flags}" STREQUAL "NOTFOUND")
    SET(swig_getdeps_extra_flags "")
  ENDIF("${swig_getdeps_extra_flags}" STREQUAL "NOTFOUND")

  IF(NOT swig_getdeps_outdir)
    SET(swig_getdeps_outdir ${CMAKE_CURRENT_BINARY_DIR})
  ENDIF(NOT swig_getdeps_outdir)
  SET(swig_getdeps_depsfile
    ${swig_getdeps_outdir}/swig_${swig_getdeps_basename}_deps.txt)
  GET_DIRECTORY_PROPERTY(swig_getdeps_include_directories INCLUDE_DIRECTORIES)
  SET(swig_getdeps_include_dirs)
  FOREACH(it ${swig_getdeps_include_directories})
    SET(swig_getdeps_include_dirs ${swig_getdeps_include_dirs} "-I${it}")
  ENDFOREACH(it)
  EXECUTE_PROCESS(
    COMMAND ${SWIG_EXECUTABLE}
    -MM -MF ${swig_getdeps_depsfile} ${swig_getdeps_extra_flags}
    ${CMAKE_SWIG_FLAGS} -${language}
    -o ${genWrapper} ${swig_getdeps_include_dirs} ${swigFile}
    RESULT_VARIABLE swig_getdeps_result
    ERROR_VARIABLE swig_getdeps_error
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  IF(NOT ${swig_getdeps_error} EQUAL 0)
    MESSAGE(SEND_ERROR "Command \"${SWIG_EXECUTABLE} -MM -MF ${swig_getdeps_depsfile} ${swig_getdeps_extra_flags} ${CMAKE_SWIG_FLAGS} -${language} -o ${genWrapper} ${swig_getdeps_include_dirs} ${swigFile}\" failed with output:\n${swig_getdeps_error}")
    SET(swig_getdeps_dependencies "")
  ELSE(NOT ${swig_getdeps_error} EQUAL 0)
    FILE(READ ${swig_getdeps_depsfile} ${DEST_VARIABLE})
    # Remove the first line
    STRING(REGEX REPLACE "^.+: +\\\\\n +" ""
      ${DEST_VARIABLE} "${${DEST_VARIABLE}}")
    # Clean the end of each line
    STRING(REGEX REPLACE " +(\\\\)?\n" "\n" ${DEST_VARIABLE}
      "${${DEST_VARIABLE}}")
    # Clean beginning of each line
    STRING(REGEX REPLACE "\n +" "\n"
      ${DEST_VARIABLE} "${${DEST_VARIABLE}}")
    # clean paths
    STRING(REGEX REPLACE "\\\\\\\\" "/" ${DEST_VARIABLE}
      "${${DEST_VARIABLE}}")
    STRING(REGEX REPLACE "\n" ";"
      ${DEST_VARIABLE} "${${DEST_VARIABLE}}")
  ENDIF(NOT ${swig_getdeps_error} EQUAL 0)
ENDMACRO(SWIG_GET_WRAPPER_DEPENDENCIES)

if (Pism_BUILD_PYTHON_BINDINGS)
    include(${SWIG_USE_FILE})

    include_directories (../verif/tests/ . ${PETSC4PY_INCLUDES} ${PYTHON_INCLUDES})
    
    set(PISM_PYTHON_BINDING_LIBRARIES "pismbase" "pismutil" "pismflowlaws" "pismverif" "pismstressbalance" )

    # Get module dependencies:
    SWIG_GET_WRAPPER_DEPENDENCIES("${CMAKE_CURRENT_SOURCE_DIR}/PISM.i" "." python PISM_Python_deps)

    SET(CMAKE_SWIG_FLAGS "")
    SET(SWIG_MODULE_NAME "cpp")

    set_source_files_properties(PISM.i PROPERTIES
      CPLUSPLUS ON
      SWIG_MODULE_NAME cpp)

    set(SWIG_MODULE_cpp_EXTRA_DEPS ${PISM_Python_deps}) # has to be set *before* calling SWIG_ADD_MODULE

    SWIG_ADD_MODULE(${SWIG_MODULE_NAME} python PISM.i pism_python.cc pism_python_signal.cpp)
    SWIG_LINK_LIBRARIES(${SWIG_MODULE_NAME}
      ${PYTHON_LIBRARY} ${PISM_PYTHON_BINDING_LIBRARIES} )

    # Disable compiler warnings:
    SET_SOURCE_FILES_PROPERTIES(${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-w")

    get_target_property(SWIG_MODULE_LOCATION ${SWIG_MODULE_cpp_REAL_NAME} LOCATION)

    # Install Python .py files
    install(FILES
      ${SWIG_MODULE_LOCATION} ${CMAKE_CURRENT_BINARY_DIR}/${SWIG_MODULE_NAME}.py
      DESTINATION ${PISM_INSTALL_PYTHON_MODULE_DIR}/PISM
      COMPONENT RuntimeLibraries
    )

    install(DIRECTORY
      ${CMAKE_SOURCE_DIR}/site-packages/PISM
      DESTINATION ${PISM_INSTALL_PYTHON_MODULE_DIR}
      USE_SOURCE_PERMISSIONS
      COMPONENT RuntimeLibraries
      PATTERN ".svn" EXCLUDE
    )

  install(CODE "MESSAGE(\"\n   ATTENTION! To complete installation of Python bindings you need to add\n\n   export PYTHONPATH=${PISM_INSTALL_PYTHON_MODULE_DIR}:\$PYTHONPATH\n\n   to your .bashrc or equivalent. \")")
    
endif(Pism_BUILD_PYTHON_BINDINGS)
