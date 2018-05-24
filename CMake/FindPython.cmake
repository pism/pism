# - Find python executable and libraries.
#
#  PYTHON_EXECUTABLE   - name of the python executable.
#  PYTHON_INCLUDES     - where to find Python.h, etc.
#  PYTHON_FOUND        - True if python is found

if(PYTHON_EXECUTABLE AND PYTHON_INCLUDES AND PYTHON_LIBRARY )
    set(PYTHON_FIND_QUIETLY TRUE)
endif()

find_program(PYTHON_EXECUTABLE python DOC "python interpreter")

if(PYTHON_EXECUTABLE)
    execute_process( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"
                     OUTPUT_VARIABLE PYTHON_INCLUDES
                     RESULT_VARIABLE PYTHON_INCLUDES_NOT_FOUND
                     OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBPL'))"
                     OUTPUT_VARIABLE PYTHON_LIBDIR
                     RESULT_VARIABLE PYTHON_LIBDIR_NOT_FOUND
                     OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_version())"
                     OUTPUT_VARIABLE PYTHON_VERSION
                     RESULT_VARIABLE PYTHON_VERSION_NOT_FOUND
                     OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

if(PYTHON_LIBDIR)
    find_library( PYTHON_LIBRARY "python${PYTHON_VERSION}" HINTS ${PYTHON_LIBDIR} )
endif()

if(PYTHON_EXECUTABLE)
  if(NOT PYTHON_FIND_QUIETLY)
    message( STATUS "Found Python executable: ${PYTHON_EXECUTABLE}")
  endif()
else()
  if(FIND_PYTHON_REQUIRED)
    message( FATAL_ERROR "Python executable missing")
  endif()
endif()

if(PYTHON_INCLUDES)
  if(NOT PYTHON_FIND_QUIETLY)
    message( STATUS "Found Python includes: ${PYTHON_INCLUDES}")
  endif()
else()
  if(FIND_PYTHON_REQUIRED)
    message( FATAL_ERROR "Python include directory missing")
  endif()
endif()

if(PYTHON_LIBRARY)
  if(NOT PYTHON_FIND_QUIETLY)
    message( STATUS "Found Python library: ${PYTHON_LIBRARY}")
  endif()
else()
  if(FIND_PYTHON_REQUIRED)
    message( FATAL_ERROR "Python library missing")
  endif()
endif()

MARK_AS_ADVANCED(PYTHON_EXECUTABLE PYTHON_INCLUDES PYTHON_LIBRARY)
