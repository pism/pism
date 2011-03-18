# - Find PETSc4Py (cribbed from the FENICs dolfin file FindNumPy.cmake)
# Find the petsc4py include directory
#
#  PETSC4PY_INCLUDES     - where to find petsc4py/petsc4py.i, etc.
#  PETSC4PY_FOUND        - True if petsc4py is found

if(PETSC4PY_INCLUDES)
  # in cache already
  set(PETSC4PY_FIND_QUIETLY TRUE)
endif(PETSC4PY_INCLUDES)

execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c "import petsc4py; from sys import stdout; stdout.write(petsc4py.get_include())"
  OUTPUT_VARIABLE PETSC4PY_INCLUDES
  RESULT_VARIABLE PETSC4PY_NOT_FOUND)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSc4Py DEFAULT_MSG PETSC4PY_INCLUDES)


# if(PETSC4PY_INCLUDES)
#   set(PETSC4PY_FOUND TRUE)
#   set(PETSC4PY_INCLUDES ${PETSC4PY_INCLUDES} CACHE STRING "PETSc4Py include path")
# else(PETSC4PY_INCLUDES)
#   set(PETSC4PY_FOUND FALSE)
# endif(PETSC4PY_INCLUDES)
# 
# if(PETSC4PY_FOUND)
#   if(NOT PETSC4PY_FIND_QUIETLY)
#     message(STATUS "Found petsc4py: ${PETSC4PY_INCLUDES}")
#   endif(NOT PETSC4PY_FIND_QUIETLY)
# else(PETSC4PY_FOUND)
#   if(PETSC4PY_FIND_REQUIRED)
#     message(FATAL_ERROR "petsc4py headers missing")
#   endif(PETSC4PY_FIND_REQUIRED)
# endif(PETSC4PY_FOUND)
# 
mark_as_advanced(PETSC4PY_INCLUDES)
