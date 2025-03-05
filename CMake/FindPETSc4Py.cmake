# - Find PETSc4Py (cribbed from the FENICs dolfin file FindNumPy.cmake)
# Find the petsc4py include directory
#
#  PETSC4PY_INCLUDES     - where to find petsc4py/petsc4py.i, etc.
#  PETSC4PY_FOUND        - True if petsc4py is found

if(PETSC4PY_INCLUDES)
  # in cache already
  set(PETSC4PY_FIND_QUIETLY TRUE)
endif(PETSC4PY_INCLUDES)

if (NOT PETSC4PY_INCLUDES)
execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import petsc4py; print(petsc4py.get_include())"
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE PETSC4PY_INCLUDES
  RESULT_VARIABLE PETSC4PY_NOT_FOUND)

set(PETSC4PY_INCLUDES ${PETSC4PY_INCLUDES} CACHE PATH "petsc4py include directory")

execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import petsc4py; print(petsc4py.__version__)"
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE PETSC4PY_VERSION
  RESULT_VARIABLE PETSC4PY_NOT_FOUND)

set(PETSC4PY_VERSION ${PETSC4PY_VERSION} CACHE PATH "petsc4py version")
endif()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSc4Py DEFAULT_MSG PETSC4PY_INCLUDES)

mark_as_advanced(PETSC4PY_INCLUDES PETSC4PY_VERSION)
