# Find latex and create doc targets:
project(PISM_Documentation "NONE")
cmake_minimum_required (VERSION 3.13)

find_package(LATEX COMPONENTS PDFLATEX)
find_package(Doxygen)

if (NOT Pism_SOURCE_DIR)
  # This assumes that the current directory ("doc") is right in the
  # PISM source code directory, i.e. pism/doc.
  get_filename_component(Pism_SOURCE_DIR "${PROJECT_SOURCE_DIR}/.." REALPATH)
  message(STATUS "Setting Pism_SOURCE_DIR to ${Pism_SOURCE_DIR}")

  list (APPEND CMAKE_MODULE_PATH "${Pism_SOURCE_DIR}/CMake")

  include ("${Pism_SOURCE_DIR}/CMake/PISM_CMake_macros.cmake")

  if(EXISTS ${Pism_SOURCE_DIR}/VERSION)
    file(STRINGS ${Pism_SOURCE_DIR}/VERSION Pism_VERSION LIMIT_COUNT 1)
  else()
    message(FATAL_ERROR "'VERSION' not found: cannot determine PISM's version")
  endif()

  pism_set_full_version()

  include (GNUInstallDirs)
endif()

# User's Manual, Climate Forcing, Installation:
add_subdirectory(sphinx)

if (DOXYGEN_EXECUTABLE AND DOXYGEN_DOT_EXECUTABLE)
  # Source code browser:
  add_subdirectory (browser EXCLUDE_FROM_ALL)
endif ()
