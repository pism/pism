# Find latex and create doc targets:
project(PISM_Documentation "NONE")
cmake_minimum_required (VERSION 2.8.0)

find_package(LATEX)
find_package(Doxygen)

if (NOT Pism_SOURCE_DIR)
  # This assumes that the current directory ("doc") is right in the
  # PISM source code directory, i.e. pism/doc.
  get_filename_component(Pism_SOURCE_DIR "${PROJECT_SOURCE_DIR}/.." REALPATH)
  message(STATUS "Setting Pism_SOURCE_DIR to ${Pism_SOURCE_DIR}")

  add_subdirectory(${Pism_SOURCE_DIR}/site-packages ${CMAKE_BINARY_DIR}/site-packages)

  list (APPEND CMAKE_MODULE_PATH "${Pism_SOURCE_DIR}/CMake")
endif()

# User's Manual, Climate Forcing, Installation:
add_subdirectory(sphinx)

if (PDFLATEX_COMPILER)
  # Source code browser:
  add_subdirectory (browser EXCLUDE_FROM_ALL)
endif ()                        # if PDFLATEX_COMPILER
