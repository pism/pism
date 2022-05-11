# - Find PROJ
# Find the native PROJ includes and library
#
#  PROJ_INCLUDES    - where to find proj.h
#  PROJ_LIBRARIES   - List of libraries when using PROJ.
#  PROJ_FOUND       - True if PROJ found.

if (PROJ_INCLUDES)
  # Already in cache, be silent
  set (PROJ_FIND_QUIETLY TRUE)
endif (PROJ_INCLUDES)

find_path (PROJ_INCLUDES proj.h
  HINTS
  "${PROJ_ROOT}/include"
  "$ENV{PROJ_ROOT}/include"
  "/opt/local/lib/proj6/include"
  "/opt/local/lib/proj7/include"
  "/opt/local/lib/proj8/include"
  "/opt/local/lib/proj9/include")

if (PROJ_INCLUDES)
  string(REGEX REPLACE "/include/?$" "/lib"
    PROJ_LIB_HINT ${PROJ_INCLUDES})
endif()

find_library (PROJ_LIBRARIES
  NAMES proj
  HINTS ${PROJ_LIB_HINT})

if ((NOT PROJ_LIBRARIES) OR (NOT PROJ_INCLUDES))
  message(STATUS "Trying to find PROJ using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(PROJ_LIBRARIES
    NAMES proj
    HINTS ${LD_LIBRARY_PATH})

  if (PROJ_LIBRARIES)
    get_filename_component(PROJ_LIB_DIR ${PROJ_LIBRARIES} PATH)
    string(REGEX REPLACE "/lib/?$" "/include"
      PROJ_H_HINT ${PROJ_LIB_DIR})

    find_path (PROJ_INCLUDES proj.h
      HINTS ${PROJ_H_HINT}
      DOC "Path to proj.h")
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set PROJ_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PROJ DEFAULT_MSG PROJ_LIBRARIES PROJ_INCLUDES)

mark_as_advanced (PROJ_LIBRARIES PROJ_INCLUDES)
