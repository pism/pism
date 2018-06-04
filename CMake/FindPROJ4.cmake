# - Find proj.4
# Find the native proj.4 includes and library
#
#  PROJ4_INCLUDES    - where to find proj_api.h
#  PROJ4_LIBRARIES   - List of libraries when using proj.4.
#  PROJ4_FOUND       - True if proj.4 found.

if (PROJ4_INCLUDES)
  # Already in cache, be silent
  set (PROJ4_FIND_QUIETLY TRUE)
endif (PROJ4_INCLUDES)

find_path (PROJ4_INCLUDES proj_api.h
  HINTS
  "${PROJ_ROOT}/include"
  "$ENV{PROJ_ROOT}/include"
  "/opt/local/lib/proj5/include"
  "/opt/local/lib/proj49/include")

if (PROJ4_INCLUDES)
  string(REGEX REPLACE "/include/?$" "/lib"
    PROJ4_LIB_HINT ${PROJ4_INCLUDES})
endif()

find_library (PROJ4_LIBRARIES
  NAMES proj
  HINTS ${PROJ4_LIB_HINT})

if ((NOT PROJ4_LIBRARIES) OR (NOT PROJ4_INCLUDES))
  message(STATUS "Trying to find proj.4 using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(PROJ4_LIBRARIES
    NAMES proj
    HINTS ${LD_LIBRARY_PATH})

  if (PROJ4_LIBRARIES)
    get_filename_component(PROJ4_LIB_DIR ${PROJ4_LIBRARIES} PATH)
    string(REGEX REPLACE "/lib/?$" "/include"
      PROJ4_H_HINT ${PROJ4_LIB_DIR})

    find_path (PROJ4_INCLUDES proj_api.h
      HINTS ${PROJ4_H_HINT}
      DOC "Path to proj_api.h")
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set PROJ4_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PROJ4 DEFAULT_MSG PROJ4_LIBRARIES PROJ4_INCLUDES)

mark_as_advanced (PROJ4_LIBRARIES PROJ4_INCLUDES)
