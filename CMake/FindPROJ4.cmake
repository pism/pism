# - Find proj.4
# Find the native proj.4 includes and library
#
#  PROJ4_INCLUDES    - where to find fftw3.h
#  PROJ4_LIBRARIES   - List of libraries when using proj.4.
#  PROJ4_FOUND       - True if proj.4 found.

if (PROJ4_INCLUDES)
  # Already in cache, be silent
  set (PROJ4_FIND_QUIETLY TRUE)
endif (PROJ4_INCLUDES)

find_path (PROJ4_INCLUDES proj_api.h)

find_library (PROJ4_LIBRARIES NAMES proj)

# handle the QUIETLY and REQUIRED arguments and set PROJ4_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PROJ4 DEFAULT_MSG PROJ4_LIBRARIES PROJ4_INCLUDES)

mark_as_advanced (PROJ4_LIBRARIES PROJ4_INCLUDES)
