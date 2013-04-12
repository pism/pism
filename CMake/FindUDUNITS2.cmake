# - Find UDUNITS2
# Find the native UDUNITS2 includes and library
#
#  UDUNITS2_INCLUDES    - where to find udunits2.h
#  UDUNITS2_LIBRARY     - The UDUNITS2 library
#  UDUNITS2_FOUND       - True if UDUNITS2 was found.

if (UDUNITS2_INCLUDES)
  # Already in cache, be silent
  set (UDUNITS2_FIND_QUIETLY TRUE)
endif (UDUNITS2_INCLUDES)

find_path (UDUNITS2_INCLUDES udunits2.h
  PATH_SUFFIXES "udunits2"
  DOC "Path to udunits2.h")

find_library (UDUNITS2_LIBRARY NAMES udunits2)

# handle the QUIETLY and REQUIRED arguments and set UDUNITS2_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (UDUNITS2 DEFAULT_MSG UDUNITS2_LIBRARY UDUNITS2_INCLUDES)

mark_as_advanced (UDUNITS2_LIBRARY UDUNITS2_INCLUDES)
