# - Find GSL
# Find the native GSL includes and library
#
#  GSL_INCLUDES    - where to find gsl/gsl_*.h, etc.
#  GSL_LIBRARIES   - List of libraries when using GSL.
#  GSL_FOUND       - True if GSL found.


if (GSL_INCLUDES)
  # Already in cache, be silent
  set (GSL_FIND_QUIETLY TRUE)
endif (GSL_INCLUDES)

find_path (GSL_INCLUDES gsl/gsl_math.h)

find_library (GSL_LIBRARIES NAMES gsl)

# handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (GSL DEFAULT_MSG GSL_LIBRARIES GSL_INCLUDES)

mark_as_advanced (GSL_LIBRARIES GSL_INCLUDES)
