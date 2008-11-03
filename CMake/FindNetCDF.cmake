# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NetCDF_INCLUDE_PATH - where to find gsl/gsl_*.h, etc.
#  NetCDF_LIBRARIES   - List of libraries when using NetCDF.
#  NetCDF_FOUND       - True if NetCDF found.

if (NetCDF_INCLUDE_PATH)
  # Already in cache, be silent
  set (NetCDF_FIND_QUIETLY TRUE)
endif (NetCDF_INCLUDE_PATH)

find_path (NetCDF_INCLUDE_PATH netcdf.h)

find_library (NetCDF_LIBRARIES NAMES netcdf)

# handle the QUIETLY and REQUIRED arguments and set NetCDF_FOUND to TRUE if 
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF DEFAULT_MSG NetCDF_LIBRARIES NetCDF_INCLUDE_PATH)

mark_as_advanced (NetCDF_LIBRARIES NetCDF_INCLUDE_PATH)
