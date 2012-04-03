# - Find PNetCDF
# Find the native PNetCDF includes and library
#
#  PNETCDF_INCLUDES    - where to find netcdf.h, etc
#  PNETCDF_LIBRARIES   - Link these libraries when using NetCDF
#  PNETCDF_FOUND       - True if PNetCDF was found
#
# Normal usage would be:
#  find_package (PNetCDF REQUIRED)
#  target_link_libraries (uses_pnetcdf ${PNETCDF_LIBRARIES})

if (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)
  # Already in cache, be silent
  set (PNETCDF_FIND_QUIETLY TRUE)
endif (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)

find_path (PNETCDF_INCLUDES pnetcdf.h HINTS $ENV{PNETCDF_DIR}/include)

find_library (PNETCDF_LIBRARIES NAMES pnetcdf HINTS $ENV{PNETCDF_DIR}/lib)

# handle the QUIETLY and REQUIRED arguments and set PNETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PNetCDF DEFAULT_MSG PNETCDF_LIBRARIES PNETCDF_INCLUDES)

mark_as_advanced (PNETCDF_LIBRARIES PNETCDF_INCLUDES)
