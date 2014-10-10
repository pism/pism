# - Find FEvoR
# Find the native FEvoR includes and library
#
#  FEvoR_INCLUDES - where to find fevor_distribution.hh, etx
#  FEvoR_LIB      - List of libraries when using FEvoR.
#  FEvoR_FOUND    - True if FEvoR was found.


if (FEvoR_INCLUDES)
  # Already in cache, be silent
  set (FEvoR_FIND_QUIETLY TRUE)
endif (FEvoR_INCLUDES)

find_path (FEvoR_INCLUDES "fevor_distribution.hh"
  HINTS "${FEVOR_ROOT}/include/" "$ENV{FEVOR_ROOT}/include/"
  PATH_SUFFIXES "FEvoR")

string(REGEX REPLACE "/include/FEvoR/?$" "/lib/FEvoR"
  FEvoR_LIB_HINT ${FEvoR_INCLUDES})

find_library (FEvoR_LIB
  NAMES FEvoR
  HINTS ${FEvoR_LIB_HINT})

# handle the QUIETLY and REQUIRED arguments and set FEvoR_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FEvoR DEFAULT_MSG FEvoR_LIB FEvoR_INCLUDES)

mark_as_advanced (FEvoR_LIB FEvoR_INCLUDES)
