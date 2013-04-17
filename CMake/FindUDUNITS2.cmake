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

set(UDUNITS2_TEST_SRC "
#include <udunits2.h>

int main(int argc, char **argv) {
  ut_system *s = ut_read_xml(NULL);
  ut_free_system(s);
  return 0;
}
")

include (CheckCSourceRuns)

set(CMAKE_REQUIRED_INCLUDES ${UDUNITS2_INCLUDES})
set(CMAKE_REQUIRED_LIBRARIES ${UDUNITS2_LIBRARY})
check_c_source_runs("${UDUNITS2_TEST_SRC}" UDUNITS2_WORKS_WITHOUT_EXPAT)

if(${UDUNITS2_WORKS_WITHOUT_EXPAT})
  message(STATUS "UDUNITS-2 does not require expat")
else()
  find_package(EXPAT REQUIRED)

  set(CMAKE_REQUIRED_INCLUDES ${UDUNITS2_INCLUDES} ${EXPAT_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${UDUNITS2_LIBRARY} ${EXPAT_LIBRARIES})
  check_c_source_runs("${UDUNITS2_TEST_SRC}" UDUNITS2_WORKS_WITH_EXPAT)

  if(NOT ${UDUNITS2_WORKS_WITH_EXPAT})
    message(FATAL_ERROR "Failed to find UDUNITS-2")
  endif()

  message(STATUS "UDUNITS-2 requires EXPAT")
  set (UDUNITS2_LIBRARY "${UDUNITS2_LIBRARY};${EXPAT_LIBRARIES}" CACHE STRING "" FORCE)
endif()

# handle the QUIETLY and REQUIRED arguments and set UDUNITS2_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (UDUNITS2 DEFAULT_MSG UDUNITS2_LIBRARY UDUNITS2_INCLUDES)

mark_as_advanced (UDUNITS2_LIBRARY UDUNITS2_INCLUDES)
