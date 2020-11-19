# - Find YAXT
# Find YAXT includes and library
#
#  YAXT_INCLUDES    - where to find pio.h
#  YAXT_LIBRARIES   - List of libraries to link with when using YAXT
#  YAXT_FOUND       - True if YAXT was found

if (YAXT_INCLUDES)
  # Already in cache, be silent
  set (YAXT_FIND_QUIETLY TRUE)
endif (YAXT_INCLUDES)

find_path (YAXT_INCLUDES yaxt.h
  HINTS "${YAXT_ROOT}/include" "$ENV{YAXT_ROOT}/include")

string(REGEX REPLACE "/include/?$" "/lib"
  YAXT_LIB_HINT ${YAXT_INCLUDES})

find_library (YAXT_LIB
  NAMES yaxt
  HINTS ${YAXT_LIB_HINT})

if ((NOT YAXT_LIB) OR (NOT YAXT_INCLUDES))
  message(STATUS "Trying to find YAXT using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(YAXT_LIB
    NAMES yaxt
    HINTS ${LD_LIBRARY_PATH})

  if (YAXT_LIB)
    get_filename_component(YAXT_LIB_DIR ${YAXT_LIB} PATH)
    string(REGEX REPLACE "/lib/?$" "/include"
      YAXT_H_HINT ${YAXT_LIB_DIR})

    find_path (YAXT_INCLUDES yaxt.h
      HINTS ${YAXT_H_HINT}
      DOC "Path to yaxt.h")
  endif()
endif()

if (YAXT_LIB AND YAXT_INCLUDES)
  set (YAXT_LIBRARIES "${YAXT_LIB}")
endif()

# handle the QUIETLY and REQUIRED arguments and set YAXT_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (YAXT DEFAULT_MSG YAXT_LIBRARIES YAXT_INCLUDES)

mark_as_advanced (YAXT_LIB YAXT_INCLUDES)
