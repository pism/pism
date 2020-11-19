# - Find YAXT_C
# Find YAXT_C includes and library
#
#  YAXT_C_INCLUDES    - where to find pio.h
#  YAXT_C_LIBRARIES   - List of libraries to link with when using YAXT_C
#  YAXT_C_FOUND       - True if YAXT_C was found

if (YAXT_C_INCLUDES)
  # Already in cache, be silent
  set (YAXT_C_FIND_QUIETLY TRUE)
endif (YAXT_C_INCLUDES)

find_path (YAXT_C_INCLUDES yaxt.h
  HINTS "${YAXT_C_ROOT}/include" "$ENV{YAXT_C_ROOT}/include")

string(REGEX REPLACE "/include/?$" "/lib"
  YAXT_C_LIB_HINT ${YAXT_C_INCLUDES})

find_library (YAXT_C_LIB
  NAMES yaxt_c
  HINTS ${YAXT_C_LIB_HINT})

if ((NOT YAXT_C_LIB) OR (NOT YAXT_C_INCLUDES))
  message(STATUS "Trying to find YAXT_C using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(YAXT_C_LIB
    NAMES yaxt_c
    HINTS ${LD_LIBRARY_PATH})

  if (YAXT_C_LIB)
    get_filename_component(YAXT_C_LIB_DIR ${YAXT_C_LIB} PATH)
    string(REGEX REPLACE "/lib/?$" "/include"
      YAXT_C_H_HINT ${YAXT_C_LIB_DIR})

    find_path (YAXT_C_INCLUDES yaxt.h
      HINTS ${YAXT_C_H_HINT}
      DOC "Path to yaxt.h")
  endif()
endif()

if (YAXT_C_LIB AND YAXT_C_INCLUDES)
  set (YAXT_C_LIBRARIES "${YAXT_C_LIB}")
endif()

# handle the QUIETLY and REQUIRED arguments and set YAXT_C_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (YAXT_C DEFAULT_MSG YAXT_C_LIBRARIES YAXT_C_INCLUDES)

mark_as_advanced (YAXT_C_LIB YAXT_C_INCLUDES)
