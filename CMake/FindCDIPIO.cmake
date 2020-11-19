# - Find CDIPIO
# Find CDIPIO includes and library
#
#  CDIPIO_INCLUDES    - where to find pio.h
#  CDIPIO_LIBRARIES   - List of libraries to link with when using CDIPIO
#  CDIPIO_FOUND       - True if CDIPIO was found

if (CDIPIO_INCLUDES)
  # Already in cache, be silent
  set (CDIPIO_FIND_QUIETLY TRUE)
endif (CDIPIO_INCLUDES)

find_path (CDIPIO_INCLUDES cdipio.h
  HINTS "${CDIPIO_ROOT}/include" "$ENV{CDIPIO_ROOT}/include")

string(REGEX REPLACE "/include/?$" "/lib"
  CDIPIO_LIB_HINT ${CDIPIO_INCLUDES})

find_library (CDIPIO_LIB
  NAMES cdipio
  HINTS ${CDIPIO_LIB_HINT})

if ((NOT CDIPIO_LIB) OR (NOT CDIPIO_INCLUDES))
  message(STATUS "Trying to find CDIPIO using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(CDIPIO_LIB
    NAMES cdipio
    HINTS ${LD_LIBRARY_PATH})

  if (CDIPIO_LIB)
    get_filename_component(CDIPIO_LIB_DIR ${CDIPIO_LIB} PATH)
    string(REGEX REPLACE "/lib/?$" "/include"
      CDIPIO_H_HINT ${CDIPIO_LIB_DIR})

    find_path (CDIPIO_INCLUDES cdipio.h
      HINTS ${CDIPIO_H_HINT}
      DOC "Path to cdipio.h")
  endif()
endif()

if (CDIPIO_LIB AND CDIPIO_INCLUDES)
  set (CDIPIO_LIBRARIES "${CDIPIO_LIB}")
endif()

# handle the QUIETLY and REQUIRED arguments and set CDIPIO_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (CDIPIO DEFAULT_MSG CDIPIO_LIBRARIES CDIPIO_INCLUDES)

mark_as_advanced (CDIPIO_LIB CDIPIO_INCLUDES)
