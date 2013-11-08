# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h
  HINTS "${FFTW_ROOT}/include" "$ENV{FFTW_ROOT}/include")

string(REGEX REPLACE "/include/?$" "/lib"
  FFTW_LIB_HINT ${FFTW_INCLUDES})

find_library (FFTW_LIBRARIES
  NAMES fftw3
  HINTS ${FFTW_LIB_HINT})

if ((NOT FFTW_LIBRARIES) OR (NOT FFTW_INCLUDES))
  message(STATUS "Trying to find FFTW3 using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(FFTW_LIBRARIES
    NAMES fftw3
    HINTS ${LD_LIBRARY_PATH})

  if (FFTW_LIBRARIES)
    get_filename_component(FFTW_LIB_DIR ${FFTW_LIBRARIES} PATH)
    string(REGEX REPLACE "/lib/?$" "/include"
      FFTW_H_HINT ${FFTW_LIB_DIR})

    find_path (FFTW_INCLUDES fftw3.h
      HINTS ${FFTW_H_HINT}
      DOC "Path to fftw3.h")
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
