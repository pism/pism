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

find_path (GSL_INCLUDES gsl/gsl_math.h
  HINTS "${GSL_ROOT}/include" "$ENV{GSL_ROOT}/include")

string(REGEX REPLACE "/include/?$" "/lib"
  GSL_LIB_HINT ${GSL_INCLUDES})

find_library (GSL_LIB
  NAMES gsl
  HINTS ${GSL_LIB_HINT})

find_library(GSL_CBLAS_LIB
  NAMES gslcblas
  HINTS ${GSL_LIB_HINT})

if ((NOT GSL_INCLUDES) OR (NOT GSL_LIB) OR (NOT GSL_CBLAS_LIB))
  message(STATUS "Trying to find GSL using 'gsl-config'...")
  find_program(GSL_CONFIG "gsl-config")
  if (GSL_CONFIG)
    execute_process(COMMAND ${GSL_CONFIG} --prefix
      OUTPUT_VARIABLE GSL_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    find_path(GSL_INCLUDES gsl/gsl_math.h
      HINTS "${GSL_PREFIX}/include")

    find_library (GSL_LIB NAMES gsl
      HINTS "${GSL_PREFIX}/lib")

    find_library(GSL_CBLAS_LIB NAMES gslcblas
      HINTS "${GSL_PREFIX}/lib")
  endif()
endif()

if ((NOT GSL_LIB) OR (NOT GSL_INCLUDES))
  message(STATUS "Trying to find GSL using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(GSL_LIB
    NAMES gsl
    HINTS ${LD_LIBRARY_PATH})

  find_library(GSL_CBLAS_LIB
    NAMES gslcblas
    HINTS ${LD_LIBRARY_PATH})

  if (GSL_LIB)
    get_filename_component(GSL_LIB_DIR ${GSL_LIB} PATH)
    string(REGEX REPLACE "/lib/?$" "/include"
      GSL_H_HINT ${GSL_LIB_DIR})

    find_path (GSL_INCLUDES gsl/gsl_math.h
      HINTS ${GSL_H_HINT}
      DOC "Path to gsl/gsl_math.h")
  endif()
endif()

if (GSL_LIB AND GSL_CBLAS_LIB)
  set (GSL_LIBRARIES "${GSL_LIB}" "${GSL_CBLAS_LIB}")
endif()

# handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (GSL DEFAULT_MSG GSL_LIBRARIES GSL_INCLUDES)

mark_as_advanced (GSL_LIB GSL_CBLAS_LIB GSL_INCLUDES)
