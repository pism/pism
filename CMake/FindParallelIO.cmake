# - Find ParallelIO
# Find ParallelIO includes and library
#
#  ParallelIO_INCLUDES    - where to find pio.h
#  ParallelIO_LIBRARIES   - List of libraries to link with when using ParallelIO
#  ParallelIO_FOUND       - True if ParallelIO was found

if (ParallelIO_INCLUDES)
  # Already in cache, be silent
  set (ParallelIO_FIND_QUIETLY TRUE)
endif (ParallelIO_INCLUDES)

find_path (ParallelIO_INCLUDES pio.h
  HINTS "${ParallelIO_ROOT}/include" "$ENV{ParallelIO_ROOT}/include")

string(REGEX REPLACE "/include/?$" "/lib"
  ParallelIO_LIB_HINT ${ParallelIO_INCLUDES})

find_library (ParallelIO_LIB
  NAMES pioc
  HINTS ${ParallelIO_LIB_HINT})

if ((NOT ParallelIO_LIB) OR (NOT ParallelIO_INCLUDES))
  message(STATUS "Trying to find ParallelIO using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(ParallelIO_LIB
    NAMES pioc
    HINTS ${LD_LIBRARY_PATH})

  if (ParallelIO_LIB)
    get_filename_component(ParallelIO_LIB_DIR ${ParallelIO_LIB} PATH)
    string(REGEX REPLACE "/lib/?$" "/include"
      ParallelIO_H_HINT ${ParallelIO_LIB_DIR})

    find_path (ParallelIO_INCLUDES pio.h
      HINTS ${ParallelIO_H_HINT}
      DOC "Path to pio.h")
  endif()
endif()

if (ParallelIO_LIB AND ParallelIO_INCLUDES)
  set (ParallelIO_LIBRARIES "${ParallelIO_LIB}")
endif()

# handle the QUIETLY and REQUIRED arguments and set ParallelIO_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (ParallelIO DEFAULT_MSG ParallelIO_LIBRARIES ParallelIO_INCLUDES)

mark_as_advanced (ParallelIO_LIB ParallelIO_INCLUDES)
