# - Find CDI
# Find CDI includes and library
#
#  CDI_INCLUDES    - where to find pio.h
#  CDI_LIBRARIES   - List of libraries to link with when using CDI
#  CDI_FOUND       - True if CDI was found

if (CDI_INCLUDES)
  # Already in cache, be silent
  set (CDI_FIND_QUIETLY TRUE)
endif (CDI_INCLUDES)

find_path (CDI_INCLUDES cdi.h
  HINTS "${CDI_ROOT}/include" "$ENV{CDI_ROOT}/include")

string(REGEX REPLACE "/include/?$" "/lib"
  CDI_LIB_HINT ${CDI_INCLUDES})

find_library (CDI_LIB
  NAMES cdi
  HINTS ${CDI_LIB_HINT})

if ((NOT CDI_LIB) OR (NOT CDI_INCLUDES))
  message(STATUS "Trying to find CDI using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(CDI_LIB
    NAMES cdi
    HINTS ${LD_LIBRARY_PATH})

  if (CDI_LIB)
    get_filename_component(CDI_LIB_DIR ${CDI_LIB} PATH)
    string(REGEX REPLACE "/lib/?$" "/include"
      CDI_H_HINT ${CDI_LIB_DIR})

    find_path (CDI_INCLUDES cdi.h
      HINTS ${CDI_H_HINT}
      DOC "Path to cdi.h")
  endif()
endif()

if (CDI_LIB AND CDI_INCLUDES)
  set (CDI_LIBRARIES "${CDI_LIB}")
endif()

# handle the QUIETLY and REQUIRED arguments and set CDI_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (CDI DEFAULT_MSG CDI_LIBRARIES CDI_INCLUDES)

mark_as_advanced (CDI_LIB CDI_INCLUDES)
