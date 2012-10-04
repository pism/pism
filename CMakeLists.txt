project (Pism)
cmake_minimum_required (VERSION 2.8.0) # 2.8.0 is almost 3 years old, so this should be OK
include ("CMake/PISM_CMake_macros.cmake")
list (APPEND CMAKE_MODULE_PATH "${Pism_SOURCE_DIR}/CMake")

set (Pism_BRANCH "stable")

# Set Pism_REVISION_TAG
pism_set_revision_tag()
message(STATUS "Configuring PISM version '${Pism_REVISION_TAG}'")

# Put executables in the build directory:
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# Allow setting a custom install prefix using the PISM_PREFIX environment variable.
string (LENGTH "$ENV{PISM_INSTALL_PREFIX}" INSTALL_PREFIX_LENGTH)
if (INSTALL_PREFIX_LENGTH)
  set (CMAKE_INSTALL_PREFIX $ENV{PISM_INSTALL_PREFIX} CACHE PATH "PISM install prefix" FORCE)
  message (STATUS "Setting PISM install prefix to ${CMAKE_INSTALL_PREFIX}.")
endif()

# Define the directory structure.
set (Pism_BIN_DIR "bin")
set (Pism_LIB_DIR "lib/pism")
set (Pism_DOC_DIR "share/doc/pism")

# Set Pism_DEFAULT_CONFIG_FILE (*after* we set the CMAKE_INSTALL_PREFIX).
pism_check_build_dir_location()
set (Pism_DEFAULT_CONFIG_FILE "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}/pism_config.nc" CACHE STRING "" FORCE)
mark_as_advanced (Pism_DEFAULT_CONFIG_FILE)
file (WRITE ${PROJECT_BINARY_DIR}/.petscrc "-config ${PROJECT_BINARY_DIR}/pism_config.nc")

if (NOT Pism_MANUAL_BUILD)

  find_package (PETSc)
  if (NOT PETSC_FOUND)
    get_filename_component(pcc ${PETSC_COMPILER} REALPATH)
    get_filename_component(cc ${CMAKE_C_COMPILER} REALPATH)
    if (NOT ${pcc} STREQUAL ${cc})
      message(WARNING
        "PETSC_COMPILER does not match CMAKE_C_COMPILER\n"
	"  PETSC_COMPILER=${PETSC_COMPILER}\n"
	"  CMAKE_C_COMPILER=${CMAKE_C_COMPILER}\n"
	"Try running \n"
	"  rm CMakeCache.txt && cmake -DCMAKE_C_COMPILER=${PETSC_COMPILER} ${CMAKE_SOURCE_DIR}")
    endif()
    message(FATAL_ERROR  "PISM configuration failed: PETSc was not found.")
  endif()

  find_package (GSL REQUIRED)
  find_package (NetCDF REQUIRED)
  find_package (PNetCDF)

  # Try to find netcdf_par.h. We assume that NetCDF was compiled with
  # parallel I/O if this header is present.
  find_file(NETCDF_PAR_H netcdf_par.h HINTS ${NETCDF_INCLUDES} NO_DEFAULT_PATH)
  mark_as_advanced(NETCDF_PAR_H)
  if (NETCDF_PAR_H)
    message(STATUS "Note: Selected NetCDF-4 library supports parallel I/O")
    option (Pism_USE_PARALLEL_NETCDF4 "Enables parallel NetCDF-4 I/O." ON)
  else()
    message(STATUS "Note: Selected NetCDF library does not support parallel I/O")
  endif()

  if (PNETCDF_FOUND)
    message(STATUS "Note: PnetCDF was found. Enabling PnetCDF-based parallel I/O.")
    option (Pism_USE_PNETCDF "Enables parallel NetCDF-3 I/O using PnetCDF." ON)
  else ()
  endif()

  # use the PETSc compiler as a hint when looking for an MPI compiler

  # FindMPI.cmake changed between 2.8.4 and 2.8.5, so we try to support both...
  if (${CMAKE_VERSION} VERSION_LESS "2.8.5")
    set (MPI_COMPILER ${PETSC_COMPILER} CACHE FILEPATH "MPI compiler. Used only to detect MPI compilation flags.")
    find_package (MPI REQUIRED)

    set (MPI_C_INCLUDE_PATH "${MPI_INCLUDE_PATH}" CACHE STRING "MPI include directories (semicolon-separated list)")
    set (MPI_C_LIBRARIES "${MPI_LIBRARY};${MPI_EXTRA_LIBRARY}" CACHE STRING "MPI libraries (semicolon-separated list)")
    mark_as_advanced(MPI_C_INCLUDE_PATH MPI_C_LIBRARIES)
    message (STATUS
      "Note: Please upgrade CMake to version 2.8.5 or later if the build fails with undefined symbols related to MPI.")
  else ()
    set (MPI_C_COMPILER ${PETSC_COMPILER} CACHE FILEPATH "MPI compiler. Used only to detect MPI compilation flags.")
    find_package (MPI REQUIRED)
  endif()

  # FindMPI.cmake uses GetPrerequisites.cmake, which sets file_cmd
  mark_as_advanced(file_cmd MPI_LIBRARY MPI_EXTRA_LIBRARY)
  # optional
  find_package (FFTW)
  find_package (PROJ4)
else ()
  # includes
  set (PETSC_INCLUDES "EDIT_THIS" CACHE STRING "PETSc include directories (semicolon-separated list)")
  set (GSL_INCLUDES   "EDIT_THIS" CACHE STRING "GSL include directories (semicolon-separated list)")
  set (NETCDF_INCLUDES "EDIT_THIS" CACHE STRING "NetCDF include directories (semicolon-separated list)")
  set (MPI_C_INCLUDE_PATH "EDIT_THIS" CACHE STRING "MPI include directories (semicolon-separated list)")
  # optional:
  set (FFTW_INCLUDES "EDIT_THIS" CACHE STRING "FFTW include directories (semicolon-separated list)")
  set (FFTW_FOUND "OFF" CACHE BOOL "Set to 'ON' to build with FFTW3.")
  set (PROJ4_INCLUDES "EDIT_THIS" CACHE STRING "proj.4 include directories (semicolon-separated list)")
  set (PROJ4_FOUND "OFF" CACHE BOOL "Set to 'ON' to build with proj.4.")
  # libraries
  set (PETSC_LIBRARIES "EDIT_THIS" CACHE STRING "PETSc libraries (semicolon-separated list)")
  set (GSL_LIBRARIES "EDIT_THIS" CACHE STRING "GSL libraries (semicolon-separated list)")
  set (NETCDF_LIBRARIES "EDIT_THIS" CACHE STRING "NetCDF libraries (semicolon-separated list)")
  set (MPI_C_LIBRARIES "EDIT_THIS" CACHE STRING "MPI libraries (semicolon-separated list)")
  # optional
  set (FFTW_LIBRARIES "EDIT_THIS" CACHE STRING "FFTW libraries (semicolon-separated list)")
  set (PROJ4_LIBRARIES "EDIT_THIS" CACHE STRING "proj.4 libraries (semicolon-separated list)")
  # programs
  set (MPIEXEC "EDIT_THIS" CACHE FILEPATH "MPI program to run parallel tasks with")
endif()


# The default options cache
option (Pism_BUILD_EXTRA_EXECS "Build extra executables (mostly testing/verification)" OFF)
option (BUILD_SHARED_LIBS      "Build shared Pism libraries"                           ON)

option (Pism_ADD_FPIC "Add -fPIC to C++ compiler flags (CMAKE_CXX_FLAGS). Try turning it off if it does not work." ON)
option (Pism_LINK_STATICALLY "Set CMake flags to try to ensure that everything is linked statically")
option (Pism_BUILD_DEBIAN_PACKAGE "Use settings appropriate for building a .deb package" OFF)

# Use rpath by default; this has to go first, because rpath settings may be overridden later.
pism_use_rpath()

if (Pism_BUILD_DEBIAN_PACKAGE)
  pism_build_debian_package()
endif ()

if (Pism_LINK_STATICALLY)
  pism_strictly_static()
endif ()

# Deal with build types
set (Pism_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
# Force the default cache entry to obey Pism_BUILD_TYPE
set (CMAKE_BUILD_TYPE ${Pism_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None (CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
mark_as_advanced (CMAKE_BUILD_TYPE) # hide the default entry
if (Pism_BUILD_TYPE MATCHES "Debug")
  set (Pism_BUILD_EXTRA_EXECS ON CACHE BOOL "Build extra executables (mostly testing/verification)" FORCE)
  option (Pism_DEBUG              "Enables extra checks in the code."          ON)
  option (Pism_PEDANTIC_WARNINGS  "Compile with pedantic warnings."            ON)
  option (Pism_GPROF_FLAGS        "Add flags necessary to profile with gprof." OFF)
endif (Pism_BUILD_TYPE MATCHES "Debug")

# Enable extra checks in the code:
if (Pism_DEBUG)
  message (STATUS "Adding -DPISM_DEBUG=1 to compiler flags.")
  add_definitions (-DPISM_DEBUG=1)
else()
  message (STATUS "Adding -DPISM_DEBUG=0 to compiler flags.")
  add_definitions (-DPISM_DEBUG=0)
endif ()

# Add -fPIC to CXX flags; note that this will not show up in CMakeCache.txt
if (Pism_ADD_FPIC)
  if (NOT CMAKE_CXX_FLAGS MATCHES "-fPIC")
    message (STATUS "Adding -fPIC to C++ compiler flags.")
    set (CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}" CACHE STRING "C++ compiler flags" FORCE)
  endif()
endif ()

if (Pism_PEDANTIC_WARNINGS)
  pism_set_pedantic_flags()
endif (Pism_PEDANTIC_WARNINGS)

if (Pism_GPROF_FLAGS)
  set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
endif ()

include_directories (${PETSC_INCLUDES} ${GSL_INCLUDES} ${NETCDF_INCLUDES} ${MPI_C_INCLUDE_PATH})

set (Pism_EXTERNAL_LIBS
  ${PETSC_LIBRARIES} ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${MPI_LIBRARY} ${MPI_C_LIBRARIES})

if (Pism_USE_PNETCDF)
  include_directories (${PNETCDF_INCLUDES})
  list (APPEND Pism_EXTERNAL_LIBS ${PNETCDF_LIBRARIES})
endif()

# Build the Lingle-Clark bed deformation model if FFTW was found
if (FFTW_FOUND)
  add_definitions (-DPISM_HAVE_FFTW=1)
  include_directories (${FFTW_INCLUDE_DIRS} ${FFTW_INCLUDES})
  list (APPEND Pism_EXTERNAL_LIBS ${FFTW_LIBRARIES})
else (FFTW_FOUND)
  add_definitions (-DPISM_HAVE_FFTW=0)
endif (FFTW_FOUND)

# Do cell area computations the right way if proj.4 was found.
if (PROJ4_FOUND)
  add_definitions (-DPISM_HAVE_PROJ4=1)
  include_directories (${PROJ4_INCLUDES})
  list (APPEND Pism_EXTERNAL_LIBS ${PROJ4_LIBRARIES})
else()
  add_definitions (-DPISM_HAVE_PROJ4=0)
endif()

add_custom_target (etags
  COMMAND find -E src -regex ".*\\.(c|cc|h|hh)" | xargs etags
  WORKING_DIRECTORY ${Pism_SOURCE_DIR}
  VERBATIM
)

add_subdirectory (src)

# Install helper scripts residing in util/
install (PROGRAMS
  util/check_stationarity.py
  util/fill_missing.py
  util/flowline.py
  util/nc2cdo.py
  util/nc2mat.py
  util/nccmp.py
  util/pism_config_editor.py
  DESTINATION ${Pism_BIN_DIR})

if (Pism_BUILD_DOCS)
  add_subdirectory (doc)
else ()
  add_subdirectory (doc EXCLUDE_FROM_ALL)
endif()

if (Pism_BUILD_BROWSER)
  add_subdirectory (doc/browser)
else ()
  add_subdirectory (doc/browser EXCLUDE_FROM_ALL)
endif()

# PISM regression testing
ENABLE_TESTING()
add_subdirectory (test/regression)
