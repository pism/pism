project (Pism)
set (Pism_BRANCH "development")

cmake_minimum_required (VERSION 2.6.2)

# Put executables in the build directory:
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# Allow setting a custom install prefix using the PISM_PREFIX environment variable.
string (LENGTH "$ENV{PISM_INSTALL_PREFIX}" INSTALL_PREFIX_LENGTH)
if (INSTALL_PREFIX_LENGTH)
  set (CMAKE_INSTALL_PREFIX $ENV{PISM_INSTALL_PREFIX} CACHE PATH "PISM install prefix" FORCE)
  message (STATUS "Setting PISM install prefix to ${CMAKE_INSTALL_PREFIX}.")
endif()

list (APPEND CMAKE_MODULE_PATH "${Pism_SOURCE_DIR}/CMake")


if (NOT Pism_MANUAL_BUILD)

  find_package (PETSc)
  if (NOT PETSC_FOUND)
    get_filename_component(pcc ${PETSC_COMPILER} REALPATH)
    get_filename_component(cc ${CMAKE_C_COMPILER} REALPATH)
    if (NOT ${pcc} STREQUAL ${cc})
      message(WARNING "PETSC_COMPILER does not match CMAKE_C_COMPILER\n"
	"  PETSC_COMPILER=${PETSC_COMPILER}\n"
	"  CMAKE_C_COMPILER=${CMAKE_C_COMPILER}\n"
	"Try running \n"
	"  rm CMakeCache.txt && cmake -DCMAKE_C_COMPILER=${PETSC_COMPILER} ${CMAKE_SOURCE_DIR}")
    endif()
    message(FATAL_ERROR  "PISM configuration failed: PETSc was not found.")
  endif()

  find_package (GSL REQUIRED)
  find_package (NetCDF REQUIRED)

  # Try to find netcdf_par.h. We assume that NetCDF was compiled with
  # parallel I/O if this header is present.
  find_file(NETCDF_PAR_H netcdf_par.h HINTS ${NETCDF_INCLUDES} NO_DEFAULT_PATH)
  mark_as_advanced(NETCDF_PAR_H)
  if (NETCDF_PAR_H)
    message(STATUS "Note: Selected NetCDF library supports parallel I/O")
    option (Pism_USE_PARALLEL_NETCDF "Enables parallel NetCDF I/O." ON)
  else()
    message(STATUS "Note: Selected NetCDF library does not support parallel I/O")
  endif()

  # use the PETSc compiler as a hint when looking for an MPI compiler
  set (MPI_COMPILER ${PETSC_COMPILER} CACHE FILEPATH "MPI compiler. Used only to detect MPI compilation flags.")
  find_package (MPI REQUIRED)
  # optional
  find_package (FFTW)
  find_package (PROJ4)
else ()
  # includes
  set (PETSC_INCLUDES "EDIT_THIS" CACHE STRING "PETSc include directories (semicolon-separated list)")
  set (GSL_INCLUDES   "EDIT_THIS" CACHE STRING "GSL include directories (semicolon-separated list)")
  set (NETCDF_INCLUDES "EDIT_THIS" CACHE STRING "NetCDF include directories (semicolon-separated list)")
  set (MPI_INCLUDE_PATH "EDIT_THIS" CACHE STRING "MPI include directories (semicolon-separated list)")
  # optional:
  set (FFTW_INCLUDES "EDIT_THIS" CACHE STRING "FFTW include directories (semicolon-separated list)")
  set (FFTW_FOUND "OFF" CACHE BOOL "Set to 'ON' to build with FFTW3.")
  set (PROJ4_INCLUDES "EDIT_THIS" CACHE STRING "proj.4 include directories (semicolon-separated list)")
  set (PROJ4_FOUND "OFF" CACHE BOOL "Set to 'ON' to build with proj.4.")
  # libraries
  set (PETSC_LIBRARIES "EDIT_THIS" CACHE STRING "PETSc libraries (semicolon-separated list)")
  set (GSL_LIBRARIES "EDIT_THIS" CACHE STRING "GSL libraries (semicolon-separated list)")
  set (NETCDF_LIBRARIES "EDIT_THIS" CACHE STRING "NetCDF libraries (semicolon-separated list)")
  set (MPI_LIBRARIES "EDIT_THIS" CACHE STRING "MPI libraries (semicolon-separated list)")
  # optional
  set (FFTW_LIBRARIES "EDIT_THIS" CACHE STRING "FFTW libraries (semicolon-separated list)")
  set (PROJ4_LIBRARIES "EDIT_THIS" CACHE STRING "proj.4 libraries (semicolon-separated list)")
  # programs
  set (MPIEXEC "EDIT_THIS" CACHE FILEPATH "MPI program to run parallel tasks with")
endif()


# The default options cache
option (Pism_BUILD_EXTRA_EXECS "Build extra executables (mostly testing/verification)" OFF)
option (BUILD_SHARED_LIBS      "Build shared Pism libraries"                           ON)
option (Pism_BUILD_PYTHON_BINDINGS "Build (experimental) python bindings" OFF)

option (Pism_ADD_FPIC "Add -fPIC to C++ compiler flags (CMAKE_CXX_FLAGS). Try turning it off if it does not work." ON)
option (Pism_RESET_LINKER_FLAGS
  "Set linker flags using CMake variables PISM_SHARED_LIBRARY_LINK_C_FLAGS and PISM_SHARED_LIBRARY_LINK_CXX_FLAGS")

set (PISM_SHARED_LIBRARY_LINK_C_FLAGS "" CACHE STRING
  "Sets CMAKE_SHARED_LIBRARY_LINK_C_FLAGS to this if Pism_RESET_LINKER_FLAGS is ON")
set (PISM_SHARED_LIBRARY_LINK_CXX_FLAGS "" CACHE STRING
  "Sets CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS to this if Pism_RESET_LINKER_FLAGS is ON")
mark_as_advanced (PISM_SHARED_LIBRARY_LINK_C_FLAGS PISM_SHARED_LIBRARY_LINK_CXX_FLAGS
  Pism_ADD_FPIC Pism_RESET_LINKER_FLAGS)


option (Pism_BUILD_DEBIAN_PACKAGE "Use settings appropriate for building a .deb package" OFF)

set (Pism_BIN_DIR "bin")
set (Pism_LIB_DIR "lib/pism")
set (Pism_DOC_DIR "share/doc/pism")

if (Pism_BUILD_DEBIAN_PACKAGE)
  set (Pism_BUILD_TYPE "Release" CACHE STRING "PISM build type" FORCE)
  set (CMAKE_INSTALL_PREFIX "/usr" CACHE STRING "Install prefix" FORCE)
  set (Pism_BUILD_DOCS OFF CACHE BOOL "Build PISM documentation" FORCE)
  set (Pism_BUILD_BROWSER OFF CACHE BOOL "Build PISM source code browsers" FORCE)
  set (Pism_BUILD_EXTRA_EXECS OFF CACHE BOOL "Build extra executables (mostly testing/verification)" FORCE)

  # RPATH handling
  set (CMAKE_SKIP_BUILD_RPATH TRUE)
  set (CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
  set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}")
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
else ()

  ## Use full RPATH, with this setting Pism libraries cannot be moved after installation
  ## but the correct libraries will always be found regardless of LD_LIBRARY_PATH
  ## in use, i.e. don't skip the full RPATH for the build tree
  set (CMAKE_SKIP_BUILD_RPATH FALSE)
  # when building, don't use the install RPATH already
  # (but later on when installing)
  set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
  # the RPATH to be used when installing
  set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}")
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Mac OS X install_name fix:
  set (CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}")
endif ()

# Deal with build types
set (Pism_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
# Force the default cache entry to obey Pism_BUILD_TYPE
set (CMAKE_BUILD_TYPE ${Pism_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None (CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
mark_as_advanced (CMAKE_BUILD_TYPE) # hide the default entry
if (Pism_BUILD_TYPE MATCHES "Debug")
  option (Pism_DEBUG              "Enables extra checks in the code."          OFF)
  option (Pism_PEDANTIC_WARNINGS  "Compile with pedantic warnings."            OFF)
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

# This resets linker options.
if (Pism_RESET_LINKER_FLAGS)
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS ${PISM_SHARED_LIBRARY_LINK_C_FLAGS})
  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS ${PISM_SHARED_LIBRARY_LINK_CXX_FLAGS})
endif()

# Set Pism_REVISION
# Git
if (EXISTS ${Pism_SOURCE_DIR}/.git)
  find_program (GIT git DOC "Git executable")
  mark_as_advanced(GIT)
  execute_process (COMMAND ${GIT} describe --always --match v?.?
    WORKING_DIRECTORY ${Pism_SOURCE_DIR}
    OUTPUT_VARIABLE Pism_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif (EXISTS ${Pism_SOURCE_DIR}/.git)

# Otherwise...
if (NOT Pism_VERSION)
  set (Pism_VERSION "unknown")
endif (NOT Pism_VERSION)

set (Pism_REVISION_TAG "${Pism_BRANCH} ${Pism_VERSION}")

# set Pism_DEFAULT_CONFIG_FILE
if (NOT Pism_DEFAULT_CONFIG_FILE)
  set (Pism_DEFAULT_CONFIG_FILE "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}/pism_config.nc")
endif ()
file (WRITE ${PROJECT_BINARY_DIR}/.petscrc "-config ${PROJECT_BINARY_DIR}/pism_config.nc")

# We could set this here, but it forces a full recompile if the version
# changes, hence we'll define PISM_REVISION only for src/base/iMutil.cc
# (see src/CMakeLists.txt)
#add_definitions (-DPISM_REVISION=\\\"${Pism_REVISION}\\\")

if (Pism_PEDANTIC_WARNINGS)
  set (DEFAULT_PEDANTIC_CFLAGS "-pedantic -Wall -Wextra -Wno-cast-qual -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wconversion -Wsign-compare -Wno-redundant-decls -Winline -Wno-long-long -Wmissing-format-attribute -Wmissing-noreturn -Wpacked -Wdisabled-optimization -Wmultichar -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wendif-labels -Winvalid-pch -Wmissing-field-initializers -Wvariadic-macros -Wstrict-aliasing -funit-at-a-time")
  set (DEFAULT_PEDANTIC_CXXFLAGS "${DEFAULT_PEDANTIC_CFLAGS} -Woverloaded-virtual")
  set (PEDANTIC_CFLAGS ${DEFAULT_PEDANTIC_CFLAGS} CACHE STRING "Compiler flags to enable pedantic warnings")
  set (PEDANTIC_CXXFLAGS ${DEFAULT_PEDANTIC_CXXFLAGS} CACHE STRING "Compiler flags to enable pedantic warnings for C++")
  mark_as_advanced (PEDANTIC_CFLAGS PEDANTIC_CXXFLAGS)
  set (CMAKE_C_FLAGS_DEBUG "-g ${PEDANTIC_CFLAGS}")
  set (CMAKE_CXX_FLAGS_DEBUG "-g ${PEDANTIC_CXXFLAGS}")
endif (Pism_PEDANTIC_WARNINGS)

if (Pism_GPROF_FLAGS)
  set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
endif ()


if (Pism_BUILD_PYTHON_BINDINGS)
  find_package(Python REQUIRED)
  find_package(PETSc4Py REQUIRED)
  find_package(SWIG REQUIRED)

  mark_as_advanced (SWIG_DIR SWIG_EXECUTABLE SWIG_VERSION)

  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import distutils.sysconfig; print distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}')"
    OUTPUT_VARIABLE PISM_INSTALL_PYTHON_MODULE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(PISM_INSTALL_PYTHON_MODULE_DIR ${PISM_INSTALL_PYTHON_MODULE_DIR}
    CACHE PATH "Python extension module installation directory." FORCE)
endif ()

include_directories (${PETSC_INCLUDES} ${GSL_INCLUDES} ${NETCDF_INCLUDES} ${MPI_INCLUDE_PATH} ${FFTW_INCLUDES})
set (Pism_EXTERNAL_LIBS ${PETSC_LIBRARIES} ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${MPI_LIBRARIES})

# Build the Lingle-Clark bed deformation model if FFTW was found
if (FFTW_FOUND)
  add_definitions (-DPISM_HAVE_FFTW=1)
  include_directories (${FFTW_INCLUDE_DIRS})
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

# Install some python-based pieces of PISM
install (PROGRAMS
  examples/python/ssa_forward.py
  examples/inverse/make_synth_ssa.py
  examples/inverse/vel2tauc.py
  DESTINATION ${Pism_BIN_DIR} )

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
