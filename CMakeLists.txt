project (Pism)
set (Pism_BRANCH "trunk")
set (Pism_VERSION_MAJOR "0")
set (Pism_VERSION_MINOR "3")

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

find_package (PETSc REQUIRED)
find_package (GSL REQUIRED)
find_package (NetCDF REQUIRED)
find_package (FFTW)

# use the PETSc compiler as a hint when looking for an MPI compiler
find_package (MPI REQUIRED)
set (MPI_COMPILER ${PETSC_COMPILER} CACHE FILEPATH "MPI compiler. Used only to detect MPI compilation flags.")

# The default options cache
#option (Pism_BUILD_TESTS       "Build tests"                                           ON)
option (Pism_BUILD_EXTRA_EXECS "Build extra executables (mostly testing/verification)" OFF)
option (BUILD_SHARED_LIBS      "Build shared Pism libraries"                           ON)
option (Pism_LOCAL_INSTALL     "Install PISM locally (not in /usr/local/ or similar). This affects installation paths." ON)
option (Pism_ADD_FPIC "Add -fPIC to C++ compiler flags (CMAKE_CXX_FLAGS). Try turning it off if it does not work." ON)

option (Pism_BUILD_DEBIAN_PACKAGE "Use settings appropriate for building a .deb package" OFF)

if (Pism_BUILD_DEBIAN_PACKAGE)
  set (Pism_LOCAL_INSTALL OFF CACHE BOOL "Install PISM locally" FORCE)
  set (Pism_BUILD_TYPE "Release" CACHE STRING "PISM build type" FORCE)
  set (CMAKE_INSTALL_PREFIX "/usr" CACHE STRING "Install prefix" FORCE)
  set (Pism_BUILD_DOCS OFF CACHE BOOL "Build PISM documentation" FORCE)
  set (Pism_BUILD_EXTRA_EXECS OFF CACHE BOOL "Build extra executables (mostly testing/verification)" FORCE)

  # RPATH handling
  set (CMAKE_SKIP_BUILD_RPATH TRUE)
  set (CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
  set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/pism")
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


set (Pism_FAST_WRITE "1" CACHE STRING "Variable order options: 1,2,3 (fastest to slowest).")

# This installs PISM in ${CMAKE_INSTALL_PREFIX}/{bin,lib,doc} if installed
# locally; otherwise in ${CMAKE_INSTALL_PREFIX}/bin/pism, .../lib/pism and
# .../share/doc/pism to avoid cluttering system directories.
if (Pism_LOCAL_INSTALL)
  set (Pism_BIN_DIR "bin")
  set (Pism_LIB_DIR "lib")
  set (Pism_DOC_DIR "doc")
else ()
  set (Pism_BIN_DIR "bin/pism")
  set (Pism_LIB_DIR "lib/pism")
  set (Pism_DOC_DIR "share/doc/pism")
endif ()

# Deal with build types
set (Pism_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
# Force the default cache entry to obey Pism_BUILD_TYPE
set (CMAKE_BUILD_TYPE ${Pism_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
mark_as_advanced (CMAKE_BUILD_TYPE) # hide the default entry
if (Pism_BUILD_TYPE MATCHES "Debug")
  option (Pism_DEBUG              "Enables extra checks in the code."          OFF)
  option (Pism_PEDANTIC_WARNINGS  "Compile with pedantic warnings."            OFF)
  option (Pism_GPROF_FLAGS        "Add flags necessary to profile with gprof." OFF)
endif (Pism_BUILD_TYPE MATCHES "Debug")

# Enable extra checks in the code:
if (Pism_DEBUG)
  message (STATUS "Adding -DPISM_DEBUG to compiler flags.")
  add_definitions (-DPISM_DEBUG=1)
endif ()

# Add -fPIC to CXX flags; note that this will not show up in CMakeCache.txt
if (Pism_ADD_FPIC)
  message (STATUS "Adding -fPIC to C++ compiler flags.")
  set (CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}")
endif ()

# Set Pism_REVISION
# Git
if (EXISTS ${Pism_SOURCE_DIR}/.git)
  find_program (GIT git DOC "Git executable")
  execute_process (COMMAND ${GIT} svn info
    WORKING_DIRECTORY ${Pism_SOURCE_DIR}
    OUTPUT_VARIABLE Pism_SVN_INFO)
  string (REGEX REPLACE "^(.*\n)?URL: ([^\n]+).*" "\\2"
    Pism_SVN_URL "${Pism_SVN_INFO}")
  string (REGEX REPLACE "^(.*\n)?Revision: ([^\n]+).*" "\\2"
    Pism_SVN_REVISION "${Pism_SVN_INFO}")
  execute_process (COMMAND ${GIT} log --pretty=format:"%h %ci %s" -1
    WORKING_DIRECTORY ${Pism_SOURCE_DIR}
    OUTPUT_VARIABLE Pism_GIT_LOG)
  string (REGEX REPLACE "\"(.*)\"" "\\1" Pism_GIT_LOG "${Pism_GIT_LOG}")
  string (REGEX REPLACE "'" "" Pism_GIT_LOG "${Pism_GIT_LOG}")
  set (Pism_REVISION "${Pism_SVN_REVISION}")
endif (EXISTS ${Pism_SOURCE_DIR}/.git)
# Subversion
if (EXISTS ${Pism_SOURCE_DIR}/.svn)
  find_package (Subversion REQUIRED)
  Subversion_WC_INFO (${Pism_SOURCE_DIR} Pism)
  set (Pism_REVISION "${Pism_WC_REVISION}")
endif (EXISTS ${Pism_SOURCE_DIR}/.svn)
# Otherwise...
if (NOT Pism_REVISION)
  set (Pism_REVISION "0")
endif (NOT Pism_REVISION)

set (Pism_REVISION_TAG "${Pism_BRANCH} ${Pism_VERSION_MAJOR}.${Pism_VERSION_MINOR}.${Pism_REVISION}")

# set Pism_DEFAULT_CONFIG_FILE
if (NOT Pism_DEFAULT_CONFIG_FILE)
  set (Pism_DEFAULT_CONFIG_FILE "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}/pism_config.nc")
endif ()

# We could set this here, but it forces a full recompile if the version
# changes, hence we'll define PISM_REVISION only for src/base/iMutil.cc
# (see src/CMakeLists.txt)
#add_definitions (-DPISM_REVISION=\\\"${Pism_REVISION}\\\")

if (Pism_PEDANTIC_WARNINGS)
  set (DEFAULT_PEDANTIC_CFLAGS "-pedantic -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wsign-compare -Wno-redundant-decls -Winline -Wno-long-long -Wmissing-format-attribute -Wmissing-noreturn -Wpacked -Wdisabled-optimization -Wmultichar -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wendif-labels -Winvalid-pch -Wmissing-field-initializers -Wvariadic-macros -Wunsafe-loop-optimizations -Wvolatile-register-var -Wstrict-aliasing -funit-at-a-time")
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

include_directories (${PETSC_INCLUDES} ${GSL_INCLUDES} ${NETCDF_INCLUDES} ${MPI_INCLUDE_PATH})
set (Pism_EXTERNAL_LIBS ${PETSC_LIBRARIES} ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${MPI_LIBRARIES})

if (FFTW_FOUND)
  add_definitions (-DPISM_HAVE_FFTW=1)
  include_directories (${FFTW_INCLUDE_DIRS})
  list (APPEND Pism_EXTERNAL_LIBS ${FFTW_LIBRARIES})
else (FFTW_FOUND)
  add_definitions (-DPISM_HAVE_FFTW=0)
endif (FFTW_FOUND)

add_custom_target (etags
  COMMAND etags `find src -name *.c`
  COMMAND etags -a `find src -name *.cc`
  COMMAND etags -a `find src -name *.h`
  COMMAND etags -a `find src -name *.hh`
  WORKING_DIRECTORY ${Pism_SOURCE_DIR}
)


add_subdirectory (pism-toolkit)
add_subdirectory (src)
add_subdirectory (doc)

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

