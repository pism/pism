cmake_minimum_required (VERSION 3.1)
project (Pism C CXX)

if (NOT (${CMAKE_VERSION} VERSION_LESS "3.13"))
  # Policy CMP0078 was introduced in 3.13. See "cmake --help-policies" for details.
  cmake_policy(SET CMP0078 OLD)       #controls names of SWIG targets
endif()

# Require C++11 compiler support.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Looks like CMAKE_CXX_STANDARD does not support Intel C++ compilers
# yet...
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND NOT CMAKE_CXX_FLAGS MATCHES "-std=c\\+\\+11")
  message (STATUS "Adding -std=c++11 to C++ compiler flags for Intel compilers.")
  set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}" CACHE STRING "C++ compiler flags" FORCE)
endif()

include ("CMake/PISM_CMake_macros.cmake")
list (APPEND CMAKE_MODULE_PATH "${Pism_SOURCE_DIR}/CMake")

set (Pism_BRANCH "stable")

# Set Pism_REVISION_TAG
pism_set_revision_tag()

# Put executables in the build directory:
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# Set the install prefix
pism_set_install_prefix()
set(CMAKE_INSTALL_MESSAGE "LAZY")

# Set Pism_DEFAULT_CONFIG_FILE (*after* we set the CMAKE_INSTALL_PREFIX above).
pism_check_build_dir_location()
set (Pism_DEFAULT_CONFIG_FILE "${CMAKE_INSTALL_PREFIX}/${Pism_SHARE_DIR}/pism_config.nc" CACHE STRING "" FORCE)
mark_as_advanced (Pism_DEFAULT_CONFIG_FILE)
file (WRITE ${PROJECT_BINARY_DIR}/.petscrc "-config ${PROJECT_BINARY_DIR}/pism_config.nc")

# The default options cache
option (Pism_BUILD_EXTRA_EXECS "Build extra executables (mostly testing/verification)" OFF)
option (BUILD_SHARED_LIBS      "Build shared Pism libraries" ON)
option (Pism_BUILD_PYTHON_BINDINGS "Build python bindings" OFF)
option (Pism_BUILD_ICEBIN "Build PISM portions of IceBin library" OFF)
option (Pism_BUILD_DOCS "Build PISM's documentation with 'make all'." OFF)
option (Pism_USE_PROJ4 "Use Proj.4 to compute longitudes and latitudes." OFF)
option (Pism_USE_PARALLEL_NETCDF4 "Enables parallel NetCDF-4 I/O." OFF)
option (Pism_USE_PNETCDF "Enables parallel NetCDF-3 I/O using PnetCDF." OFF)
option (Pism_ENABLE_DOCUMENTATION "Enable targets building PISM's documentation." ON)

# PISM will eventually use Jansson to read configuration files.
set (Pism_USE_JANSSON OFF)
# option (Pism_USE_JANSSON "Use Jansson to read configuration files." OFF)

option (Pism_TEST_USING_VALGRIND "Add extra regression tests using valgrind" OFF)
mark_as_advanced (Pism_TEST_USING_VALGRIND)

option (Pism_ADD_FPIC "Add -fPIC to C++ compiler flags (CMAKE_CXX_FLAGS). Try turning it off if it does not work." ON)
option (Pism_CODE_COVERAGE "Add compiler options for code coverage testing." OFF)
option (Pism_LINK_STATICALLY "Set CMake flags to try to ensure that everything is linked statically")
option (Pism_LOOK_FOR_LIBRARIES "Specifies whether PISM should look for libraries. (Disable this on Crays.)" ON)
option (Pism_USE_EVERYTRACE "Use the Everytrace library to provide stacktraces on crashes." OFF)

# Use rpath by default; this has to go first, because rpath settings may be overridden later.
pism_use_rpath()

if (Pism_LINK_STATICALLY)
  pism_strictly_static()
endif ()

# Deal with build types
mark_as_advanced(CLEAR CMAKE_BUILD_TYPE)
if (CMAKE_BUILD_TYPE MATCHES "Debug")
  set (Pism_BUILD_EXTRA_EXECS ON CACHE BOOL "Build extra executables (mostly testing/verification)" FORCE)
  option (Pism_DEBUG              "Enables extra checks in the code."          ON)
  option (Pism_PEDANTIC_WARNINGS  "Compile with pedantic warnings."            ON)
  option (Pism_GPROF_FLAGS        "Add flags necessary to profile with gprof." OFF)
endif (CMAKE_BUILD_TYPE MATCHES "Debug")

# Enable extra checks in the code:
if (Pism_DEBUG)
  message (STATUS "Adding -DPISM_DEBUG=1 to compiler flags.")
  add_definitions (-DPISM_DEBUG=1)
else()
  message (STATUS "Adding -DPISM_DEBUG=0 to compiler flags.")
  add_definitions (-DPISM_DEBUG=0)
endif ()

# Add -fPIC to C and CXX flags.
if (Pism_ADD_FPIC)
  if (NOT CMAKE_C_FLAGS MATCHES "-fPIC")
    message (STATUS "Adding -fPIC to C compiler flags.")
    set (CMAKE_C_FLAGS "-fPIC ${CMAKE_C_FLAGS}" CACHE STRING "C compiler flags" FORCE)
  endif()

  if (NOT CMAKE_CXX_FLAGS MATCHES "-fPIC")
    message (STATUS "Adding -fPIC to C++ compiler flags.")
    set (CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}" CACHE STRING "C++ compiler flags" FORCE)
  endif()
endif ()

if (Pism_CODE_COVERAGE)
  if (NOT CMAKE_CXX_FLAGS MATCHES "-fprofile-arcs")
    message (STATUS "Adding -fprofile-arcs -ftest-coverage to C++ compiler flags.")
    set (CMAKE_CXX_FLAGS "-fprofile-arcs -ftest-coverage -g -O0 ${CMAKE_CXX_FLAGS}" CACHE STRING "C++ compiler flags" FORCE)
  endif()
  add_custom_target (coverage_report
    # remove coverage data from src/pythonbindings
    COMMAND lcov --directory ${Pism_BINARY_DIR}/src/pythonbindings -z
    COMMAND lcov --base-directory ${Pism_SOURCE_DIR} --directory src --quiet --no-external --capture -o pism-coverage.info
    COMMAND genhtml -t "PISM Coverage report" -o cover --demangle-cpp --legend pism-coverage.info
    WORKING_DIRECTORY ${Pism_BINARY_DIR}
    VERBATIM
    )
  add_custom_target (coverage_reset
    COMMAND lcov -d ${Pism_BINARY_DIR} -z
    WORKING_DIRECTORY ${Pism_BINARY_DIR}
    VERBATIM
    )
endif ()

if (Pism_PEDANTIC_WARNINGS)
  pism_set_pedantic_flags()
endif (Pism_PEDANTIC_WARNINGS)

if (Pism_GPROF_FLAGS)
  set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
endif ()

if (Pism_USE_JANSSON)
  add_definitions (-DPISM_USE_JANSSON=1)
endif()

# Look for libraries using find_package(...), etc. Run CMake with -DPism_LOOK_FOR_LIBRARIES=OFF
# to build on systems that rely on the module system to set all compiler and linker flags.
if (Pism_LOOK_FOR_LIBRARIES)
  pism_find_prerequisites()
endif()

# Set Pism_EXTERNAL_LIBS and include directories.
pism_set_dependencies()

# Make sure that PetscScalar is double (not complex<double>.)
pism_check_petsc_scalar_type()
pism_set_version_info()

if (Pism_USE_EVERYTRACE)
  find_package(Everytrace REQUIRED)
  list (APPEND Pism_EXTERNAL_LIBS ${EVERYTRACE_LIBRARY})
endif()

if (Pism_BUILD_PYTHON_BINDINGS)
  find_package(Python REQUIRED)
  find_package(PETSc4Py REQUIRED)
  find_package(SWIG REQUIRED)

  if (DEFINED PETSC4PY_VERSION)
    # FindPETSc4Py.cmake does not put PETSC4PY_VERSION into the CMake cache,
    # so we save it here.
    set(Pism_PETSC4PY_VERSION ${PETSC4PY_VERSION} CACHE STRING "PETSc4Py version")
    mark_as_advanced(Pism_PETSC4PY_VERSION)
  endif()

  mark_as_advanced (SWIG_DIR SWIG_EXECUTABLE SWIG_VERSION)

  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}'))"
    OUTPUT_VARIABLE PISM_INSTALL_PYTHON_MODULE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(PISM_INSTALL_PYTHON_MODULE_DIR ${PISM_INSTALL_PYTHON_MODULE_DIR}
    CACHE PATH "Python extension module installation directory." )
endif ()

# Use PROJ.4 to compute latitudes and longitudes of grid points (if available).
if (Pism_USE_PROJ4)
  add_definitions (-DPISM_USE_PROJ4=1)
else()
  add_definitions (-DPISM_USE_PROJ4=0)
endif()

add_custom_target (etags
  COMMAND find . -regex ".*\\.(c|cc|h|hh)" | xargs etags --class-qualify --no-defines
  WORKING_DIRECTORY ${Pism_SOURCE_DIR}
  VERBATIM
)

# re-run tests that failed
add_custom_target (retest
  COMMAND ${CMAKE_CTEST_COMMAND} --rerun-failed --output-on-failure
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  VERBATIM
)

# run Python tests
add_custom_target (test-python
  COMMAND ${CMAKE_CTEST_COMMAND} -R "Python:"
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  VERBATIM
)

# Install helper scripts residing in util/
install(DIRECTORY ${PROJECT_SOURCE_DIR}/util/ DESTINATION ${Pism_BIN_DIR}
  USE_SOURCE_PERMISSIONS
  FILES_MATCHING PATTERN "*.py")


install(DIRECTORY ${PROJECT_SOURCE_DIR}/examples
  DESTINATION ${Pism_SHARE_DIR}
  USE_SOURCE_PERMISSIONS)

add_subdirectory (src)
add_subdirectory (site-packages)

if (Pism_ENABLE_DOCUMENTATION)
  if (Pism_BUILD_DOCS)
    add_subdirectory (doc)
  else()
    add_subdirectory (doc EXCLUDE_FROM_ALL)
  endif()
endif()

# PISM regression testing
ENABLE_TESTING()
include(CTest)
add_subdirectory (test)
add_subdirectory (test/regression)
