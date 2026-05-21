cmake_minimum_required (VERSION 3.20)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Get PISM's "short" X.Y.Z version from a file
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION)
  file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION Pism_VERSION LIMIT_COUNT 1)
else()
  message(FATAL_ERROR "'VERSION' not found: cannot determine PISM's version")
endif()

project (Pism
  VERSION ${Pism_VERSION}
  HOMEPAGE_URL "https://www.pism.io"
  LANGUAGES C CXX)

include ("CMake/PISM_CMake_macros.cmake")
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

pism_set_full_version()

# Set install directories for scikit-build-core (pip install .):
if (SKBUILD)
  set(CMAKE_INSTALL_BINDIR      "${SKBUILD_SCRIPTS_DIR}"      CACHE PATH "" FORCE)
  set(CMAKE_INSTALL_LIBDIR      "${SKBUILD_DATA_DIR}/lib"     CACHE PATH "" FORCE)
  set(CMAKE_INSTALL_DATAROOTDIR "${SKBUILD_DATA_DIR}/share"   CACHE PATH "" FORCE)
  set(CMAKE_INSTALL_DATADIR     "${SKBUILD_DATA_DIR}/share"   CACHE PATH "" FORCE)
  set(CMAKE_INSTALL_INCLUDEDIR  "${SKBUILD_DATA_DIR}/include" CACHE PATH "" FORCE)
  set(PISM_PYTHON_MODULE_DIR    "${SKBUILD_PLATLIB_DIR}"      CACHE PATH "" FORCE)
endif()

include (GNUInstallDirs)
# doc directory: override capitalization ("pism" instead of "Pism" from the "project(...)"
# call above)
set(CMAKE_INSTALL_DOCDIR ${CMAKE_INSTALL_DATAROOTDIR}/doc/pism)

# Deal with build types
mark_as_advanced(CLEAR CMAKE_BUILD_TYPE)
# Set the build type to "Release" so that users get optimized code by default:
#
# This code block comes from the Kokkos project.
if (NOT CMAKE_BUILD_TYPE)
  set(DEFAULT_BUILD_TYPE "Release")
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: Debug, Release, RelWithDebInfo and MinSizeRel."
    FORCE)
endif()

# ``find_package()`` uses ``<PackageName>_ROOT`` variables.
cmake_policy(SET CMP0074 NEW)

# Require C++11 compiler support.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add a flag that makes it easier to use LLDB to debug code compiled with Clang (if
# supported by C and C++ compilers).
include(CheckCompilerFlag)
check_compiler_flag(C -fstandalone-debug PISM_C_FSTANDALONE_DEBUG)
check_compiler_flag(CXX -fstandalone-debug PISM_CXX_FSTANDALONE_DEBUG)
if (PISM_C_FSTANDALONE_DEBUG AND PISM_CXX_FSTANDALONE_DEBUG)
  add_compile_options(-fstandalone-debug)
  message(STATUS "Adding -fstandalone-debug to compiler flags")
endif()

# Put executables in the build directory:
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

set(CMAKE_INSTALL_MESSAGE "LAZY")
set(CMAKE_COLOR_MAKEFILE OFF)

pism_check_build_dir_location()

# Set Pism_CONFIG_FILE (*after* we set the CMAKE_INSTALL_PREFIX above).
if (Pism_CONFIG_FILE)
  set (Pism_CONFIG_FILE "${Pism_CONFIG_FILE}" CACHE STRING "" FORCE)
else()
  set(Pism_CONFIG_FILE "${CMAKE_INSTALL_FULL_DATADIR}/pism/pism_config.nc" CACHE STRING "" FORCE)
  if (SKBUILD)
    if (DEFINED ENV{CONDA_PREFIX})
      message(STATUS "CONDA_PREFIX is set: setting Pism_CONFIG_FILE accordingly")
      set(Pism_CONFIG_FILE "$ENV{CONDA_PREFIX}/share/pism/pism_config.nc" CACHE STRING "" FORCE)
    elseif (DEFINED ENV{VIRTUAL_ENV})
      message(STATUS "VIRTUAL_ENV is set: setting Pism_CONFIG_FILE accordingly")
      set(Pism_CONFIG_FILE "$ENV{VIRTUAL_ENV}/share/pism/pism_config.nc" CACHE STRING "" FORCE)
    endif()
  endif(SKBUILD)
endif()
mark_as_advanced (Pism_CONFIG_FILE)

# The default options cache
option (BUILD_SHARED_LIBS      "Build shared Pism libraries" ON)
option (Pism_DEBIAN_SYSTEMWIDE "Use settings appropriate when installing in a system-wide location on Debian/Ubuntu or building a package" OFF)
option (Pism_BUILD_PYTHON_BINDINGS "Build python bindings" OFF)
option (Pism_BUILD_ICEBIN "Build PISM portions of IceBin library" OFF)
option (Pism_BUILD_DOCS "Build PISM's documentation with 'make all'." OFF)
option (Pism_USE_PROJ "Use PROJ to compute longitudes and latitudes." OFF)
option (Pism_USE_YAC "Use YAC (for interpolation and async output)" OFF)
option (Pism_USE_PARALLEL_NETCDF4 "Enables parallel NetCDF-4 I/O." OFF)
option (Pism_USE_PNETCDF "Enables parallel NetCDF-3 I/O using PnetCDF." OFF)
option (Pism_ENABLE_DOCUMENTATION "Enable targets building PISM's documentation." ON)

option (Pism_TEST_USING_VALGRIND "Add extra regression tests using valgrind" OFF)
mark_as_advanced (Pism_TEST_USING_VALGRIND)

option (Pism_PKG_CONFIG_STATIC "Use the --static pkg-config flag to get linker flags needed by prerequisites" OFF)
option (Pism_ADD_FPIC "Add -fPIC to C & C++ compiler flags. Try turning it off if it does not work." ON)
option (Pism_CODE_COVERAGE "Add compiler options for code coverage testing." OFF)
option (Pism_LINK_STATICALLY "Set CMake flags to try to ensure that everything is linked statically")
option (Pism_LOOK_FOR_LIBRARIES "Specifies whether PISM should look for libraries. (Disable this on Crays.)" ON)
option (Pism_USE_EVERYTRACE "Use the Everytrace library to provide stacktraces on crashes." OFF)
option (Pism_DEBUG "Enables extra checks in the code." OFF)
option (Pism_PEDANTIC_WARNINGS  "Compile with pedantic warnings." ON)
option (Pism_GPROF_FLAGS        "Add flags necessary to profile with gprof." OFF)

# Use rpath by default; this has to go first, because rpath settings may be overridden later.
if (Pism_DEBIAN_SYSTEMWIDE)
  pism_dont_use_rpath()
elseif (SKBUILD)
  # Use CMAKE_INSTALL_RPATH and CMAKE_BUILD_RPATH_USE_ORIGIN passed via pyproject.toml.
  #
  # Avoid baking the wheel-staging absolute libdir into the install RPATH or the macOS
  # install_name_dir, which would break the binary once pip moves files out of the staging
  # area.
  set (CMAKE_SKIP_BUILD_RPATH FALSE)
  set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set (CMAKE_MACOSX_RPATH 1)
else()
  pism_use_rpath()
endif()

if (Pism_LINK_STATICALLY)
  pism_strictly_static()
endif ()

# Add -fPIC to C and CXX flags.
if (Pism_ADD_FPIC)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

if (Pism_CODE_COVERAGE)
  add_compile_options(-fprofile-arcs -ftest-coverage -g -O0)
  add_link_options(--coverage)

  add_custom_target (coverage_report
    # remove coverage data from src/pythonbindings
    COMMAND lcov --directory ${Pism_BINARY_DIR}/src/pythonbindings -z
    COMMAND lcov --directory ${Pism_BINARY_DIR}/src/external -z
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

# Look for libraries using find_package(...), etc. Run CMake with -DPism_LOOK_FOR_LIBRARIES=OFF
# to build on systems that rely on the module system to set all compiler and linker flags.
if (Pism_LOOK_FOR_LIBRARIES)
  pism_find_prerequisites()
endif()

# Set include directories.
pism_set_dependencies()

if (Pism_USE_EVERYTRACE)
  find_package(Everytrace REQUIRED)
endif()

if (Pism_BUILD_PYTHON_BINDINGS)
  find_package(Python3 COMPONENTS Interpreter Development)
  find_package(PETSc4Py REQUIRED)
  find_package(SWIG REQUIRED)

  if (DEFINED PETSC4PY_VERSION)
    # FindPETSc4Py.cmake does not put PETSC4PY_VERSION into the CMake cache,
    # so we save it here.
    set(Pism_PETSC4PY_VERSION ${PETSC4PY_VERSION} CACHE STRING "PETSc4Py version")
    mark_as_advanced(Pism_PETSC4PY_VERSION)
  endif()

  mark_as_advanced (SWIG_DIR SWIG_EXECUTABLE SWIG_VERSION)

  if (NOT PISM_PYTHON_MODULE_DIR)
    set (pism_python_path_py "${CMAKE_CURRENT_BINARY_DIR}/pism_python_path.py")
    if (Pism_DEBIAN_SYSTEMWIDE)
      set (pism_sysconfig_scheme "deb_system")
    else()
      set (pism_sysconfig_scheme "posix_prefix")
    endif()
    # Note: this is almost the same as Python3_SITEARCH provided by CMake's FindPython3,
    # but relative to the selected install prefix (CMAKE_INSTALL_PREFIX) and including a
    # special case needed to build the .deb package.
    file (WRITE "${pism_python_path_py}" "
# Get module installation path from Python:
import sysconfig
print(sysconfig.get_path('platlib', '${pism_sysconfig_scheme}', vars={'platbase': '${CMAKE_INSTALL_PREFIX}'}))
")
    unset(pism_sysconfig_scheme)

    execute_process(
      COMMAND ${Python3_EXECUTABLE} ${pism_python_path_py}
      OUTPUT_VARIABLE PISM_PROBED_PYTHON_MODULE_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    # convert to a path relative to the install prefix
    file(RELATIVE_PATH PISM_PYTHON_MODULE_DIR
      ${CMAKE_INSTALL_PREFIX} ${PISM_PROBED_PYTHON_MODULE_DIR})
    unset(pism_python_path_py)
    unset(PISM_PROBED_PYTHON_MODULE_DIR)
  endif(NOT PISM_PYTHON_MODULE_DIR)
  mark_as_advanced(PISM_PYTHON_MODULE_DIR)
endif (Pism_BUILD_PYTHON_BINDINGS)

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
install(PROGRAMS
  util/pism_adjust_timeline
  util/pism_async_writer
  util/pism_check_stationarity
  util/pism_create_timeline
  util/pism_fill_missing
  util/pism_fill_missing_petsc
  util/pism_flowline
  util/pism_ismip7_writer
  util/pism_nc2cdo
  util/pism_nccmp
  util/pism_plot_profiling
  DESTINATION ${CMAKE_INSTALL_BINDIR})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/examples
  DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/pism
  DIRECTORY_PERMISSIONS
  OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
  PATTERN .gitignore EXCLUDE
  PATTERN __pycache__ EXCLUDE
  PATTERN "*.pyc" EXCLUDE
  PATTERN "*.nc" EXCLUDE
  PATTERN "*.npz" EXCLUDE
  PATTERN "*~" EXCLUDE)

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
list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")
list(REMOVE_DUPLICATES CMAKE_CTEST_ARGUMENTS)
enable_testing()
include(CTest)

add_subdirectory (test)
add_subdirectory (test/regression)

add_subdirectory (docker)

add_subdirectory (package)
