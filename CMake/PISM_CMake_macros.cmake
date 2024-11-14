# This file contains CMake macros used in the root CMakeLists.txt

# Set CMake variables to enable rpath
macro(pism_use_rpath)
  ## Use full RPATH, with this setting Pism libraries cannot be moved after installation
  ## but the correct libraries will always be found regardless of LD_LIBRARY_PATH
  ## in use, i.e. don't skip the full RPATH for the build tree
  set (CMAKE_SKIP_BUILD_RPATH FALSE)
  # when building, don't use the install RPATH already
  # (but later on when installing)
  set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  # the RPATH to be used when installing
  set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Mac OS X install_name fix:
  set (CMAKE_MACOSX_RPATH 1)
  set (CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}")
endmacro(pism_use_rpath)

# Set CMake variables to disable rpath
macro(pism_dont_use_rpath)
  set (CMAKE_SKIP_INSTALL_RPATH TRUE)
endmacro(pism_dont_use_rpath)

# Set CMake variables to ensure that everything is static
macro(pism_strictly_static)

  if (BUILD_SHARED_LIBS)
    message(FATAL_ERROR "Please set BUILD_SHARED_LIBS to OFF.")
  endif()

  set (CMAKE_SKIP_RPATH ON CACHE BOOL "Disable RPATH completely")
  set (CMAKE_FIND_LIBRARY_SUFFIXES .a)

  set (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "") # get rid of -rdynamic
  set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "") # ditto

  set_property(GLOBAL PROPERTY LINK_SEARCH_END_STATIC 1)
  set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
  set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)

  pism_dont_use_rpath()
endmacro(pism_strictly_static)

# Set Pism_VERSION_LONG
macro(pism_set_full_version)
  if (EXISTS ${Pism_SOURCE_DIR}/.git)
    # get version information from the repository
    find_program (GIT_EXECUTABLE git DOC "Git executable")
    mark_as_advanced(GIT_EXECUTABLE)
    execute_process (COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      WORKING_DIRECTORY ${Pism_SOURCE_DIR}
      OUTPUT_VARIABLE Pism_COMMIT_HASH
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process (COMMAND ${GIT_EXECUTABLE} --no-pager log -1 "--pretty=format:%an"
      WORKING_DIRECTORY ${Pism_SOURCE_DIR}
      OUTPUT_VARIABLE Pism_COMMIT_AUTHOR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process (COMMAND ${GIT_EXECUTABLE} --no-pager log -1 "--pretty=format:%cs"
      WORKING_DIRECTORY ${Pism_SOURCE_DIR}
      OUTPUT_VARIABLE Pism_COMMIT_DATE
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(Pism_VERSION_LONG "${Pism_VERSION}-${Pism_COMMIT_HASH} committed by ${Pism_COMMIT_AUTHOR} on ${Pism_COMMIT_DATE}")
  else()
    set(Pism_VERSION_LONG "${Pism_VERSION}")
  endif()

  message(STATUS "PISM version: '${Pism_VERSION_LONG}'")
endmacro(pism_set_full_version)

# Set pedantic compiler flags
macro(pism_set_pedantic_flags)
  set (DEFAULT_PEDANTIC_FLAGS "-pedantic -Wall -Wextra -Wno-cast-qual -Wundef -Wshadow -Wpointer-arith -Wno-cast-align -Wwrite-strings -Wno-conversion -Wsign-compare -Wno-redundant-decls -Wno-inline -Wno-long-long -Wmissing-format-attribute -Wpacked -Wdisabled-optimization -Wmultichar -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wendif-labels -Winvalid-pch -Wmissing-field-initializers -Wvariadic-macros -Wstrict-aliasing -Wno-unknown-pragmas -Wno-system-headers")
  set (DEFAULT_PEDANTIC_CFLAGS "${DEFAULT_PEDANTIC_FLAGS} -std=c99")
  set (DEFAULT_PEDANTIC_CXXFLAGS "${DEFAULT_PEDANTIC_FLAGS} -Woverloaded-virtual")
  set (PEDANTIC_CFLAGS ${DEFAULT_PEDANTIC_CFLAGS} CACHE STRING "Compiler flags to enable pedantic warnings")
  set (PEDANTIC_CXXFLAGS ${DEFAULT_PEDANTIC_CXXFLAGS} CACHE STRING "Compiler flags to enable pedantic warnings for C++")
  mark_as_advanced (PEDANTIC_CFLAGS PEDANTIC_CXXFLAGS)
  set (CMAKE_C_FLAGS_DEBUG "-g ${PEDANTIC_CFLAGS}")
  set (CMAKE_CXX_FLAGS_DEBUG "-g ${PEDANTIC_CXXFLAGS}")
endmacro(pism_set_pedantic_flags)

# Make sure that we don't create .petscrc in $HOME, because this would affect
# all PISM runs by the current user.
macro(pism_check_build_dir_location)
  if (DEFINED ENV{HOME})
    # Don't assume that HOME env var is set.
    file (TO_CMAKE_PATH $ENV{HOME} home_dir)
    file (TO_CMAKE_PATH ${PROJECT_BINARY_DIR} build_dir)

    if (${home_dir} STREQUAL ${build_dir})
      message (FATAL_ERROR
        "\n"
        "The build directory is the same as your $HOME!\n"
        "Buiding PISM here would result in a big mess. "
        "Please create a special build directory and run cmake from there.\n")
    endif()
  endif()
endmacro()

function(pism_find_library PREFIX SPEC)
  # Find a library using pkgconfig
  find_package(PkgConfig REQUIRED)
  # Trick pkg-config into using "Libs.private" instead of "Libs" to get the full list of
  # libraries (e.g. not just PETSc, but PETSc and its dependencies). This is needed to
  # build PISM on some systems. See also: https://gitlab.kitware.com/cmake/cmake/-/issues/21714
  if (Pism_PKG_CONFIG_STATIC)
    if (CMAKE_VERSION VERSION_LESS "3.22")
      list(APPEND PKG_CONFIG_EXECUTABLE "--static")
    else()
      set(PKG_CONFIG_ARGN "--static" CACHE INTERNAL "command-line arguments for pkg-config")
    endif()
  endif()
  pkg_search_module(${PREFIX} REQUIRED IMPORTED_TARGET ${SPEC})
  if (${${PREFIX}_FOUND})
    message(STATUS "Found ${PREFIX}: ${${PREFIX}_PREFIX} (found version \"${${PREFIX}_VERSION}\")")
  endif()
endfunction()

macro(pism_find_petsc)
  # Add PETSc directories to CMAKE_PREFIX_PATH (the first on for the
  # "in place" build, the second for the "prefix" build):
  list(APPEND CMAKE_PREFIX_PATH "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")
  list(APPEND CMAKE_PREFIX_PATH "$ENV{PETSC_DIR}")

  message(STATUS
    "Looking for PETSc (PETSC_DIR='$ENV{PETSC_DIR}', PETSC_ARCH='$ENV{PETSC_ARCH}')...")

  pism_find_library(PETSC "PETSc>=3.7.0")
endmacro()

macro(pism_find_prerequisites)

  pism_find_petsc()

  # MPI
  #
  # Note: PISM does not use the MPI's C++ API, but Open MPI may require linking to its
  # library
  find_package (MPI REQUIRED COMPONENTS C CXX)

  # Other required libraries
  pism_find_library (NETCDF "netcdf>=4.4")
  pism_find_library (GSL "gsl>=1.15")
  pism_find_library (FFTW "fftw3>=3.1")

  # UDUNITS does not support pkg-config
  find_package (UDUNITS2)

  # Optional libraries
  if (Pism_USE_PNETCDF)
    pism_find_library(PNETCDF "pnetcdf")
  endif()

  if (Pism_USE_PROJ)
    pism_find_library(PROJ "proj>=6.0")
  endif()

  if (Pism_USE_PIO)
    find_package (ParallelIO REQUIRED)
  endif()

  if (Pism_USE_YAC_INTERPOLATION)
    if (NOT Pism_USE_PROJ)
      message(FATAL_ERROR "Please build PISM with PROJ to use YAC for interpolation")
    endif()

    pism_find_library(YAC "yac-mci>=3.4.0")
    pism_find_library(YAXT "yaxt>=0.11.0")
  endif()

  if (Pism_USE_PARALLEL_NETCDF4)
    # Try to find netcdf_par.h. We assume that NetCDF was compiled with
    # parallel I/O if this header is present.
    find_file(NETCDF_PAR_H netcdf_par.h
      HINTS ${NETCDF_INCLUDE_DIRS} ${NETCDF_INCLUDEDIR}
      NO_DEFAULT_PATH)

    if (NOT NETCDF_PAR_H)
      message(FATAL_ERROR
        "Selected NetCDF library (include: ${NETCDF_INCLUDEDIR}, lib: ${NETCDF_LIBRARIES}) does not support parallel I/O.")
    endif()
  endif()

  if (Pism_USE_JANSSON)
    pism_find_library(JANSSON "jansson>=2.7")
  endif()

endmacro()

macro(pism_set_dependencies)

  # Set include and library directories for *required* libraries.
  include_directories (BEFORE SYSTEM
    ${PETSC_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
    ${UDUNITS2_INCLUDE_DIRS}
    ${NETCDF_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_PATH}
  )

  if (Pism_USE_PROJ)
    include_directories (BEFORE SYSTEM ${PROJ_INCLUDE_DIRS})
  endif()

  if (Pism_USE_PIO)
    include_directories (BEFORE SYSTEM ${ParallelIO_INCLUDES})
  endif()

  if (Pism_USE_PNETCDF)
    include_directories (BEFORE SYSTEM ${PNETCDF_INCLUDE_DIRS})
  endif()

  # Hide distracting CMake variables
  mark_as_advanced(file_cmd MPI_LIBRARY MPI_EXTRA_LIBRARY
    CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT
    MAKE_EXECUTABLE NETCDF_PAR_H)

endmacro()

# Create a list of subdirectories.
# See https://stackoverflow.com/questions/7787823/cmake-how-to-get-the-name-of-all-subdirectories-of-a-directory
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()
