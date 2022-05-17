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
  set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}")
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Mac OS X install_name fix:
  set (CMAKE_MACOSX_RPATH 1)
  set (CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}")
endmacro(pism_use_rpath)

# Set CMake variables to disable rpath
macro(pism_dont_use_rpath)
  set (CMAKE_SKIP_BUILD_RPATH TRUE)
  set (CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}")
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
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

# Set the revision tag if PISM was checked out using Git.
macro(pism_set_revision_tag_git)
  if (NOT Pism_VERSION)
    if (EXISTS ${Pism_SOURCE_DIR}/.git)
      find_program (GIT_EXECUTABLE git DOC "Git executable")
      mark_as_advanced(GIT_EXECUTABLE)
      if (${Pism_BRANCH} MATCHES "stable")
        execute_process (COMMAND ${GIT_EXECUTABLE} describe --always --match v?.?*
          WORKING_DIRECTORY ${Pism_SOURCE_DIR}
          OUTPUT_VARIABLE Pism_VERSION
          OUTPUT_STRIP_TRAILING_WHITESPACE)
      else()
        execute_process (COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
          WORKING_DIRECTORY ${Pism_SOURCE_DIR}
          OUTPUT_VARIABLE Pism_VERSION
          OUTPUT_STRIP_TRAILING_WHITESPACE)
      endif()
      execute_process (COMMAND ${GIT_EXECUTABLE} --no-pager log -1 "--pretty=format:committed by %an on %ci"
        WORKING_DIRECTORY ${Pism_SOURCE_DIR}
        OUTPUT_VARIABLE Pism_COMMIT_INFO
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      set(Pism_VERSION "${Pism_VERSION} ${Pism_COMMIT_INFO}")
    endif (EXISTS ${Pism_SOURCE_DIR}/.git)
  endif(NOT Pism_VERSION)
endmacro(pism_set_revision_tag_git)

# Set the PISM revision tag
macro(pism_set_revision_tag)
  # Git
  pism_set_revision_tag_git()

  # Otherwise...
  if (NOT Pism_VERSION)
    set (Pism_VERSION "no-version-control")
  endif (NOT Pism_VERSION)

  set (Pism_REVISION_TAG "${Pism_BRANCH} ${Pism_VERSION}")

  message(STATUS "Configuring PISM version '${Pism_REVISION_TAG}'")
endmacro(pism_set_revision_tag)

macro(pism_set_install_prefix)
  # Allow setting a custom install prefix using the PISM_INSTALL_PREFIX environment variable.
  string (LENGTH "$ENV{PISM_INSTALL_PREFIX}" INSTALL_PREFIX_LENGTH)
  if (INSTALL_PREFIX_LENGTH)
    set (CMAKE_INSTALL_PREFIX $ENV{PISM_INSTALL_PREFIX} CACHE PATH "PISM install prefix" FORCE)
    message (STATUS "Setting PISM install prefix to ${CMAKE_INSTALL_PREFIX}.")
  endif()

  # Define the directory structure.
  set (Pism_BIN_DIR "bin")
  set (Pism_LIB_DIR "lib")
  set (Pism_SHARE_DIR "share/pism")
  set (Pism_DOC_DIR "share/doc/pism")
endmacro()

# Set pedantic compiler flags
macro(pism_set_pedantic_flags)
  set (DEFAULT_PEDANTIC_FLAGS "-pedantic -Wall -Wextra -Wno-cast-qual -Wundef -Wshadow -Wpointer-arith -Wno-cast-align -Wwrite-strings -Wno-conversion -Wsign-compare -Wno-redundant-decls -Wno-inline -Wno-long-long -Wmissing-format-attribute -Wpacked -Wdisabled-optimization -Wmultichar -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wendif-labels -Winvalid-pch -Wmissing-field-initializers -Wvariadic-macros -Wstrict-aliasing -funit-at-a-time -Wno-unknown-pragmas")
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

macro(pism_find_prerequisites)
  find_package(PkgConfig REQUIRED)
  set(OLD_PKG_CONFIG_PATH $ENV{PKG_CONFIG_PATH})

  # PETSc
  # set root of location to find PETSc's pkg-config

  set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
  set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

  pkg_search_module(PETSC REQUIRED IMPORTED_TARGET "PETSc>=3.7.0")
  message(STATUS "Found PETSc ${PETSC_VERSION} in ${PETSC_PREFIX}")
  # restore old PKG_CONFIG_PATH
  set(ENV{PKG_CONFIG_PATH} ${OLD_PKG_CONFIG_PATH})
  unset(OLD_PKG_CONFIG_PATH)

  # MPI
  find_package (MPI REQUIRED COMPONENTS C)

  # Other required libraries
  find_package (UDUNITS2 REQUIRED)
  find_package (GSL REQUIRED)
  find_package (NetCDF REQUIRED)
  find_package (FFTW REQUIRED)
  find_package (HDF5 COMPONENTS C HL)

  # Optional libraries
  if (Pism_USE_PNETCDF)
    find_package (PNetCDF REQUIRED)
  endif()

  if (Pism_USE_PROJ)
    find_package (PROJ REQUIRED)
  endif()

  if (Pism_USE_PIO)
    find_package (ParallelIO REQUIRED)
  endif()

  if (Pism_USE_PARALLEL_NETCDF4)
    # Try to find netcdf_par.h. We assume that NetCDF was compiled with
    # parallel I/O if this header is present.
    find_file(NETCDF_PAR_H netcdf_par.h HINTS ${NETCDF_INCLUDES} NO_DEFAULT_PATH)

    # Set default values for build options
    if (NOT NETCDF_PAR_H)
      message(FATAL_ERROR
        "Selected NetCDF library (include: ${NETCDF_INCLUDES}, lib: ${NETCDF_LIBRARIES}) does not support parallel I/O.")
    endif()
  endif()

  if (Pism_USE_JANSSON)
    find_package(Jansson REQUIRED)

    if (NOT JANSSON_FOUND)
      set(pism_jansson_dir ${Pism_BINARY_DIR}/jansson)
      include(ExternalProject)
      ExternalProject_Add(pism_jansson
        GIT_REPOSITORY https://github.com/akheron/jansson.git
        GIT_TAG 2.7
        TIMEOUT 10
        PREFIX ${Pism_BINARY_DIR} # install with PISM
        INSTALL_DIR ${pism_jansson_dir}
        CMAKE_ARGS -DJANSSON_BUILD_DOCS=OFF -DCMAKE_INSTALL_PREFIX=${pism_jansson_dir}
        LOG_DOWNLOAD ON
        LOG_BUILD ON
        LOG_CONFIGURE ON
        )
      set(JANSSON_INCLUDE_DIRS ${pism_jansson_dir}/include CACHE STRING "Jansson include directory" OFRCE)
      set(JANSSON_LIBRARIES "-L${pism_jansson_dir}/lib -ljansson" CACHE STRING "Jansson library" FORCE)
      set(Pism_BUILD_JANSSON ON CACHE BOOL "ON if we are using our own Jansson build." FORCE)
      message(WARNING "
Jansson was not found.
We will try to download and build it automatically. If it does not work, please install it manually and try again.
")
    endif()
  endif()

endmacro()

macro(pism_set_dependencies)

  # Set include and library directories for *required* libraries.
  #
  # Note: PISM does not use HDF5 directly, but we still need to be able to include hdf5.h
  # to record its version.
  include_directories (BEFORE
    ${PETSC_INCLUDE_DIRS}
    ${FFTW_INCLUDES}
    ${GSL_INCLUDES}
    ${UDUNITS2_INCLUDES}
    ${HDF5_INCLUDE_DIRS}
    ${NETCDF_INCLUDES}
    ${MPI_C_INCLUDE_PATH})

  # Use option values to set compiler and linker flags
  set (Pism_EXTERNAL_LIBS "")

  # required libraries
  list (APPEND Pism_EXTERNAL_LIBS
    PkgConfig::PETSC
    ${UDUNITS2_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${GSL_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_C_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES})

  # optional libraries
  if (Pism_USE_JANSSON)
    include_directories (${JANSSON_INCLUDE_DIRS})
    list (APPEND Pism_EXTERNAL_LIBS ${JANSSON_LIBRARIES})
  endif()

  if (Pism_USE_PROJ)
    include_directories (${PROJ_INCLUDES})
    list (APPEND Pism_EXTERNAL_LIBS ${PROJ_LIBRARIES})
  endif()

  if (Pism_USE_PIO)
    include_directories (${ParallelIO_INCLUDES})
    list (APPEND Pism_EXTERNAL_LIBS ${ParallelIO_LIBRARIES})
  endif()

  if (Pism_USE_PNETCDF)
    include_directories (${PNETCDF_INCLUDES})
    list (APPEND Pism_EXTERNAL_LIBS ${PNETCDF_LIBRARIES})
  endif()

  # Hide distracting CMake variables
  mark_as_advanced(file_cmd MPI_LIBRARY MPI_EXTRA_LIBRARY
    HDF5_C_LIBRARY_dl HDF5_C_LIBRARY_hdf5 HDF5_C_LIBRARY_hdf5_hl HDF5_C_LIBRARY_m HDF5_C_LIBRARY_z
    CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT
    MAKE_EXECUTABLE HDF5_DIR NETCDF_PAR_H)

endmacro()

macro (pism_petsc_get_variable name var)
  set (pism_petsc_config_makefile "${PROJECT_BINARY_DIR}/Makefile.pism_petsc")
  file (WRITE "${pism_petsc_config_makefile}"
"## This file was autogenerated by FindPETSc.cmake
# PETSC_DIR  = ${PETSC_DIR}
# PETSC_ARCH = ${PETSC_ARCH}
include ${petsc_conf_rules}
include ${petsc_conf_variables}
show :
\t-@echo -n \${\${VARIABLE}}
")
  set (${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
  execute_process (COMMAND ${MAKE_EXECUTABLE} --no-print-directory -f ${pism_petsc_config_makefile} show VARIABLE=${name}
    OUTPUT_VARIABLE ${var}
    RESULT_VARIABLE petsc_return)
  string(CONFIGURE "\${${var}}" ${var} ESCAPE_QUOTES)
  file (REMOVE ${pism_petsc_config_makefile})
endmacro (pism_petsc_get_variable)

# Make sure that PetscScalar is double and not complex.
macro(pism_check_petsc_scalar_type)
  pism_petsc_get_variable("PETSC_SCALAR" PISM_PETSC_SCALAR)

  if (${PISM_PETSC_SCALAR} MATCHES "complex")
    message(FATAL_ERROR
      "PETSc configured with --with-scalar-type=complex cannot be used to build PISM.")
  endif()
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
