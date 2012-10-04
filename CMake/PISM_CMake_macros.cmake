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
  set (CMAKE_SKIP_RPATH ON CACHE BOOL "Disable RPATH completely")
  set (CMAKE_FIND_LIBRARY_SUFFIXES .a)
  set (BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared Pism libraries" FORCE)
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "") # get rid of -rdynamic
  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "") # ditto
  set_property(GLOBAL PROPERTY LINK_SEARCH_END_STATIC 1) # get rid of -Bdynamic

  pism_dont_use_rpath()
endmacro(pism_strictly_static)

# Set CMake variables appropriate for building a .deb package
macro(pism_build_debian_package)
  set (Pism_BUILD_TYPE "Release" CACHE STRING "PISM build type" FORCE)
  set (CMAKE_INSTALL_PREFIX "/usr" CACHE STRING "Install prefix" FORCE)
  set (Pism_BUILD_DOCS OFF CACHE BOOL "Build PISM documentation" FORCE)
  set (Pism_BUILD_BROWSER OFF CACHE BOOL "Build PISM source code browsers" FORCE)
  set (Pism_BUILD_EXTRA_EXECS OFF CACHE BOOL "Build extra executables (mostly testing/verification)" FORCE)

  # RPATH handling
  pism_dont_use_rpath()
endmacro(pism_build_debian_package)

# Set the revision tag if PISM was checked out using Git.
macro(pism_set_revision_tag_git)
  if (NOT Pism_VERSION)
    if (EXISTS ${Pism_SOURCE_DIR}/.git)
      find_program (GIT_EXECUTABLE git DOC "Git executable")
      mark_as_advanced(GIT_EXECUTABLE)
      execute_process (COMMAND ${GIT_EXECUTABLE} describe --always --match v*
        WORKING_DIRECTORY ${Pism_SOURCE_DIR}
        OUTPUT_VARIABLE Pism_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif (EXISTS ${Pism_SOURCE_DIR}/.git)
  endif(NOT Pism_VERSION)
endmacro(pism_set_revision_tag_git)

# Set the revision tag if PISM was checked out using Subversion.
macro(pism_set_revision_tag_svn)
  if (NOT Pism_VERSION)
    if (EXISTS ${Pism_SOURCE_DIR}/.svn)
      find_package(Subversion)
      if (SUBVERSION_FOUND)
        Subversion_WC_INFO(${Pism_SOURCE_DIR}/src "Pism")
        set(Pism_VERSION "${Pism_WC_LAST_CHANGED_DATE}")
      endif(SUBVERSION_FOUND)
    endif(EXISTS ${Pism_SOURCE_DIR}/.svn)
  endif(NOT Pism_VERSION)
endmacro(pism_set_revision_tag_svn)

# Set the PISM revision tag
macro(pism_set_revision_tag)
  # Git
  pism_set_revision_tag_git()

  # Subversion
  pism_set_revision_tag_svn()

  # Otherwise...
  if (NOT Pism_VERSION)
    set (Pism_VERSION "unknown")
  endif (NOT Pism_VERSION)

  set (Pism_REVISION_TAG "${Pism_BRANCH} ${Pism_VERSION}")
endmacro(pism_set_revision_tag)

# Set pedantic compiler flags
macro(pism_set_pedantic_flags)
  set (DEFAULT_PEDANTIC_CFLAGS "-pedantic -Wall -Wextra -Wno-cast-qual -Wundef -Wshadow -Wpointer-arith -Wno-cast-align -Wwrite-strings -Wno-conversion -Wsign-compare -Wno-redundant-decls -Winline -Wno-long-long -Wmissing-format-attribute -Wmissing-noreturn -Wpacked -Wdisabled-optimization -Wmultichar -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wendif-labels -Winvalid-pch -Wmissing-field-initializers -Wvariadic-macros -Wstrict-aliasing -funit-at-a-time")
  set (DEFAULT_PEDANTIC_CXXFLAGS "${DEFAULT_PEDANTIC_CFLAGS} -Woverloaded-virtual")
  set (PEDANTIC_CFLAGS ${DEFAULT_PEDANTIC_CFLAGS} CACHE STRING "Compiler flags to enable pedantic warnings")
  set (PEDANTIC_CXXFLAGS ${DEFAULT_PEDANTIC_CXXFLAGS} CACHE STRING "Compiler flags to enable pedantic warnings for C++")
  mark_as_advanced (PEDANTIC_CFLAGS PEDANTIC_CXXFLAGS)
  set (CMAKE_C_FLAGS_DEBUG "-g ${PEDANTIC_CFLAGS}")
  set (CMAKE_CXX_FLAGS_DEBUG "-g ${PEDANTIC_CXXFLAGS}")
endmacro(pism_set_pedantic_flags)

# Make sure that we don't create .petscrc in $HOME, because this would affect
# all PISM runs by the current user.
macro(pism_check_build_dir_location)
  file (TO_CMAKE_PATH $ENV{HOME} home_dir)
  file (TO_CMAKE_PATH ${PROJECT_BINARY_DIR} build_dir)

  if (${home_dir} STREQUAL ${build_dir})
    message (FATAL_ERROR
      "\n"
      "The build directory is the same as your $HOME!\n"
      "Buiding PISM here would result in a big mess. "
      "Please create a special build directory and run cmake from there.\n")
  endif()
endmacro()
