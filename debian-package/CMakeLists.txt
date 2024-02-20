# CPack setup

set (CPACK_PACKAGE_NAME "pism")
set (CPACK_PACKAGE_CONTACT "uaf-pism@alaska.edu")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Parallel Ice Sheet Model (PISM)")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/description.txt")
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set (CPACK_DEBIAN_PACKAGE_DEPENDS "python3")

# "libjs-mathjax"

set (CPACK_GENERATOR "DEB")

# use all available CPUs
set (CPACK_THREADS "0")

set (CPACK_DEBIAN_PACKAGE_NAME "pism")

set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")

set (CPACK_DEBIAN_PACKAGE_MAINTAINER "PISM developers at UAF <uaf-pism@alaska.edu>")

set (CPACK_DEBIAN_PACKAGE_SECTION "science")

set (CPACK_DEBIAN_PACKAGE_PRIORITY "optional")

set (CPACK_STRIP_FILES YES)

# work-around for CPack's lack of proper ldconfig support
# https://gitlab.kitware.com/cmake/cmake/-/issues/21834
set (CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_SOURCE_DIR}/triggers")
set (CPACK_DEBIAN_PACKAGE_CONTROL_STRICT_PERMISSION TRUE)

include(CPack)
