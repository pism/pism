

# CPack setup

set (CPACK_PACKAGE_NAME "pism")
set (CPACK_PACKAGE_CONTACT "uaf-pism@alaska.edu")

set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set (CPACK_DEBIAN_PACKAGE_DEPENDS "python3")
# "libjs-mathjax"
set (CPACK_GENERATOR "DEB")

# use all available CPUs
set (CPACK_THREADS "0")

set (CPACK_DEBIAN_PACKAGE_NAME "pism")

set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")

set (CPACK_DEBIAN_PACKAGE_MAINTAINER "PISM developers at UAF <uaf-pism@alaska.edu>")

set (CPACK_DEBIAN_PACKAGE_DESCRIPTION "Parallel Ice Sheet Model (PISM)")

set (CPACK_DEBIAN_PACKAGE_SECTION "science")

set (CPACK_DEBIAN_PACKAGE_PRIORITY "optional")

set (CPACK_STRIP_FILES YES)

set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_SOURCE_DIR}/triggers")

include(CPack)
