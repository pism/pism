find_program(CUT_EXECUTABLE cut)
find_program(SED_EXECUTABLE sed)
find_program(XARGS_EXECUTABLE xargs)
find_program(GREP_EXECUTABLE grep)

mark_as_advanced(CUT_EXECUTABLE SED_EXECUTABLE XARGS_EXECUTABLE GREP_EXECUTABLE)

# Get required prerequisites from the manual and remove ones that will be replaced by
# custom versions in Docker images. In particular:
#
# - libx11-dev is not needed
# - petsc-dev uses OpenMPI (it is easier to use MPICH for development)
# - petsc-dev uses OpenMPI-based parallel HDF5 (we build our own HDF5)
# - libnetcdf-dev is replaced by a parallel version built from sources
execute_process(
  COMMAND ${GREP_EXECUTABLE} "required" ${Pism_SOURCE_DIR}/doc/sphinx/installation/debian-packages.csv
  COMMAND ${CUT_EXECUTABLE} -f1 -d,
  COMMAND ${XARGS_EXECUTABLE}
  COMMAND ${SED_EXECUTABLE} "s/``//g"
  COMMAND ${SED_EXECUTABLE} "s/libx11-dev//"
  COMMAND ${SED_EXECUTABLE} "s/libnetcdf-dev//"
  COMMAND ${SED_EXECUTABLE} "s/petsc-dev//"
  OUTPUT_VARIABLE PISM_DEBIAN_PACKAGE_LIST
  OUTPUT_STRIP_TRAILING_WHITESPACE)
mark_as_advanced(PISM_DEBIAN_PACKAGE_LIST)

set(DOCKER_USERNAME "ckhrulev" CACHE STRING "Docker user name to use when building images")
mark_as_advanced(DOCKER_USERNAME)

add_subdirectory(ubuntu-ci)
add_subdirectory(ubuntu-ci-intel)
add_subdirectory(ubuntu-deb-package)
add_subdirectory(archlinux)
