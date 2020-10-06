find_program(CUT_EXECUTABLE cut)
find_program(SED_EXECUTABLE sed)
find_program(XARGS_EXECUTABLE xargs)
mark_as_advanced(CUT_EXECUTABLE SED_EXECUTABLE XARGS_EXECUTABLE)

execute_process(
  COMMAND ${CUT_EXECUTABLE} -f1 -d, ${Pism_SOURCE_DIR}/doc/sphinx/installation/debian-packages.csv
  COMMAND ${XARGS_EXECUTABLE}
  COMMAND ${SED_EXECUTABLE} "s/``//g;s/libx11-dev//;s/ncview//;s/cdo//"
  OUTPUT_VARIABLE PISM_DEBIAN_PACKAGE_LIST
  OUTPUT_STRIP_TRAILING_WHITESPACE)
mark_as_advanced(PISM_DEBIAN_PACKAGE_LIST)

set(DOCKER_USERNAME "ckhrulev" CACHE STRING "Docker user name to use when building images")
mark_as_advanced(DOCKER_USERNAME)

add_subdirectory(ubuntu)
add_subdirectory(archlinux)
