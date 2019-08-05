
execute_process(
  COMMAND ${CUT_EXECUTABLE} -f1 -d, ${Pism_SOURCE_DIR}/doc/sphinx/installation/debian-packages.csv
  COMMAND ${XARGS_EXECUTABLE}
  OUTPUT_VARIABLE PISM_DEBIAN_PACKAGE_LIST)
mark_as_advanced(PISM_DEBIAN_PACKAGE_LIST)

configure_file(Dockerfile.in Dockerfile)

add_custom_target(docker-debian
  COMMAND docker build -t ckhrulev/pism-debian .
  VERBATIM)
