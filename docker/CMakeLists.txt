find_program(CUT_EXECUTABLE cut)
find_program(XARGS_EXECUTABLE xargs)
mark_as_advanced(CUT_EXECUTABLE XARGS_EXECUTABLE)

execute_process(
  COMMAND ${CUT_EXECUTABLE} -f1 -d, ${Pism_SOURCE_DIR}/doc/sphinx/installation/debian-packages.csv
  COMMAND ${XARGS_EXECUTABLE}
  OUTPUT_VARIABLE PISM_DEBIAN_PACKAGE_LIST)
mark_as_advanced(PISM_DEBIAN_PACKAGE_LIST)

add_subdirectory(debian)
add_subdirectory(archlinux)
