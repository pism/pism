configure_file(Dockerfile.in Dockerfile)

foreach(script hdf5.sh netcdf.sh pnetcdf.sh parallelio.sh)
  configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/code/${script} . COPYONLY)
endforeach()

set(PISM_DOCKER_DEBIAN_VERSION 0.1.2)

add_custom_target(docker-debian
  COMMAND docker build -t ${DOCKER_USERNAME}/pism-debian:${PISM_DOCKER_DEBIAN_VERSION} .
  VERBATIM)
