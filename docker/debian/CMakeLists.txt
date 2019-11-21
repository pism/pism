configure_file(Dockerfile.in Dockerfile)
configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/code/hdf5.sh . COPYONLY)
configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/code/netcdf.sh . COPYONLY)
configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/code/pnetcdf.sh . COPYONLY)
configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/code/parallelio.sh . COPYONLY)

set(PISM_DOCKER_DEBIAN_VERSION 0.1.2)

add_custom_target(docker-debian
  COMMAND docker build -t ${DOCKER_USERNAME}/pism-debian:${PISM_DOCKER_DEBIAN_VERSION} .
  VERBATIM)
