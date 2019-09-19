configure_file(Dockerfile.in Dockerfile)
configure_file(build_netcdf.sh build_netcdf.sh COPYONLY)

set(PISM_DOCKER_DEBIAN_VERSION 0.1.0)

add_custom_target(docker-debian
  COMMAND docker build -t ckhrulev/pism-debian:${PISM_DOCKER_DEBIAN_VERSION} .
  VERBATIM)
