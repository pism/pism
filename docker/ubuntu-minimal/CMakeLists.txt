configure_file(Dockerfile.in Dockerfile)

foreach(script install_libraries.sh pism.sh)
  configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/code/${script} . COPYONLY)
endforeach()

set(PISM_DOCKER_UBUNTU_MIN_VERSION 0.1.0)

add_custom_target(docker-ubuntu-minimal
  COMMAND docker build -t ${DOCKER_USERNAME}/pism-ubuntu-minimal:${PISM_DOCKER_UBUNTU_MIN_VERSION} .
  VERBATIM)
