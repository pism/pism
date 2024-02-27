configure_file(Dockerfile . COPYONLY)

configure_file(${Pism_SOURCE_DIR}/package/debian/build.sh . COPYONLY)

set(PISM_DOCKER_UBUNTU_DEB_VERSION 0.1.1)

add_custom_target(docker-ubuntu-deb
  COMMAND docker build -t ${DOCKER_USERNAME}/pism-ubuntu-deb-package:${PISM_DOCKER_UBUNTU_DEB_VERSION} .
  VERBATIM)
