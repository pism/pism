foreach(file Dockerfile hdf5.sh netcdf.sh petsc.sh)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${file} . COPYONLY)
endforeach()

set(PISM_DOCKER_UBUNTU_INTEL_VERSION 0.2.0)

add_custom_target(docker-ubuntu-intel
  COMMAND docker build -t ${DOCKER_USERNAME}/pism-ubuntu-intel:${PISM_DOCKER_UBUNTU_INTEL_VERSION} .
  VERBATIM)
