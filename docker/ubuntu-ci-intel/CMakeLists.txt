foreach(file Dockerfile petsc.sh)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${file} . COPYONLY)
endforeach()

foreach(script hdf5.sh netcdf.sh pnetcdf.sh yac.sh)
  configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/code/${script} . COPYONLY)
endforeach()

set(PISM_DOCKER_UBUNTU_INTEL_VERSION 0.2.0)

add_custom_target(docker-ubuntu-intel
  COMMAND docker build -t ${DOCKER_USERNAME}/pism-ubuntu-intel:${PISM_DOCKER_UBUNTU_INTEL_VERSION} .
  VERBATIM)
