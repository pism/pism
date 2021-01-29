configure_file(Dockerfile.in Dockerfile)

foreach(script petsc.sh hdf5.sh netcdf.sh pnetcdf.sh parallelio.sh)
  configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/code/${script} . COPYONLY)
endforeach()

set(PISM_DOCKER_UBUNTU_VERSION 0.1.2)

add_custom_target(docker-ubuntu
  COMMAND docker build -t ${DOCKER_USERNAME}/pism-ubuntu:${PISM_DOCKER_UBUNTU_VERSION} .
  VERBATIM)
