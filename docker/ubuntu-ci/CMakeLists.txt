configure_file(Dockerfile.in Dockerfile)

foreach(script petsc.sh hdf5.sh netcdf.sh pnetcdf.sh parallelio.sh)
  configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/code/${script} . COPYONLY)
endforeach()

configure_file(cmake-3_7_0.sh . COPYONLY)
configure_file(old-petsc.sh . COPYONLY)

set(PISM_DOCKER_UBUNTU_VERSION 0.1.7)

add_custom_target(docker-ubuntu
  COMMAND docker build -t ${DOCKER_USERNAME}/pism-ubuntu:${PISM_DOCKER_UBUNTU_VERSION} .
  VERBATIM)
