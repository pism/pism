set (Pism_MPICH_CONTAINER_TAG "pism/pism-mpich" CACHE STRING
  "Tag to use when building the PISM+MPICH container" )

add_custom_target(docker-ubuntu-mpich
  docker build --file ${Pism_SOURCE_DIR}/docker/ubuntu-mpich/Dockerfile --tag ${Pism_MPICH_CONTAINER_TAG} ${Pism_SOURCE_DIR}
  VERBATIM
  WORKING_DIRECTORY ${Pism_SOURCE_DIR})
