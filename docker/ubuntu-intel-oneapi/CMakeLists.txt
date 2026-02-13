set (Pism_CONTAINER_TAG "pism/pism" CACHE STRING
  "Tag to use when building the PISM container using Intel oneAPI" )

add_custom_target(docker-intel-oneapi
  docker build --file ${Pism_SOURCE_DIR}/docker/ubuntu-intel-oneapi/Dockerfile --tag ${Pism_CONTAINER_TAG} ${Pism_SOURCE_DIR}
  VERBATIM
  WORKING_DIRECTORY ${Pism_SOURCE_DIR})
