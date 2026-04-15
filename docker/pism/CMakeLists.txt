set (Pism_CONTAINER_TAG "pism/pism" CACHE STRING
  "Tag to use when building the PISM container using Intel oneAPI" )

add_custom_target(docker-pism
  docker build --file ${Pism_SOURCE_DIR}/docker/pism/Dockerfile --tag ${Pism_CONTAINER_TAG} --build-arg=BASE=pism/pism-base --build-arg=RUNTIME=pism/pism-runtime ${Pism_SOURCE_DIR}
  VERBATIM
  WORKING_DIRECTORY ${Pism_SOURCE_DIR})

add_dependencies(docker-pism docker-pism-base docker-pism-runtime)

add_subdirectory(base)
