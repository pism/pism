add_custom_target(docker-pism-base
  docker build --tag ${Pism_CONTAINER_TAG}-base --target=pism-base ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
  WORKING_DIRECTORY ${Pism_BINARY_DIR})

add_custom_target(docker-pism-runtime
  docker build --tag ${Pism_CONTAINER_TAG}-runtime --target=pism-runtime ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
  WORKING_DIRECTORY ${Pism_BINARY_DIR})
