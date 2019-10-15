configure_file(Dockerfile Dockerfile COPYONLY)

# No image version tag: this Docker image should not be pushed to the Docker Hub.
add_custom_target(docker-archlinux
  COMMAND docker build -t ${DOCKER_USERNAME}/pism-archlinux .
  VERBATIM)
