configure_file(Dockerfile.in Dockerfile)

add_custom_target(docker-debian
  COMMAND docker build -t ckhrulev/pism-debian .
  VERBATIM)
