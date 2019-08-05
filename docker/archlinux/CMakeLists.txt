configure_file(Dockerfile.in Dockerfile)
configure_file(${Pism_SOURCE_DIR}/doc/sphinx/installation/arch-install-yay.sh arch-install-yay.sh)

add_custom_target(docker-archlinux
  COMMAND docker build -t ckhrulev/pism-arch .
  VERBATIM)
