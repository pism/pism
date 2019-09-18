configure_file(Dockerfile.in Dockerfile)
configure_file(build_netcdf.sh build_netcdf.sh COPYONLY)

add_custom_target(docker-debian
  COMMAND docker build -t ckhrulev/pism-debian .
  VERBATIM)
