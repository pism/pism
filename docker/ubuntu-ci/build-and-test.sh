#!/bin/bash

set -x
set -e

if [[ -z ${CC} || -z ${CXX} ]];
then
  compiler=${compiler:-gcc}
  version=${version:-14}

  if [ ${compiler} == "gcc" ];
  then
    CC=gcc-${version}
    CXX=g++-${version}
  elif [ ${compiler} == "clang" ];
  then
    CC=clang-${version}
    CXX=clang++-${version}
  fi
fi

set -u

N=${N:-8}

pism_dir=${pism_dir:?Please set pism_dir}

container_pism_dir=/home/builder/project

docker run --rm -it \
       -v ${pism_dir}:${container_pism_dir} \
       -e CC=${CC} \
       -e CXX=${CXX} \
       -e N=${N} \
       -e source_dir=${container_pism_dir} \
       ckhrulev/pism-ubuntu:0.1.12 \
       bash -c "cd ${container_pism_dir}/docker/ubuntu-ci/ && ./build.sh && ./run-tests.sh"
