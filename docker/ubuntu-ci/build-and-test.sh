#!/bin/bash

pism_dir=${pism_dir:?Please set pism_dir}

set -x
set -e

if [[ -z ${CC} || -z ${CXX} ]];
then
  compiler=${compiler:-gcc}

  if [ ${compiler} == "gcc" ];
  then
    # available versions: 9 10 11 12 13 14
    version=${version:-14}
    CC=gcc-${version}
    CXX=g++-${version}
  elif [ ${compiler} == "clang" ];
  then
    # available versions: 14 15 16 17 18
    version=${version:-18}
    CC=clang-${version}
    CXX=clang++-${version}
  fi
fi

set -u

N=${N:-8}

# The Make target to build
target=${target:-all}

container_pism_dir=/home/builder/project

cmd="cd ${container_pism_dir}/docker/ubuntu-ci/ && ./build.sh"

if [ ${target} == "all" ];
then
  cmd="${cmd} && ./run-tests.sh"
fi

docker run --rm -it \
       -v ${pism_dir}:${container_pism_dir} \
       -e CC=${CC} \
       -e CXX=${CXX} \
       -e N=${N} \
       -e source_dir=${container_pism_dir} \
       -e target="${target}" \
       ckhrulev/pism-ubuntu:0.1.13 \
       bash -c "${cmd}"
