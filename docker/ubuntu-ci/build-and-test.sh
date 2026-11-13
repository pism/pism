#!/bin/bash

pism_dir=${pism_dir:?Please set pism_dir}

set -x
set -e

if [[ -z ${CC} || -z ${CXX} ]];
then
  compiler=${compiler:-gcc}

  if [ ${compiler} == "gcc" ];
  then
    # available versions: 11 12 13 14 15 16
    version=${version:-16}
    CC=gcc-${version}
    CXX=g++-${version}
  elif [ ${compiler} == "clang" ];
  then
    # available versions: 17 18 19 20 21 22
    version=${version:-22}
    CC=clang-${version}
    CXX=clang++-${version}
  fi
fi

set -u

N=${N:-8}

# The Make target to build
target=${target:-all}

# set "minimal" to a non-zero-length string to build the minimal PISM configuration
minimal=${minimal:-}

# set "run_analyzer" to run Clang's static analyzer
run_analyzer=${run_analyzer:-}

container_pism_dir=/home/builder/project

cmd="cd ${container_pism_dir}/docker/ubuntu-ci/ && ./build.sh"

if [ ${target} == "all" && ! -z ${run_analyzer} ];
then
  cmd="${cmd} && ./run-tests.sh"
fi

docker run --rm -it \
       --shm-size=2g \
       -v ${pism_dir}:${container_pism_dir}:ro \
       -e CC=${CC} \
       -e CXX=${CXX} \
       -e N=${N} \
       -e source_dir=${container_pism_dir} \
       -e target="${target}" \
       -e minimal="${minimal}" \
       -e run_analyzer="${run_analyzer}" \
       ckhrulev/pism-ubuntu:0.1.19 \
       bash -c "${cmd}"
