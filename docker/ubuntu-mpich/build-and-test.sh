#!/bin/bash

pism_dir=${pism_dir:?Please set pism_dir}

set -x
set -e
set -u

N=${N:-8}

# The Make target to build
target=${target:-all}

container_pism_dir=/home/builder/project

cmd="cd ${container_pism_dir}/docker/ubuntu-mpich/ && ./build.sh"

if [ ${target} == "all" ];
then
  cmd="${cmd} && ./run-tests.sh"
fi

docker run --rm -it \
       --shm-size=2g \
       -v ${pism_dir}:${container_pism_dir}:ro \
       -e N=${N} \
       -e source_dir=${container_pism_dir} \
       -e target="${target}" \
       pism/pism-mpich:latest \
       bash -c "${cmd}"
