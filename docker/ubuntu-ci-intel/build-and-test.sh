#!/bin/bash

pism_dir=${pism_dir:?Please set pism_dir}

set -x
set -e
set -u

N=${N:-8}

container_pism_dir=/var/tmp/pism

cmd="cd ${container_pism_dir}/docker/ubuntu-ci-intel/ && ./build.sh && ./run-tests.sh"

docker run --rm -it \
       -v ${pism_dir}:${container_pism_dir} \
       -e N=${N} \
       -e source_dir=${container_pism_dir} \
       ckhrulev/pism-ubuntu-intel:0.2.0 \
       bash -c "${cmd}"
