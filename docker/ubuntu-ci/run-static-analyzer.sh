#!/bin/bash

set -e
set -u
set -x

pism_dir=${pism_dir:?Please set pism_dir}

container_pism_dir=/home/worker/project

cmd="cd ${container_pism_dir}/docker/ubuntu-ci/ && ./static-analyzer.sh"

docker run --rm -it \
       -v ${pism_dir}:${container_pism_dir} \
       -e N=${N:-8} \
       -e source_dir=${container_pism_dir} \
       ckhrulev/pism-ubuntu:0.1.17 \
       bash -c "${cmd}"
