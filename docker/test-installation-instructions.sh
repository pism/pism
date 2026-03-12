#!/bin/bash

# Tests if PISM in ${pism_dir} can be built on current Ubuntu after installing packages
# listed in the manual.

set -x
set -e
set -u

pism_dir=${pism_dir:?Please set pism_dir}

container_pism_dir=/var/tmp/pism

packages=`grep required ${pism_dir}/doc/sphinx/installation/debian-packages.csv | cut -f1 -d, | sed 's/\`//g' | xargs`

cmd="
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends ${packages}

git config --global --add safe.directory ${container_pism_dir}

cmake -S ${container_pism_dir} -B /tmp/pism-build

make -j -C /tmp/pism-build all
"

docker run --rm -it \
       -v ${pism_dir}:${container_pism_dir} \
       -e source_dir=${container_pism_dir} \
       -e packages="${packages}" \
       ubuntu:rolling \
       bash -x -u -e -c "${cmd}"
