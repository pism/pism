#!/bin/bash

. /opt/intel/oneapi/setvars.sh

set -u
set -e
set -x

yaxt_prefix=${yaxt_prefix:-/opt/yaxt}
build_dir=${build_dir:-/var/tmp/build/yaxt}

mkdir -p ${build_dir}
cd ${build_dir}

yaxt_version=${yaxt_version:-0.11.4}
git clone -b release-${yaxt_version} \
    https://gitlab.dkrz.de/dkrz-sw/yaxt.git

cd yaxt

autoreconf -i

./configure --prefix=${yaxt_prefix} FC=no \
            --with-pic \
            CC=mpiicx \
            CFLAGS="-O3 -fp-model=precise" || cat config.log && echo $PATH

make all && make install
