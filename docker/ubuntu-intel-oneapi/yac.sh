#!/bin/bash

. /opt/intel/oneapi/setvars.sh

set -u
set -e
set -x

yaxt_prefix=${yaxt_prefix:-/opt/yaxt}
yac_prefix=${yac_prefix:-/opt/yac}
build_dir=${build_dir:-/var/tmp/build/yac}

mkdir -p ${build_dir}
cd ${build_dir}

yac_version=${yac_version:-3.13.2}
git clone -b release-${yac_version} \
    https://gitlab.dkrz.de/dkrz-sw/yac.git

cd yac

test -f ./configure || ./autogen.sh

./configure --prefix=${yac_prefix} \
            --with-yaxt-root=${yaxt_prefix} \
            --disable-netcdf \
            --disable-examples \
            --disable-tools \
            --disable-deprecated \
            --disable-fortran-bindings \
            --enable-python-bindings \
            --with-pic \
            --with-external-lapack=mkl \
            --with-mkl-root=${MKLROOT} \
            LDFLAGS="-qmkl=sequential" \
	    CC=mpiicx \
	    CFLAGS="-O3 -fp-model=precise"

make all && make install
