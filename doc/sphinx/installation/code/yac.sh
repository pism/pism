#!/bin/bash

set -u
set -e
set -x

build_dir=${build_dir:-.}

rm -rf ${build_dir}/yaxt
rm -rf ${build_dir}/yac

mkdir -p ${build_dir}
cd ${build_dir}

# manual-begin
export CC=${CC:-mpicc}
export FC=${FC:-mpifort}
export CFLAGS="-O3 -g -march=native"

prefix=${prefix:-$HOME/local/yac}

yaxt_version=0.11.3
git clone -b release-${yaxt_version} \
    https://gitlab.dkrz.de/dkrz-sw/yaxt.git

cd yaxt

autoreconf -i

./configure --prefix=${prefix} \
            --with-pic

make all && make install

cd -

yac_version=3.5.2
git clone -b release-${yac_version} \
    https://gitlab.dkrz.de/dkrz-sw/yac.git

cd yac

test -f ./configure || ./autogen.sh

./configure --prefix=${prefix} \
            --with-yaxt-root=${prefix} \
            --disable-netcdf \
            --disable-examples \
            --disable-tools \
            --disable-deprecated \
            --with-pic

make all && make install
# manual-end
