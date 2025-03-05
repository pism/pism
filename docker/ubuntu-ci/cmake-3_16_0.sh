#!/bin/bash

set -e
set -u
set -x

# Install CMake 3.16 (the oldest supported version).

version=3.16.0
build_dir=${build_dir:-/var/tmp/build/cmake-${version}}
prefix=${prefix:-~/local/cmake-${version}}

mkdir -p ${build_dir}
cd ${build_dir}

wget -nc \
     https://github.com/Kitware/CMake/releases/download/v${version}/cmake-${version}.tar.gz
rm -rf cmake-${version}
tar xzf cmake-${version}.tar.gz

cd cmake-${version}

./configure --prefix=${prefix} -- -DCMAKE_USE_OPENSSL=NO

make -j all
make install

# clean up
cd ~ && rm -rf ${build_dir}
