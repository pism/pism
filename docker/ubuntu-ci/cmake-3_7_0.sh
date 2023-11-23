#!/bin/bash

set -e
set -u
set -x

# Install CMake 3.7 (the oldest supported version).

version=3.7.0
build_dir=${build_dir:-/var/tmp/build/cmake-${version}}
prefix=${prefix:-~/local/cmake-${version}}

mkdir -p ${build_dir}
cd ${build_dir}

wget -nc \
     https://github.com/Kitware/CMake/releases/download/v${version}/cmake-${version}.tar.gz
rm -rf cmake-${version}
tar xzf cmake-${version}.tar.gz

cd cmake-${version}

# Use GCC 10 to build old CMake: it requires modifications (adding #include <limits>) to
# build with more recent versions.
export CC=gcc-10
export CXX=g++-10
./configure --prefix=${prefix}

make -j8 all
make install
