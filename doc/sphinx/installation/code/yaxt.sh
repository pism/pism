#!/bin/bash

set -e
set -u
set -x

# Install YAXT 0.9.0 in ~/local/yaxt,
# using ~/local/build/yaxt as the build directory.

prefix=$HOME/local/yaxt
version=0.9.0
build_dir=~/local/build/yaxt
CFLAGS="-std=gnu99 -O3"
url=https://www.dkrz.de/redmine/attachments/498/yaxt-0.9.0.tar.gz

mkdir -p ${build_dir}
pushd ${build_dir}

wget -nc ${url}
tar xzf yaxt-${version}.tar.gz

pushd yaxt-${version}

./configure --prefix=${prefix} --enable-shared --enable-static CFLAGS="-std=gnu99 -O3"
make
make install

popd
popd