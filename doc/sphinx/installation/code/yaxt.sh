#!/bin/bash

set -u
set -e
set -x

rm -rf yaxt

# manual-begin
prefix=$HOME/local/yac
version=0.11.1

git clone -b release-${version} https://gitlab.dkrz.de/dkrz-sw/yaxt.git

cd yaxt

autoreconf -i

./configure --prefix=${prefix} \
            --with-pic

make all && make install
# manual-end
