#!/bin/bash

# Build and test using all available compilers

set -e
set -u
set -x

export compiler=gcc

for v in 9 10 11 12 13 14;
do
  version=${v} ./build-and-test.sh
done

export compiler=clang

for v in 14 15 16 17 18;
do
  version=${v} ./build-and-test.sh
done
