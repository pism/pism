#!/bin/bash

set -e
set -u
set -x

source_dir=${source_dir:-.}

cmake -S ${source_dir} -B build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/usr \
      -DPism_DEBIAN_SYSTEMWIDE=YES \
      -DPism_BUILD_DOCS=YES \
      -DPism_BUILD_EXTRA_EXECS=YES \
      -DPism_BUILD_PYTHON_BINDINGS=NO \
      -DPism_USE_PROJ=YES

make --no-print-directory -C build -j ${N:-4} package
