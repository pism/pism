#!/usr/bin/env bash
# Compile solistis.c -> libsolstis.dylib (Apple clang + OpenMP).
set -euo pipefail
cd "$(dirname "$0")"

P="${CONDA_PREFIX:-/Users/andy/miniforge3/envs/pism}"

/usr/bin/cc -O3 -ffast-math -fno-finite-math-only -std=c11 -shared -fPIC \
    -Xpreprocessor -fopenmp -I"$P/include" \
    solstis.c \
    -L"$P/lib" -lomp -Wl,-rpath,"$P/lib" \
    -o libsolstis.dylib

echo "built libsolstis.dylib"
