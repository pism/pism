FROM ubuntu:rolling
ARG DEBIAN_FRONTEND=noninteractive

RUN <<EOF
echo "Install packages using apt-get"

    packages="
    autoconf
    automake
    build-essential
    ca-certificates
    cdo
    cmake
    cython3
    g++
    gfortran
    git
    libfftw3-dev
    libfyaml-dev
    libgsl-dev
    libmpich-dev
    libproj-dev
    libpython3.13-dev
    libtool
    libtool
    libudunits2-dev
    make
    nco
    pkg-config
    python3
    python3-dev
    python3-pip
    python3-setuptools
    python3-venv
    swig
    wget
    "

    apt-get update
    apt-get install --no-install-recommends -y ${packages}
EOF

RUN useradd --create-home --system --shell=/bin/false builder && usermod --lock builder
USER builder

ARG BUILD_SCRIPTS=doc/sphinx/installation/code/
ARG HOME=/home/builder
ARG PREFIX=$HOME/local
# Create a venv for Python stuff that we cannot install using apt-get:
RUN python3 -m venv --system-site-packages $PREFIX

RUN <<EOF
    echo "Install mpi4py, netCDF4, NumPy, SciPy"
    set -eux
    . $PREFIX/bin/activate
    python3 -m pip install "mpi4py>=4.0.0" netCDF4 numpy scipy
EOF

# Build the current PETSc version
COPY ${BUILD_SCRIPTS}/petsc.sh /tmp/
RUN <<EOF
    echo "Install PETSc"
    set -eux
    . $PREFIX/bin/activate
    cd /tmp
    prefix=$PREFIX petsc4py=1 ./petsc.sh
    rm -rf ~/local/build
EOF

COPY ${BUILD_SCRIPTS}/yac.sh /tmp/
RUN <<EOF
    echo "Install YAXT and YAC"
    set -eux
    . $PREFIX/bin/activate
    cd /tmp
    prefix=$PREFIX ./yac.sh
    rm -rf yaxt yac
EOF

# Install parallel HDF5, NetCDF and PnetCDF by building them from sources.
#
# Uses wget installed above.
#
# This step uses scripts that are used to describe how to install these libraries in the
# manual. Re-using them here makes it easier to keep this part of the manual up to date.

COPY ${BUILD_SCRIPTS}/hdf5.sh /tmp/
RUN <<EOF
    echo "Install HDF5"
    set -eux
    . $PREFIX/bin/activate
    cd /tmp
    prefix=$PREFIX ./hdf5.sh
    rm -rf ~/local/build
EOF

COPY ${BUILD_SCRIPTS}/pnetcdf.sh /tmp/
RUN <<EOF
    echo "Install PnetCDF"
    set -eux
    . $PREFIX/bin/activate
    cd /tmp
    prefix=$PREFIX ./pnetcdf.sh
    rm -rf ~/local/build
EOF

COPY ${BUILD_SCRIPTS}/netcdf.sh /tmp/
RUN <<EOF
    echo "Install NetCDF"
    set -eux
    . $PREFIX/bin/activate
    cd /tmp
    prefix=$PREFIX ./netcdf.sh
    rm -rf ~/local/build
EOF

ENV PETSC_DIR="$PREFIX"
ENV PYTHONPATH="$PREFIX/lib:$PREFIX/lib/python3.13/site-packages"
ENV PATH=$PREFIX/bin:$PATH
