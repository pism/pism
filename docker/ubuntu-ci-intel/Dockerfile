FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

RUN <<EOF
    echo "Set up the repository"

    set -x
    set -e
    set -u

    apt-get update
    apt-get install -y --no-install-recommends wget gpg sudo ca-certificates

# See https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2024-1/apt.html
#
# download the key to system keyring
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
    | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

# add signed entry to apt sources and configure the APT client to use Intel repository:
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

    apt-get update
    apt-get upgrade -y

    rm -rf /var/lib/apt/lists/*
EOF

RUN <<EOF
    echo "Install required oneAPI packages (development)"

    apt-get update

    apt-get install -y --no-install-recommends \
    intel-oneapi-compiler-dpcpp-cpp-2024.1 \
    intel-oneapi-mpi-devel-2021.12 \
    ""

    rm -rf /var/lib/apt/lists/*
EOF

RUN <<EOF
    echo "Install PISM's prerequisite tools and libraries"
    apt-get update

    apt-get install -y --no-install-recommends \
    cmake \
    libfftw3-dev \
    libgsl-dev \
    libproj-dev \
    libudunits2-dev \
    m4 \
    make \
    pkg-config \
    python3 \
    ""

    rm -rf /var/lib/apt/lists/*
EOF

# See
# https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2024-1/use-the-setvars-and-oneapi-vars-scripts-with-linux.html
#
# for the right way to initialize Intel oneAPI.

COPY hdf5.sh /tmp/
RUN <<EOF
    echo "Install HDF5"
    . /opt/intel/oneapi/setvars.sh

    set -eux
    export MPICC=mpiicx
    export build_dir=/var/tmp/build/hdf5

    /tmp/hdf5.sh
    rm -rf ${build_dir}
    rm -rf /opt/hdf5/share/hdf5_examples
EOF

COPY netcdf.sh /tmp/
RUN <<EOF
    echo "Install NetCDF"
    . /opt/intel/oneapi/setvars.sh

    set -eux
    export MPICC=mpiicx
    export build_dir=/var/tmp/build/netcdf

    /tmp/netcdf.sh
    rm -rf ${build_dir}
EOF

COPY petsc.sh /tmp/
RUN <<EOF
    echo "Install PETSc"
    . /opt/intel/oneapi/setvars.sh

    set -eux
    export MPICC=mpiicx
    export MPICXX=mpiicpx
    export build_dir=/var/tmp/build/petsc
    export COPTFLAGS="-O3 -fp-model=precise -march=native -mtune=native"
    export version=3.19.6

    /tmp/petsc.sh

    rm -rf ${build_dir}
    rm -rf /opt/petsc/share/petsc
EOF

ENV PATH="/opt/hdf5/bin:/opt/netcdf/bin:/opt/petsc/bin:$PATH"

RUN ln -s /opt/intel/oneapi/setvars.sh /etc/profile.d/02-intel-oneapi-setvars.sh

# Add a non-privileged user
RUN useradd --create-home --system --shell=/bin/false worker && usermod --lock worker
USER worker

# Important parts of ". /opt/intel/oneapi/setvars.sh"
ENV ONEAPI_ROOT=/opt/intel/oneapi
ENV LD_LIBRARY_PATH=${ONEAPI_ROOT}/tbb/2021.12/env/../lib/intel64/gcc4.8:${ONEAPI_ROOT}/mpi/2021.12/opt/mpi/libfabric/lib:${ONEAPI_ROOT}/mpi/2021.12/lib:${ONEAPI_ROOT}/mkl/2024.1/lib:${ONEAPI_ROOT}/compiler/2024.1/opt/compiler/lib:${ONEAPI_ROOT}/compiler/2024.1/lib:${ONEAPI_ROOT}/mpi/2021.12/opt/mpi/libfabric/lib:${ONEAPI_ROOT}/mpi/2021.12/lib:${ONEAPI_ROOT}/compiler/2024.1/opt/compiler/lib:${ONEAPI_ROOT}/compiler/2024.1/lib
ENV I_MPI_ROOT=${ONEAPI_ROOT}/mpi/2021.12
ENV FI_PROVIDER_PATH=${ONEAPI_ROOT}/mpi/2021.12/opt/mpi/libfabric/lib/prov:/usr/lib/x86_64-linux-gnu/libfabric
ENV PATH=${ONEAPI_ROOT}/mpi/2021.12/bin:${ONEAPI_ROOT}/mkl/2024.1/bin/:${ONEAPI_ROOT}/compiler/2024.1/bin:${PATH}
