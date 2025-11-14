FROM ubuntu:24.04 AS intel-runtime
ARG DEBIAN_FRONTEND=noninteractive
ARG ONEAPI_VERSION=2025.3
ARG IMPI_VERSION=2021.17

RUN <<EOF
    echo "Set up the Intel oneAPI repository"

    set -x
    set -e
    set -u

    apt-get update
    apt-get install -y --no-install-recommends wget gpg sudo ca-certificates

# See https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2025-0/apt-005.html
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
    echo "Install required Intel oneAPI packages (runtime)"
    apt-get update

    apt-get install -y --no-install-recommends \
    intel-oneapi-compiler-shared-runtime-${ONEAPI_VERSION} \
    intel-oneapi-compiler-fortran-runtime-${ONEAPI_VERSION} \
    intel-oneapi-mkl-core-${ONEAPI_VERSION} \
    intel-oneapi-mpi-${IMPI_VERSION} \
    ""

    rm -rf /var/lib/apt/lists/*
EOF

FROM intel-runtime AS pism-runtime

RUN <<EOF
    echo "Install PISM's runtime dependencies (ones from packages)"
    apt-get update

    apt-get install -y --no-install-recommends \
    libfftw3-double3 \
    libgsl27 \
    libproj25 \
    libudunits2-0 \
    libfyaml0 \
    ""

    rm -rf /var/lib/apt/lists/*
EOF

FROM pism-runtime AS intel-devel

RUN <<EOF
    echo "Install required Intel oneAPI packages (development)"

    apt-get update

    apt-get install -y --no-install-recommends \
    intel-oneapi-compiler-dpcpp-cpp-${ONEAPI_VERSION} \
    intel-oneapi-compiler-fortran-${ONEAPI_VERSION} \
    intel-oneapi-mkl-core-devel-${ONEAPI_VERSION} \
    intel-oneapi-mpi-devel-${IMPI_VERSION} \
    ""

    rm -rf /var/lib/apt/lists/*
EOF

FROM intel-devel AS pism-devel

RUN <<EOF
    echo "Install PISM's build tools and build dependencies"
    apt-get update

    apt-get install -y --no-install-recommends \
    autoconf \
    automake \
    cmake \
    coreutils \
    git \
    libfftw3-dev \
    libfyaml-dev \
    libgsl-dev \
    libproj-dev \
    libtool \
    libudunits2-dev \
    m4 \
    make \
    pkg-config \
    python3 \
    ssh \
    ""

    rm -rf /var/lib/apt/lists/*
EOF

# See
# https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2025-1/use-the-setvars-and-oneapi-vars-scripts-with-linux.html
#
# for the right way to initialize Intel oneAPI.
COPY docker/ubuntu-intel-oneapi/hdf5.sh /tmp/
RUN <<EOF
    echo "Install HDF5"

    set -eux
    export build_dir=/var/tmp/build/hdf5
    export prefix=/opt/hdf5

    bash /tmp/hdf5.sh

    rm -rf ${build_dir}
    rm -rf /opt/hdf5/share/hdf5_examples
EOF

COPY docker/ubuntu-intel-oneapi/netcdf.sh /tmp/
RUN <<EOF
    echo "Install NetCDF"

    set -eux
    export hdf5_prefix=/opt/hdf5
    export build_dir=/var/tmp/build/netcdf
    export prefix=/opt/netcdf

    bash /tmp/netcdf.sh
    rm -rf ${build_dir}
EOF

COPY docker/ubuntu-intel-oneapi/yaxt.sh /tmp/
RUN <<EOF
    echo "Install YAXT"

    set -eux
    bash /tmp/yaxt.sh

    rm -rf /var/tmp/build/yaxt
EOF

COPY docker/ubuntu-intel-oneapi/yac.sh /tmp/
RUN <<EOF
    echo "Install YAC"

    set -eux
    bash /tmp/yac.sh

    rm -rf /var/tmp/build/yac
EOF

COPY docker/ubuntu-intel-oneapi/petsc.sh /tmp/
RUN <<EOF
    echo "Install PETSc"

    set -eux
    export build_dir=/var/tmp/build/petsc
    export prefix=/opt/petsc

    bash /tmp/petsc.sh

    rm -rf ${build_dir}
    rm -rf /opt/petsc/share/petsc
EOF

COPY docker/ubuntu-intel-oneapi/pism.sh /tmp/
# Copy the Docker context to /var/tmp/pism-src
COPY . /var/tmp/pism-src
RUN <<EOF
    echo "Install PISM"

    set -eux
    export src_dir=/var/tmp/pism-src
    export build_dir=/var/tmp/build/pism
    export prefix=/opt/pism

    bash /tmp/pism.sh

    rm -rf ${build_dir}
    rm -rf ${prefix}/include ${prefix}/share/doc/pism
EOF

FROM pism-runtime AS pism-container

# copy PISM's runtime dependencies that were built from sources:
COPY --from=pism-devel /opt/hdf5 /opt/hdf5/
COPY --from=pism-devel /opt/netcdf /opt/netcdf/
COPY --from=pism-devel /opt/petsc /opt/petsc/
COPY --from=pism-devel /opt/pism /opt/pism/
COPY --from=pism-devel /opt/yac /opt/yac/
COPY --from=pism-devel /opt/yaxt /opt/yaxt/

ENV PATH="/opt/hdf5/bin:/opt/netcdf/bin:/opt/petsc/bin:/opt/pism/bin:$PATH"

RUN ln -s /opt/intel/oneapi/setvars.sh /etc/profile.d/02-intel-oneapi-setvars.sh

# Add a non-privileged user
RUN useradd --create-home --system --shell=/bin/false worker && usermod --lock worker
USER worker

# Important parts of ". /opt/intel/oneapi/setvars.sh"
ENV ONEAPI_ROOT=/opt/intel/oneapi
ENV LD_LIBRARY_PATH=${ONEAPI_ROOT}/tbb/${IMPI_VERSION}/env/../lib/intel64/gcc4.8:${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/opt/mpi/libfabric/lib:${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/lib:${ONEAPI_ROOT}/mkl/${ONEAPI_VERSION}/lib:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/opt/compiler/lib:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/lib:${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/opt/mpi/libfabric/lib:${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/lib:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/opt/compiler/lib:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/lib
ENV I_MPI_ROOT=${ONEAPI_ROOT}/mpi/${IMPI_VERSION}
ENV FI_PROVIDER_PATH=${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/opt/mpi/libfabric/lib/prov:/usr/lib/x86_64-linux-gnu/libfabric
ENV PATH=${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/bin:${ONEAPI_ROOT}/mkl/${ONEAPI_VERSION}/bin/:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/bin:${PATH}
