FROM debian:buster
RUN apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        environment-modules \
        git \
        libgsl-dev \
        libnetcdf-dev \
        libopenmpi-dev \
        netcdf-bin \
        pkgconf \
        procps \
        python

RUN useradd --create-home --system --shell=/bin/false builder && usermod --lock builder
USER builder
WORKDIR /home/builder/

RUN mkdir -p ~/.spack
COPY packages.yaml .spack/

RUN git clone --depth=1 https://github.com/spack/spack.git

RUN . spack/share/spack/setup-env.sh && \
        spack install pism ^petsc~metis~hdf5~hypre~superlu-dist ^fftw~mpi precision=double
