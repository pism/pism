ARG BASE=ghcr.io/pism/pism-base:latest
ARG RUNTIME=ghcr.io/pism/pism-runtime:latest

FROM ${BASE} AS pism-base

COPY docker/pism/pism.sh /tmp/
# Copy the Docker context to /var/tmp/pism-src
COPY . /var/tmp/pism-src
RUN <<EOF
    echo "Install PISM"

    set -eux
    export src_dir=/var/tmp/pism-src
    export build_dir=/var/tmp/build/pism
    export prefix=/opt/pism

    /tmp/pism.sh

    rm -rf ${build_dir}
    rm -rf ${prefix}/include ${prefix}/share/doc/pism
EOF

FROM ${RUNTIME} AS pism-container

# copy PISM's runtime dependencies that were built from sources:
COPY --from=pism-base /opt/hdf5 /opt/hdf5/
COPY --from=pism-base /opt/netcdf /opt/netcdf/
COPY --from=pism-base /opt/petsc /opt/petsc/
COPY --from=pism-base /opt/pism /opt/pism/
COPY --from=pism-base /opt/yac /opt/yac/
COPY --from=pism-base /opt/yaxt /opt/yaxt/

# Add HDF5, NetCDF, PETSc to PATH:
ENV PATH="/opt/hdf5/bin:/opt/netcdf/bin:/opt/petsc/bin:$PATH"

# "Activate" the virtual environment in /opt/pism.
#
# See https://docs.python.org/3/library/venv.html#how-venvs-work for details.
ENV VIRTUAL_ENV=/opt/pism
# Note that this puts the "pism" executable in the path as well
ENV PATH=/opt/pism/bin:$PATH

ENV PYTHONPATH=""
# Add YAC to PYTHONPATH
ENV PYTHONPATH=/opt/yac/local/lib/python3.12/dist-packages/:$PYTHONPATH
# Add PETSc to PYTHONPATH
ENV PYTHONPATH=/opt/petsc/lib/:$PYTHONPATH
# Add PISM to PYTHONPATH
ENV PYTHONPATH=/opt/pism/lib/python3.12/site-packages/:$PYTHONPATH

RUN ln -s /opt/intel/oneapi/setvars.sh /etc/profile.d/02-intel-oneapi-setvars.sh

# Add a non-privileged user
RUN useradd --create-home --system --shell=/bin/false worker && usermod --lock worker
USER worker

ARG ONEAPI_VERSION=2025.3
ARG IMPI_VERSION=2021.17

# Important parts of "source /opt/intel/oneapi/setvars.sh"
ENV ONEAPI_ROOT=/opt/intel/oneapi
ENV LD_LIBRARY_PATH=${ONEAPI_ROOT}/tbb/${IMPI_VERSION}/env/../lib/intel64/gcc4.8:${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/opt/mpi/libfabric/lib:${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/lib:${ONEAPI_ROOT}/mkl/${ONEAPI_VERSION}/lib:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/opt/compiler/lib:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/lib:${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/opt/mpi/libfabric/lib:${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/lib:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/opt/compiler/lib:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/lib
ENV I_MPI_ROOT=${ONEAPI_ROOT}/mpi/${IMPI_VERSION}
ENV FI_PROVIDER_PATH=${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/opt/mpi/libfabric/lib/prov:/usr/lib/x86_64-linux-gnu/libfabric
ENV PATH=${ONEAPI_ROOT}/mpi/${IMPI_VERSION}/bin:${ONEAPI_ROOT}/mkl/${ONEAPI_VERSION}/bin/:${ONEAPI_ROOT}/compiler/${ONEAPI_VERSION}/bin:${PATH}
