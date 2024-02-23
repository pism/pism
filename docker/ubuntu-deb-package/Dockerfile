FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

RUN <<EOF
    echo "Install build dependencies"
    apt-get update
    apt-get install -y \
    cmake \
    g++ \
    git \
    latexmk \
    libfftw3-dev \
    libgsl-dev \
    libnetcdf-dev \
    libproj-dev \
    libudunits2-dev \
    lintian \
    netcdf-bin \
    petsc-dev \
    python3-netcdf4 \
    python3-sphinx \
    python3-sphinxcontrib.bibtex \
    texlive-latex-extra \
    texlive-latex-recommended \
    ""
    rm -rf /var/lib/apt/lists/*
EOF

ENV PETSC_DIR="/usr/lib/petsc"
