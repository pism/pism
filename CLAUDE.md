# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PISM (Parallel Ice Sheet Model) is an open-source, parallel, high-resolution ice sheet model written primarily in C++ (C++11). It uses MPI for parallelism, PETSc for numerical solvers, and NetCDF for file I/O. Developed jointly at UAF and PIK.

## Build Commands

PISM uses CMake with out-of-source builds:

```bash
# Configure (from build directory)
cmake /path/to/pism/source \
  -DCMAKE_INSTALL_PREFIX=/install/path \
  -DPetsc_DIR=/path/to/petsc

# Build
make -j$(nproc)

# Install
make install
```

Key CMake options:
- `Pism_BUILD_PYTHON_BINDINGS=ON` - Build SWIG Python bindings
- `Pism_DEBUG=ON` - Enable runtime checks
- `Pism_USE_PROJ=ON` - PROJ for coordinate transforms
- `Pism_USE_PARALLEL_NETCDF4=ON` - Parallel NetCDF-4 I/O
- `Pism_PEDANTIC_WARNINGS=ON` (default) - Pedantic compiler warnings

Default build type is Release. The build generates `pism_config.nc` from `src/pism_config.cdl` using `ncgen`.

## Running Tests

```bash
# Run all tests (from build directory)
make test
# or
ctest --output-on-failure

# Run tests in parallel
CTEST_PARALLEL_LEVEL=4 ctest

# Re-run only failed tests
make retest

# Run only Python tests (requires nose and Python bindings)
make test-python
```

Tests are organized in `test/` and `test/regression/`. Shell-script regression tests (test_XX.sh) take the build directory and mpiexec path as arguments. Python tests use `nose` and require `Pism_BUILD_PYTHON_BINDINGS=ON`.

## Required Dependencies

- PETSc >= 3.7.0
- NetCDF >= 4.4
- GSL >= 1.15
- FFTW3 >= 3.1
- UDUNITS2
- MPI (C and C++)
- ncgen (for building pism_config.nc)

## Architecture

### Main Executable

`src/pism.cc` - Entry point supporting multiple modes:
- Standard ice sheet evolution (default)
- EISMINT II experiments (`-eisII`)
- Verification tests (`-test`)
- Regional modeling (`-regional`)

### Core Library (`libpism`)

All source lives under `src/` and compiles into a single shared library. Major subsystems:

**Ice model orchestration** (`src/icemodel/`): `IceModel` is the central class that owns all components and manages the time-stepping loop. Implementation split across multiple files by concern (energy, diagnostics, output, timestepping, etc.).

**Stress balance** (`src/stressbalance/`): Hierarchy of stress balance solvers:
- `sia/` - Shallow Ice Approximation
- `ssa/` - Shallow Shelf Approximation (FD and FEM variants)
- `blatter/` - Blatter-Pattyn full 3D solver

**Couplers** (`src/coupler/`): Modular atmosphere/ocean/surface forcing via composable model chains:
- `atmosphere/` - Atmospheric forcing
- `surface/` - Surface mass balance
- `ocean/` - Ocean interaction (includes PICO/PICOP)
- `frontalmelt/` - Frontal melt parameterizations

**Energy** (`src/energy/`): Enthalpy-based energy conservation with bedrock thermal unit (`BedThermalUnit`) having full and minimal implementations.

**Geometry** (`src/geometry/`): Ice sheet geometry evolution, grounding line handling, flux limiters, UNO advection scheme.

**Utilities** (`src/util/`): Core infrastructure:
- `array/` - Distributed array types (Scalar, Vector)
- `io/` - NetCDF I/O layer
- `petscwrappers/` - C++ RAII wrappers for PETSc objects
- `fem/` - Finite element infrastructure
- `Component.hh` - Base class for model components
- `Context.hh` - Shared runtime context (grid, config, MPI communicator)

**Other subsystems:**
- `hydrology/` - Basal hydrology models
- `basalstrength/` - Yield stress models (Mohr-Coulomb, constant, optimal)
- `earth/` - Glacial isostatic adjustment (Lingle-Clark)
- `frontretreat/` - Calving laws and front retreat
- `rheology/` - Ice flow laws
- `age/` - Age tracking and isochrones
- `inverse/` - Inverse modeling
- `fracturedensity/` - Fracture density evolution
- `regional/` - Regional/outlet glacier mode
- `verification/` - Exact solutions for code verification

### Component Pattern

Model components inherit from `Component` (in `src/util/Component.hh`) which provides access to the shared `Context` (grid, config, MPI communicator, unit system). Components declare diagnostics and are initialized/updated through a uniform interface.

### Configuration

Parameters are defined in `src/pism_config.cdl` (NetCDF CDL format) and compiled to `pism_config.nc`. Runtime access is through the `Config` class. During testing, `.petscrc` in the build directory points to the local `pism_config.nc`.

### Header Include Convention

Headers are included as `#include "pism/subsystem/Header.hh"`. The build system creates symlinks in `${BUILD_DIR}/include/pism/` pointing to source directories, allowing the same include paths for both in-tree builds and installed headers.

### External Libraries

`src/external/` contains vendored code: `calcalcs` (calendar calculations), `cubature` (numerical integration), `nlohmann` (JSON).
