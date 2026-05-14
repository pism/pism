# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build, install, run

```bash
cmake -B build [options]            # configure
cmake --build build -j              # compile
cmake --install build --prefix DIR  # install (optional)
```

The CLI binary is `pism` (renamed from `pismr` in 2.1; `pismv` is gone — use `pism -test X`).

Frequently used CMake options (off by default unless noted):

- `-DPism_BUILD_PYTHON_BINDINGS=ON` — build SWIG-generated `PISM.cpp` Python module. **Requires `petsc4py` importable by the host Python at configure time** (`CMake/FindPETSc4Py.cmake` runs `python -c "import petsc4py"`).
- `-DPism_DEBUG=ON` — enables extra runtime checks in the C++ code.
- `-DPism_USE_PROJ=ON`, `-DPism_USE_YAC=ON`, `-DPism_USE_PARALLEL_NETCDF4=ON`, `-DPism_USE_PNETCDF=ON` — opt into the corresponding optional dependency.
- `-DCMAKE_BUILD_TYPE=Release|Debug|RelWithDebInfo` — defaults to `Release` if unset.

`Pism_DEBIAN_SYSTEMWIDE=ON` switches off RPATH and changes the Python install scheme — do not enable it for normal developer builds.

## Tests

PISM uses CTest. From the build directory:

```bash
ctest -j N                                 # run everything
ctest -R "<regex>"                         # run tests whose names match
ctest --rerun-failed --output-on-failure   # also exposed as: cmake --build build --target retest
ctest -R "Python:"                         # also exposed as: cmake --build build --target test-python
```

Test-name conventions used by `test/CMakeLists.txt`:

- `Config:*` — config-file metadata and ordering checks (always run; required fixture for most others).
- `Python:unittest:<name>` — Python regression tests; only registered when `Pism_BUILD_PYTHON_BINDINGS=ON`. They run `python -m unittest -v <file>` with `PYTHONPATH` pointed at `${PROJECT_BINARY_DIR}/site-packages` and `PETSC_OPTIONS=-config <build>/pism_config.nc`.
- Verification suites and regression scripts live under `test/` and `test/regression/`.

To run a single Python regression test against an in-tree build without going through CTest:

```bash
PYTHONPATH=$(pwd)/build/site-packages PETSC_OPTIONS="-config $(pwd)/build/pism_config.nc" \
    python -m unittest -v test/<name>.py
```

## Architecture (big picture)

- **Single C++ library + single CLI.** `src/CMakeLists.txt` builds `libpism` (most of `src/`) and the `pism` executable (`src/pism.cc`). Optional auxiliary executables: `pism_btutest`; the `pism_test_*` SSAFE/SSAFD test harnesses under `src/stressbalance/ssa/tests/`. `IceModel` (in `src/icemodel/`) is the time-stepping orchestrator.
- **Physics is partitioned by subdirectory under `src/`** — each is a self-contained component with its own CMake target that gets linked into `libpism`: `stressbalance/` (SIA, SSA, Blatter), `energy/`, `hydrology/`, `coupler/{atmosphere,ocean,surface,frontalmelt}`, `geometry/`, `earth/` (bed deformation), `frontretreat/` (calving + front retreat), `age/`, `tracer/`, `inverse/` (SSA inversion / `pismi.py`), `regional/`, `rheology/`, `verification/`.
- **PETSc is the parallel backbone.** Most state lives on `DMDA`-managed distributed arrays; the higher-level wrapper is `pism::array::Array` (`src/util/array/`). Use these instead of raw PETSc objects when adding code.
- **`src/pism_config.cdl` is the single source of truth for all configuration parameters** (~2900 lines, alphabetically sorted by parameter name). Each parameter declares `_doc`, `_type`, optionally `_units` and `_option`. **Two CI checks enforce this**: `Config:metadata_structure` (every param has the required attributes) and `Config:parameters_are_sorted` (alphabetical order). When adding or renaming a parameter, edit `pism_config.cdl` and keep it sorted; CMake compiles it into `pism_config.nc` which ships with the install. Read parameters via `Config::get_*()`.
- **Python layer has two parts.** SWIG generates `_cpp.so` + `cpp.py` from 28 `.i` files in `src/pythonbindings/` — the master is `PISM.i`, the C++ glue is `pism_python.cc`, and the resulting module is exposed as `PISM.cpp`. Pure-Python helpers (`util.py`, `vec.py`, `sia.py`, `ssa.py`, `model.py`, `logging.py`, `invert/`) live in `site-packages/PISM/`; `site-packages/siple/` is a co-shipped optimization library used by the SSA inversion code. Both trees get installed into the active interpreter's `site-packages`.
- **`util/pism_*` Python scripts** in `util/` are user-facing CLI helpers (timeline manipulation, async writer, NetCDF utilities); they ship in `bin/` after install.

## Change-log discipline

`CHANGES.rst` is part of the contribution checklist (see `CONTRIBUTING.rst`). For any user-visible behavior change, parameter rename, new feature, or bug fix, add a bullet under the topmost `Changes since vX.Y.Z` heading. Keep entries short, present-tense, and use single-backtick `literal` markup (the file sets `.. default-role:: literal` so backticks render as code).

## External docs

Authoritative installation, usage, and contributing docs live at https://www.pism.io/docs/ and the source is under `doc/sphinx/`. The README at the repo root mostly points to those.
