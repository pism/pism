.. include:: ../global.txt

.. _sec-install-common-problems:

Common build problems and solutions
-----------------------------------

We recommend using ``ccmake``, the text-based CMake interface to adjust PISM’s build
parameters. One can also set CMake cache variables using the ``-D`` command-line option
(``cmake -Dvariable=value``) or by editing ``CMakeCache.txt`` in the build directory.

Here are some issues and workarounds we know about.

- The PISM build system uses ``ncgen`` from the NetCDF package to generate the
  configuration file |config-file|. This means that a working NetCDF installation is
  required on both the "host" and the "target" systems when cross-compiling PISM.

- Some systems support static libraries only. To build PISM statically and tell CMake not
  to try to link to shared libraries, set the CMake variable ``Pism_LINK_STATICALLY`` to
  ``ON``.

- You can set ``Pism_LOOK_FOR_LIBRARIES`` to "``OFF``" to disable all heuristics and set
  compiler flags by hand. See `PISM builds <pism-builds_>`_ for examples.

- When linking PISM to *shared* `prerequisite libraries <sec-install-prerequisites>`_ it
  is usually sufficient to link to a library (e.g. NetCDF) and the linker will
  automatically include its dependencies (for NetCDF: HDF5, ``libz``, ``libm``, MPI) using
  the ``NEEDED`` and ``RUNPATH`` headers of its object file. On some systems these headers
  (especially ``RUNPATH``) may not be set, requiring a build system to explicitly list all
  of PISM's dependencies, the dependencies' dependencies, etc, just as when linking to
  *static* libraries.

  To work around this issue, set the CMake variable ``Pism_PKG_CONFIG_STATIC`` to ``YES``.
  This tells PISM to use ``pkg-config``'s ``--static`` flag when looking for PISM's
  dependencies.

  If this proves insufficient and you need to add custom linker flags, set CMake variables
  ``CMAKE_SHARED_LINKER_FLAGS`` (flags used to link shared libraries) and
  ``CMAKE_EXE_LINKER_FLAGS`` (flags used to link executables).
