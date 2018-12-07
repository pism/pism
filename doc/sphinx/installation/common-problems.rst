.. include:: ../global.txt

.. _sec-install-common-problems:

Common build problems and solutions
-----------------------------------

We recommend using ``ccmake``, the text-based CMake interface to adjust PISM’s build
parameters. One can also set CMake cache variables using the ``-D`` command-line option
(``cmake -Dvariable=value``) or by editing ``CMakeCache.txt`` in the build directory.

Here are some issues we know about.

- If you are compiling PISM on a system using a cross-compiler, you will need to disable
  CMake’s tests trying to determine if PETSc is installed properly. To do this, set
  ``PETSC_EXECUTABLE_RUNS`` to "yes".

  To tell CMake where to look for libraries for the target system, see `CMake cross
  compiling <CMake-cross-compiling_>`_ and the paragraph about ``CMAKE_FIND_ROOT_PATH``
  in particular.

- Note that the PISM build system uses ``ncgen`` from the NetCDF package to generate the
  configuration file |config-file|. This means that a working NetCDF installation is
  required on both the "host" and the "target" systems when cross-compiling PISM.

- Some systems support static libraries only. To build PISM statically and tell CMake not
  to try to link to shared libraries, set ``Pism_LINK_STATICALLY`` to ``ON`` using
  ``ccmake``.

- You can set ``Pism_LOOK_FOR_LIBRARIES`` to "``OFF``" to disable all heuristics and set
  compiler flags by hand. See `PISM builds <pism-builds_>`_ for examples.
