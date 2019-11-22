.. include:: ../global.txt

.. _sec-install-parallel-io-libs:

Installing parallel I/O libraries
=================================

High-resolution simulations (e.g. modeling the entire Greenland ice sheet using the 900 m
resolution and a grid of 1700 by 3000 points) can produce huge amounts of data and I/O can
become a bottleneck.

PISM supports several parallel I/O approaches that take advantage of parallel NetCDF_,
PnetCDF_, and NCAR ParallelIO_. The administrators of your HPC system should be able to
help you install these libraries, but it may be easier to install them in the "home"
directory instead.

This section describes the steps needed to build

- NetCDF with parallel I/O based on HDF5 (needed to use PISM's option :opt:`-o_format
  netcdf4_parallel`),
- PNetCDF (needed to use PISM's option :opt:`-o_format pnetfdf`),
- ParallelIO (needed to use options :opt:`-o_format pio_netcdf4p`,
  :opt:`-o_format pio_netcdf4c`, :opt:`-o_format pio_pnetcdf`, :opt:`-o_format pio_netcdf`).

Scripts below install libraries in ``~/local/library_name``, using
``~/local/build/library_name`` to build them.

Section :ref:`sec-install-local-libraries` explains how build PISM with these libraries.

Section :ref:`sec-pism-io-performance` explains how to use them in PISM.

.. _sec-install-parallel-netcdf:

Installing HDF5-based parallel NetCDF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _sec-install-hdf5:

Installing  HDF5
~~~~~~~~~~~~~~~~

.. literalinclude:: code/hdf5.sh
   :language: bash
   :linenos:
   :lines: 7-

To compile parallel HDF5 one should use the MPI compiler wrapper ``mpicc`` and run
``configure`` with the option ``--enable-parallel``. The flag ``-w`` is not important: it
hides numerous compiler warnings emitted when building HDF5.

.. _sec-install-netcdf:

Installing NetCDF
~~~~~~~~~~~~~~~~~

.. literalinclude:: code/netcdf.sh
   :language: bash
   :linenos:
   :lines: 7-

Here we use the same compiler wrapper and set ``CPPFLAGS`` and ``LDFLAGS`` to select the
parallel HDF5 library installed earlier. The option ``--enable-netcdf4`` is required for
parallel I/O; ``--disable-dap`` is not required (it disables a NetCDF feature not used by
PISM).

.. _sec-install-pnetcdf:

Installing PnetCDF
^^^^^^^^^^^^^^^^^^

.. literalinclude:: code/pnetcdf.sh
   :language: bash
   :linenos:
   :lines: 7-

Here we disable PnetCDF's C++ and Fortran APIs and build the shared library.

.. _sec-install-parallelio:

Installing NCAR ParallelIO
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: code/parallelio.sh
   :language: bash
   :linenos:
   :lines: 7-

Here we use CMake's variable ``CMAKE_FIND_ROOT_PATH`` to tell CMake to use libraries in
``~/local/netcdf`` and ``~/local/pnetcdf``, to install in ``~/local/parallelio``, and to
disable ParallelIO features not used by PISM.
