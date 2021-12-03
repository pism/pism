.. include:: ../../global.txt

.. _sec-input-output:

Input and output
----------------

PISM is a program that reads NetCDF files and then outputs NetCDF files. Table
:numref:`tab-input-output-options` summarizes command-line options controlling the most
basic ways to input and output NetCDF files when starting and ending PISM runs.

.. list-table:: Basic NetCDF input and output options
   :name: tab-input-output-options
   :header-rows: 1
   :widths: 1,2

   * - Option
     - Description

   * - :opt:`-i`
     - Chooses a PISM output file (NetCDF format) to initialize or restart from. See
       section :ref:`sec-initboot`.

   * - :opt:`-bootstrap`
     - Bootstrap from the file set using :opt:`-i` using heuristics to "fill in" missing
       fields. See section :ref:`sec-initboot`.

   * - :opt:`-ssa_read_initial_guess false`
     - Turns off reading the ``ubar_ssa`` and ``vbar_ssa`` velocities saved by a previous
       run using the ``ssa`` or ``ssa+sia`` stress balance (see section
       :ref:`sec-stressbalance`).

   * - :opt:`-o`
     - Chooses the output file name.  Default name is ``unnamed.nc``.

   * - :opt:`-o_size` ``size_keyword``
     - Chooses the size of the output file to produce. Possible sizes are

       - ``none`` (*no* output file at all),
       - ``small`` (only variables necessary to restart PISM),
       - ``medium`` (the default, includes diagnostic quantities listed in the
         configuration parameter :config:`output.sizes.medium`, if they are available in
         the current PISM setup),
       - ``big_2d`` (same as ``medium``, plus variables listed in
         :config:`output.sizes.big_2d`), and
       - ``big`` (same as ``big_2d``, plus variables listed in
         :config:`output.sizes.big`).

:numref:`tab-stdout` lists the controls on what is printed to the standard output.
Note the ``-help`` and ``-usage`` options for getting help at the command line.

.. list-table:: Options controlling PISM's standard output
   :header-rows: 1
   :name: tab-stdout
   :widths: 1,2

   * - Option
     - Description

   * - :opt:`-help`
     - Brief descriptions of the many PISM and PETSc options. The run occurs as usual
       according to the other options. (The option documentation does not get listed if
       the run didn't get started properly.) Use with a pipe into ``grep`` to get
       usefully-filtered information on options, for example ``pismr -help | grep cold``.

   * - :opt:`-info`
     - Gives information about PETSc operations during the run.

   * - :opt:`-list_diagnostics`
     - Prints a list of all available diagnostic outputs (time series and spatial) for the
       run with the given options. Stops run after printing the list.

   * - :opt:`-log_summary`
     - At the end of the run gives a performance summary and also a synopsis of the PETSc
       configuration in use.

   * - :opt:`-options_left`
     - At the end of the run shows an options table which will indicate if a user option
       was not read or was misspelled.

   * - :opt:`-usage`
     - Short summary of PISM executable usage, without listing all the options, and
       without doing the run.

   * - :opt:`-verbose`
     - Increased verbosity of standard output. Usually given with an integer level;
       0,1,2,3,4,5 are allowed. If given without argument then sets level 3, while
       ``-verbose 2`` is the default (i.e. equivalent to no option). At the extremes,
       ``-verbose 0`` produces no stdout at all, ``-verbose 1`` prints only warnings and a
       few high priority messages, and ``-verbose 5`` spews a lot of usually-undesirable
       stuff. ``-verbose 3`` output regarding initialization may be useful.

   * - :opt:`-version`
     - Show version numbers of PETSc and PISM.

The following sections describe more input and output options, especially related to
saving quantities during a run, or adding to the "diagnostic" outputs of PISM.

.. _sec-pism-io-performance:

PISM's I/O performance
^^^^^^^^^^^^^^^^^^^^^^

When working with fine grids (resolutions of 2km and higher on the whole-Greenland scale,
for example), the time PISM spends writing output files, spatially-varying diagnostic
files, or backup files can become significant.

For fast file I/O the order of dimensions of a NetCDF variable in an output file has to
match the order used by PISM in memory, so we use the ``time,y,x,z`` storage order instead of
the more convenient (e.g. for NetCDF tools) order ``time,z,y,x``.

To transpose dimensions in an existing file, use the ``ncpdq`` ("permute dimensions
quickly") tool from the NCO_ suite. For example, run

.. code-block:: none

   ncpdq -a time,z,zb,y,x bad.nc good.nc

to turn ``bad.nc`` (with any inconvenient storage order) into ``good.nc`` using the
``time,z,y,x`` order.

PISM also supports parallel I/O using parallel NetCDF_, PnetCDF_, or ParallelIO_, which
can give better performance in high-resolution runs.

Use the command-line option :opt:`-o_format` (parameter :config:`output.format`) to choose
the approach to use when writing to output files (see :numref:`tab-output-format`). The
``netcdf4_parallel`` requires parallel NetCDF, ``pnetcdf`` requires PnetCDF, and
``pio_...`` require ParallelIO build with parallel NetCDF and PnetCDF. Section
:ref:`sec-install-pism-cmake-options`) explains how to select these libraries when
building PISM.

.. note::

   When built with parallel NetCDF or PnetCDF (or both) PISM attempts to choose the best
   way to *read* from input files and this logic appears to work well. This is why there
   is no ``-i_format``.

.. csv-table:: Methods of writing to output files
   :name: tab-output-format
   :header: ``-o_format`` argument, Description

   ``netcdf3``, (default); serialized I/O from rank 0 (NetCDF-3 file)
   ``netcdf4_parallel``, parallel I/O using NetCDF (HDF5-based NetCDF-4 file)
   ``pnetcdf``, parallel I/O using PnetCDF (CDF5 file)
   ``pio_pnetcdf``,  parallel I/O using ParallelIO (CDF5 file)
   ``pio_netcdf4p``, parallel I/O using ParallelIO (HDF5-based NetCDF-4 file)
   ``pio_netcdf4c``, serial I/O using ParallelIO (*compressed* HDF5-based NetCDF-4 file)
   ``pio_netcdf``,   serial I/O using ParallelIO (using data aggregation in ParallelIO)

The ParallelIO library can aggregate data in a subset of processes used by PISM. To choose
a subset, set

- :config:`output.pio.n_writers` number of "writers"
- :config:`output.pio.base` the index of the first writer
- :config:`output.pio.stride` interval between writers

.. note::

   The CDF5 file format is a large-variable extension of the NetCDF-3 file format
   developed by the authors of PnetCDF. This format is supported by NetCDF since version
   4.4.

We recommend performing a number of test runs to determine the best choice for your
simulations.

In our test runs on 120 cores (whole Greenland setup on a 900m grid) ``pio_pnetcdf`` with
:config:`output.pio.n_writers` set to the number of cores used by PISM (120) gave the best
performance.

.. note::

   It is important to make sure that PISM's output files are written to a parallel file
   system and this file system is configured to achieve optimal performance.

   On Lustre_ (a common parallel file systems) the theoretical throughput when writing to
   a file depends on the number of *object storage targets* used to store it: if a target
   can write 500 MiB/s, a file spread over 2 could be written at 1000 MiB/s assuming that
   we are writing to both of them at the same time, and so on.

   **For maximum speed we want to distribute an output file over all available targets.**

   To do this:

   1. Create a directory that will contain PISM output files (``output_directory`` below).
   2. Run

      .. code-block:: bash

         lfs setstripe -c -1 output_directory

      This sets the "stripe count" to ``-1``, which means "all".

      Now all files in ``output_directory`` and all its sub-directories can use all
      available targets.
