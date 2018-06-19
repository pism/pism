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
       usefully-filtered information on options, for example ``pisms -help | grep cold``.

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

PISM file I/O performance
^^^^^^^^^^^^^^^^^^^^^^^^^

When working with fine grids\ [#]_, the time PISM spends writing output files,
spatially-varying diagnostic files, or backup files can become significant.

It turns out that it is a lot faster to read and write files using the ``t,y,x,z`` storage
order, as opposed to the more convenient (e.g. for NetCDF tools) ``t,z,y,x`` order. The
reason is that PISM uses the ``y,x,z`` order internally,\ [#]_ and therefore writing an
array in a different order is an inherently-expensive operation.

You can, however, choose any one of the three supported output orders using the
:opt:`-o_order` option with one of ``xyz``, ``yxz``, and ``zyx`` as the argument.

To transpose dimensions in an existing file, use the ``ncpdq`` ("permute dimensions
quickly") tool from the NCO_ suite. For example, run

.. code-block:: none

   ncpdq -a t,z,zb,y,x bad.nc good.nc


to turn ``bad.nc`` (with any inconvenient storage order) into ``good.nc`` using the
``t,z,y,x`` order.

PISM also supports NetCDF-4 parallel I/O, which gives better performance in
high-resolution runs and avoids NetCDF-3 file format limitations. (In a NetCDF-3 file a
variable record cannot exceed 4 gigabytes.) Build PISM with parallel NetCDF-4 and use
:opt:`-o_format` ``netcdf4_parallel`` to enable this code.

In addition to ``-o_format netcdf4_parallel`` and ``netcdf3`` (default) modes, PISM can be
built with PnetCDF for best I/O performance. The option ``-o_format pnetcdf`` turns "on"
PnetCDF I/O code. (PnetCDF seems to be somewhat fragile, though, so use at your own risk.)

.. rubric:: Footnotes

.. [#] For example, resolutions of 2km and higher on the whole-Greenland scale.
.. [#] This is not likely to change.
