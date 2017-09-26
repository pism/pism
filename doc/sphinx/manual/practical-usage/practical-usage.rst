.. include:: ../../global.rst

.. |diagnostics| replace:: **FIXME: link to the list of diagnostics**

.. contents::

.. _sec-practical-usage:

Practical usage
===============

.. _sec-input-output:

Input and output
----------------

PISM is a program that reads NetCDF files and then outputs NetCDF files. Table
:numref:`tab-input-output-options` summarizes command-line options controlling the most
basic ways to input and output NetCDF files when starting and ending PISM runs.

.. list-table:: Basic NetCDF input and output options
   :name: tab-input-output-options
   :header-rows: 1

   * - Option
     - Description

   * - :opt:`-i`
     - Chooses a PISM output file (NetCDF format) to initialize or restart from. See
       section :ref:`sec-initboot`.

   * - :opt:`-bootstrap`
     - Bootstrap from the file set using :opt:`-i` using heuristics to "fill in" missing
       fields. See section :ref:`sec-initboot`.

   * - :opt:`-dontreadSSAvels`
     - Turns off reading the ``ubar_ssa`` and ``vbar_ssa`` velocities saved by a previous
       run using the ``ssa`` or ``ssa+sia`` stress balance (see section
       :ref:`sec-stressbalance`).

   * - :opt:`-o`
     - Chooses the output file name.  Default name is ``unnamed.nc``.

   * - :opt:`-o_size` ``size_keyword``
     - Chooses the size of the output file to produce. Possible sizes are ``none`` (*no*
       output file at all), ``small`` (only variables necessary to restart PISM),
       ``medium`` (the default, includes a few diagnostic quantities), ``big`` (writes all
       the variables mentioned in |diagnostics|), and ``big_2d`` writes all 2D variables
       but only 3D variables that are model state. Configuration variables
       :config:`output.sizes.medium`, :config:`output.sizes.big`, and
       :config:`output.sizes.big_2d` list the written variables for those sizes.

:numref:`tab-stdout` lists the controls on what is printed to the standard output.
Note the ``-help`` and ``-usage`` options for getting help at the command line.

.. list-table:: Options controlling PISM's standard output
   :header-rows: 1
   :name: tab-stdout

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

When working with fine grids [#]_, the time PISM spends writing output files,
spatially-varying diagnostic files, or backup files can become significant.

It turns out that it is a lot faster to read and write files using the ``t,y,x,z`` storage
order, as opposed to the more convenient (e.g. for NetCDF tools) ``t,z,y,x`` order. The
reason is that PISM uses the ``y,x,z`` order internally, [#]_ and therefore writing an
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


.. _sec-saving-time-series:

Saving time series of scalar diagnostic quantities
--------------------------------------------------

It is also possible to save time-series of certain scalar diagnostic quantities using a
combination of the options ``-ts_file``, ``-ts_times``, and ``-ts_vars``. For example,

.. code-block:: none

   pismr -i foo.nc -y 1e4 -o output.nc -ts_file time-series.nc \
         -ts_times 0:1:1e4 -ts_vars volume_glacierized,area_glacierized_grounded


will run for 10000 years, saving total ice volume and grounded ice area to
``time-series.nc`` *yearly*. See tables :numref:`tab-time-series-opts` for the list of
options and tables |diagnostics| for the full list of supported time-series.

Note that, similarly to the snapshot-saving code (section :ref:`sec-snapshots`), this
mechanism does not affect adaptive time-stepping. Here, however, PISM will save exactly
the number of time-series records requested, *linearly interpolated onto requested times*.

Omitting the ``-ts_vars`` option makes PISM save *all* available variables, as listed in
tables |diagnostics|. Because scalar time-series take minimal storage space, compared to
spatially-varying data, this is usually a reasonable choice. Run PISM with the
:opt:`-list_diagnostics` option to see the list of all available time-series.

If the file ``foo.nc``, specified by ``-ts_file foo.nc``, already exists then by default
the existing file will be moved to ``foo.nc~`` and the new time series will go into
``foo.nc``. To append the time series onto the end of the existing file, use option
``-ts_append``.

PISM buffers time-series data and writes it at the end of the run, once 10000 values are
stored, or when an ``-extra_file`` is saved, whichever comes first. Sending an ``USR1``
(or ``USR2``) signal to a PISM process flushes these buffers, making it possible to
monitor the run. (See section :ref:`sec-signal` for more about PISM's signal handling.)

.. list-table:: Command-line options controlling saving scalar time-series
   :name: tab-time-series-opts
   :header-rows: 1

   * - Option
     - Description

   * - :opt:`-ts_file`
     - Specifies the file to save to.

   * - :opt:`-ts_times`
     - Specifies times to save at as a MATLAB-style range :math:`a:\Delta t:b`, a
       comma-separated list, or a keyword (``hourly``, ``daily``, ``monthly``,
       ``yearly``). See section :ref:`sec-saving-diagnostics`.

   * - :opt:`-ts_vars`
     - Comma-separated list of variables. Omitting this option is equivalent to listing
       the *all* variables.

   * - :opt:`-ts_append`
     - Append time series to file if it already exists. No effect if file does not yet
       exist.

Besides the above information on usage, here are comments on the physical significance of
several scalar diagnostics:

- For each variable named ``..._flux``, positive values mean ice sheet mass gain.

- PISM reports ice volume, ice mass, and several other quantities for "glacierized" areas.
  These quantities do not include contributions from areas where the ice thickness is
  equal to or below the value of the configuration parameter
  ``output.ice_free_thickness_standard`` (in meters). Corresponding "nonglacierized"
  quantities *do* include areas with a thin, "seasonal" ice cover.

- The ``sub_shelf_ice_flux`` may be non-zero even if ``area_glacierized_shelf`` (floating
  ice area) is zero. This is due to the fact that during time-stepping fluxes are computed
  before calving is applied, and the ice area is computed *after* calving. Hence ice that
  is calved off experiences top-surface and basal fluxes, but does not contribute to the
  reported area. This is a small error that approaches zero as the grid is refined. In
  this case ``sub_shelf_ice_flux`` should be added to the calving flux during
  post-processing. [#]_

- Ice volume and area are computed and then split among floating and grounded portions:
  ``volume_glacierized`` :math:`\mapsto` (``volume_glacierized_shelf``,
  ``volume_glacierized_grounded``) while ``area_glacierized`` :math:`\mapsto`
  (``area_glacierized_shelf``, ``area_glacierized_grounded``). The volumes have units
  :math:`m^3` and the areas have units :math:`m^2`.

- The thermodynamic state of the ice sheet can be assessed, in part, by the amount of cold
  or temperate ("``temp``") ice. Thus there is another splitting: ``volume_glacierized``
  :math:`\mapsto` (``volume_glacierized_cold``, ``volume_glacierized_temperate``) and
  ``area_glacierized`` :math:`\mapsto`
  (``area_glacierized_cold_base``, ``area_glacierized_temperate_base``).

- If a PISM input file contains the ``proj4`` global attribute with a PROJ.4 string
  defining the projection then PISM computes corrected cell areas using this information,
  grid parameters, and the WGS84 reference ellipsoid. This yields areas and volumes with
  greater accuracy.

- The sea-level-relevant ice volume ``slvol`` is the total grounded ice volume minus the
  amount of ice, that, in liquid form, would fill up the regions with bedrock below sea
  level, if this ice were removed. That is, ``slvol`` is the sea level rise potential of
  the ice sheet at that time. The result is reported in sea-level equivalent, i.e. meters
  of sea level rise.

- Fields ``max_diffusivity`` and ``max_hor_vel`` relate to PISM time-stepping. These
  quantities appear in per-time-step form in the standard output from PISM (i.e. at
  default verbosity). ``max_diffusivity`` determines the length of the mass continuity
  sub-steps for the SIA stress balance (sub-)model. ``max_hor_vel`` determines the
  CFL-type restriction for mass continuity and conservation of energy contributions of the
  SSA stress balance (i.e. sliding) velocity.

.. note:: Document "Scalar time-series supported by PISM"

          (with or without the hydrology model)

.. _sec-saving-diagnostics:

Saving time series of spatially-varying diagnostic quantities
-------------------------------------------------------------

Sometimes it is useful to have PISM save a handful of diagnostic *maps* at some interval
like every 10 years or even every month. One can use snapshots (section
:ref:`sec-snapshots`), but doing so can easily fill your hard-drive because snapshots are
complete (i.e. re-startable) model states. Sometimes you want a *subset* of model
variables saved frequently in an output file.

Use options ``-extra_file``, ``-extra_times``, and ``-extra_vars`` for this. For example,

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc -extra_file extras.nc \
         -extra_times 0:10:1e4 -extra_vars velsurf_mag,velbase_mag


will run for 10000 years, saving the magnitude of horizontal velocities at the ice surface
and at the base of ice every 10 years. Times are specified using a comma-separated list or
a MATLAB-style range. See :numref:`tab-extras` for all the options controlling this
feature. Tables |diagnostics| list all the variable choices.

Note that options :opt:`-extra_times`, :opt:`-save_times`, :opt:`-ts_times` take *dates*
if a non-trivial calendar is selected. For example,

.. code-block:: bash

   pismr ... -extra_times 10       # every 10 years
   pismr ... -extra_times 2days    # every 2 days
   pismr ... -calendar gregorian -extra_times 1-1-1:daily:11-1-1 # daily for 10 years
   pismr ... -calendar gregorian -extra_times daily -ys 1-1-1 -ye 11-1-1
   pismr ... -calendar gregorian -extra_times 2hours -ys 1-1-1 -ye 1-2-1


The step in the range specification can have the form ``Nunit``, for example ``5days``.
Units based on "months" and "years" are not supported if a non-trivial calendar is
selected.

In addition to specifying a constant step in ``-extra_times a:step:b`` one can save every
hour, day, month, or every year by using ``hourly``, ``daily``, ``monthly`` or ``yearly``
instead of a number; for example

.. code-block:: none

   pismr -i foo.nc -y 100 -o output.nc -extra_file extras.nc \
         -extra_times 0:monthly:100 -extra_vars dHdt

will save the rate of change of the ice thickness every month for 100 years. With
``-calendar none`` (the default), "monthly" means "every :math:`\frac 1 {12}` of the
year", and "yearly" is "every :math:`3.14\dots\times10^7`" seconds, otherwise PISM uses
month lengths computed using the selected calendar.

It is frequently desirable to save diagnostic quantities at regular intervals for the
whole duration of the run; options :opt:`-extra_times`, :opt:`-ts_times`, and
:opt:`-save_times` provide a shortcut. For example, use ``-extra_times yearly`` to save at
the end of every year.

This is especially useful when using a climate forcing file to set run duration:

.. code-block:: none

   pismr -i foo.nc -surface given -surface_given_file climate.nc \
         -calendar gregorian -time_file climate.nc \
         -extra_times monthly -extra_file ex.nc -extra_vars thk


will save ice thickness at the end of every month while running PISM for the duration of
climate forcing data in ``climate.nc``.

Times given using ``-extra_times`` describe the reporting intervals by giving the
endpoints of these reporting intervals. The save itself occurs at the end of each
interval. This implies, for example, that ``0:1:10`` will produce 10 records at times
1,...,10 and *not* 11 records.

If the file ``foo.nc``, specified by ``-extra_file foo.nc``, already exists then by
default the existing file will be moved to ``foo.nc~`` and the new time series will go
into ``foo.nc``. To append the time series onto the end of the existing file, use option
``-extra_append``.

The list of available diagnostic quantities depends on the model setup. For example, a run
with only one vertical grid level in the bedrock thermal layer will not be able to save
``litho_temp``, an SIA-only run does not use a basal yield stress model and so will not
provide ``tauc``, etc. To see which quantities are available in a particular setup, use
the :opt:`-list_diagnostics` option, which prints the list of diagnostics and stops.

The ``-extra_file`` mechanism modifies PISM's adaptive time-stepping scheme so as to step
to, and save at, *exactly* the times requested. By contrast, as noted in subsection
:ref:`sec-saving-time-series`, the ``-ts_file`` mechanism does not alter PISM's time-steps
and instead uses linear interpolation to save at the requested times in between PISM's
actual time-steps.

.. list-table:: Command-line options controlling extra diagnostic output
   :name: tab-extras
   :header-rows: 1

   * - Option
     - Description

   * - :opt:`-extra_file`
     - Specifies the file to save to; should be different from the output (:opt:`-o`)
       file.

   * - :opt:`-extra_times`
     - Specifies times to save at either as a MATLAB-style range :math:`a:\Delta t:b` or a
       comma-separated list.

   * - :opt:`-extra_vars`
     - Comma-separated list of variables

   * - :opt:`-extra_split`
     - Save to separate files, similar to :opt:`-save_split`.

   * - :opt:`-extra_append`
     - Append variables to file if it already exists. No effect if file does not yet
       exist, and no effect if :opt:`-extra_split` is set.

.. note:: Document "Scalar 3D diagnostic quantities"

.. note:: Document "Vector 3D diagnostic quantities"

.. note:: Document "Scalar 2D diagnostic quantities"

.. note:: Document "Vector 2D diagnostic quantities"

.. _sec-snapshots:

Saving re-startable snapshots of the model state
------------------------------------------------

Sometimes you want to check the model state every 1000 years, for example. One possible
solution is to run PISM for a thousand years, have it save all the fields at the end of
the run, then restart and run for another thousand, and etc. This forces the adaptive
time-stepping mechanism to stop *exactly* at multiples of 1000 years, which may be
desirable in some cases.

If saving exactly at specified times is not critical, then use the ``-save_file`` and
``-save_times`` options. For example,

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc -save_file snapshots.nc \
         -save_times 1000:1000:10000

starts a PISM evolution run, initializing from ``foo.nc``, running for 10000 years and
saving snapshots to ``snapshots.nc`` at the first time-step after each of the years 1000,
2000, ..., 10000.

We use a MATLAB-style range specification, :math:`a:\Delta t:b`, where :math:`a,\Delta
t,b` are in years. The time-stepping scheme is not affected, but as a consequence we do
not guarantee producing the exact number of snapshots requested if the requested save
times have spacing comparable to the model time-steps. This is not a problem in the
typical case in which snapshot spacing is much greater than the length of a typical time
step.

It is also possible to save snapshots at intervals that are not equally-spaced by giving
the ``-save_times`` option a comma-separated list. For example,

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc -save_file snapshots.nc \
         -save_times 1000,1500,2000,5000

will save snapshots on the first time-step after years 1000, 1500, 2000 and 5000. The
comma-separated list given to the ``-save_times`` option can be at most 200 numbers long.

If ``snapshots.nc`` was created by the command above, running

.. code-block:: none

   pismr -i snapshots.nc -y 1000 -o output_2.nc

will initialize using the last record in the file, at about :math:`5000` years. By
contrast, to restart from :math:`1500` years (for example) it is necessary to extract the
corresponding record using ``ncks``

.. code-block:: none

   ncks -d t,1500years snapshots.nc foo.nc

and then restart from ``foo.nc``. Note that ``-d t,N`` means "extract the :math:`N`-th
record" (counting from zero). So, this command is equivalent to

.. code-block:: none

   ncks -d t,1 snapshots.nc foo.nc

Also note that the second snapshot will probably be *around* :math:`1500` years and
``ncks`` handles this correctly: it takes the record closest to :math:`1500` years.

By default re-startable snapshots contain only the variables needed for restarting PISM.
Use the command-line option ``-save_size`` to change what is saved.

Another possible use of snapshots is for restarting runs on a batch system which kills
jobs which go over their allotted time. Running PISM with options ``-y 1500``
``-save_times 1000:100:1400`` would mean that if the job is killed before completing the
whole 1500 year run, we can restart from near the last multiple of :math:`100` years.
Restarting with option ``-ye`` would finish the run on the desired year.

When running PISM on such a batch system it is also possible to save re-startable
snapshots at equal wall-clock time (as opposed to model time) intervals by adding the
":opt:`-backup_interval` (hours)" option.

.. caution::

   If the wall-clock limit is equal to :math:`N` times backup interval for a whole number
   :math:`N` PISM will likely get killed while writing the last backup.

It is also possible to save snapshots to separate files using the ``-save_split`` option.
For example, the run above can be changed to

.. code-block:: none

   pismr -i foo.nc -y 10000 -o output.nc -save_file snapshots \
         -save_times 1000,1500,2000,5000 -save_split

for this purpose. This will produce files called ``snapshots-year.nc``. This option is
generally faster if many snapshots are needed, apparently because of the time necessary to
reopen a large file at each snapshot when ``-save_split`` is not used. Note that tools
like NCO and ``ncview`` usually behave as desired with wildcards like
"``snapshots-*.nc``".

:numref:`tab-snapshot-opts` lists the options related to saving snapshots of the
model state.

.. list-table:: Command-line options controlling saving snapshots of the model state.
   :name: tab-snapshot-opts
   :header-rows: 1

   * - Option
     - Description
   * - :opt:`-save_file`
     - Specifies the file to save to.
   * - :opt:`-save_times`
     - Specifies times at which to save snapshots, by either a MATLAB-style range
       :math:`a:\Delta t:b` or a comma-separated list.
   * - :opt:`-save_split`
     - Separate the snapshot output into files named ``snapshots-year.nc``. Faster if you
       are saving more than a dozen or so snapshots.
   * - :opt:`-save_size` ``[none,small,medium,big,big_2d]``
     - Similar to ``o_size``, changes the "size" of the file (or files) written; the
       default is "small"

.. _sec-diagnostic-viewers:

Run-time diagnostic viewers
---------------------------

Basic graphical views of the changing state of a PISM ice model are available at the
command line by using options listed in :numref:`tab-diag-viewers`. All the
quantities listed in tables |diagnostics| are available. Additionally, a couple of
diagnostic quantities are *only* available as run-time viewers; these are shown in table
:numref:`tab-special-diag-viewers`.

.. list-table:: Options controlling run-time diagnostic viewers
   :name: tab-diag-viewers
   :header-rows: 1

   * - Option
     - Description
   * - :opt:`-view`
     - Turns on map-plane views of one or several variables, see tables FIXME
   * - :opt:`-view_size` (number)
     - desired viewer size, in pixels
   * - :opt:`-display`
     - The option ``-display :0`` seems to frequently be needed to let PETSc use Xwindows
       when running multiple processes. It must be given as a *final* option, after all
       the others.

The option ``-view`` shows map-plane views of 2D fields and surface and basal views of 3D
fields (see tables |diagnostics|); for example:

.. code-block:: none

   pismr -i input.nc -y 1000 -o output.nc -view thk,tempsurf

shows ice thickness and ice temperature at the surface.

.. list-table:: Special run-time-only diagnostic viewers
   :name: tab-special-diag-viewers
   :header-rows: 1

   * - Option
     - Description
   * - :opt:`-ssa_view_nuh`
     - log base ten of ``nuH``, only available if the finite-difference SSA solver is
       active.
   * - :opt:`-ssa_nuh_viewer_size` (number)
     - Adjust the viewer size.
   * - :opt:`-ksp_monitor_draw`
     - Iteration monitor for the Krylov subspace routines (KSP) in PETSc. Residual norm
       versus iteration number.

.. _sec-pism-defaults:

PISM's configuration flags and parameters, and how to change them
-----------------------------------------------------------------

PISM's behavior depends on values of many flags and physical parameters (see `PISM Source
Code Browser <pism-browser_>`_ for details). Most of parameters have default values [#]_
which are read from the configuration file ``pism_config.nc`` in the ``lib``
sub-directory.

It is possible to run PISM with an alternate configuration file using the :opt:`-config`
command-line option:

.. code-block:: none

   pismr -i foo.nc -y 1000 -config my_config.nc

The file ``my_config.nc`` has to contain *all* of the flags and parameters present in
``pism_config.nc``.

The list of parameters is too long to include here; please see the `PISM Source Code
Browser`_ for an automatically-generated table describing them.

Some command-line options *set* configuration parameters; some PISM executables have
special parameter defaults. To examine what parameters were used in a particular run, look
at the attributes of the ``pism_config`` variable in a PISM output file.

.. _sec-parameter-studies:

Managing parameter studies
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Keeping all PISM output files in a parameter study straight can be a challenge. If the
parameters of interest were controlled using command-line options then one can use
``ncdump -h`` and look at the ``history`` global attribute.

Alternatively, one can change parameter values by using an "overriding" configuration
file. The :opt:`-config_override` command-line option provides this alternative. A file
used with this option can have a subset of the configuration flags and parameters present
in ``pism_config.nc``. Moreover, PISM adds the ``pism_config`` variable with values used
in a run to the output file, making it easy to see which parameters were used.

Here's an example. Suppose we want to compare the dynamics of an ice-sheet on Earth to the
same ice-sheet on Mars, where the only physical change was to the value of the
acceleration due to gravity. Running

.. code-block:: none

   pismr -i input.nc -y 1e5 -o earth.nc <other PISM options>

produces the "Earth" result, since PISM's defaults correspond to this planet. Next, we
create ``mars.cdl`` containing the following:

.. code-block:: none

   netcdf mars {
       variables:
       byte pism_overrides;
       pism_overrides:constants.standard_gravity = 3.728;
       pism_overrides:constants.standard_gravity_doc = "m s-2; standard gravity on Mars";
   }


Notice that the variable name is ``pism_overrides`` and not ``pism_config`` above. Now

.. code-block:: none

   ncgen -o mars_config.nc mars.cdl
   pismr -i input.nc -y 1e5 -config_override mars_config.nc -o mars.nc <other PISM options>

will create ``mars.nc``, the result of the "Mars" run. Then we can use ``ncdump`` to see
what was different about ``mars.nc``:

.. code-block:: diff

   ncdump -h earth.nc | grep pism_config: > earth_config.txt
   ncdump -h mars.nc | grep pism_config: > mars_config.txt
   diff -U 1 earth_config.txt mars_config.txt
   --- earth_config.txt	2015-05-08 12:44:43.000000000 -0800
   +++ mars_config.txt	2015-05-08 12:44:51.000000000 -0800
   @@ -734,3 +734,3 @@
                   pism_config:ssafd_relative_convergence_units = "1" ;
   -               pism_config:constants.standard_gravity_doc = "acceleration due to gravity on Earth geoid" ;
   +               pism_config:constants.standard_gravity_doc = "m s-2; standard gravity on Mars" ;
                   pism_config:constants.standard_gravity_type = "scalar" ;
   @@ -1057,3 +1057,3 @@
                   pism_config:ssafd_relative_convergence = 0.0001 ;
   -               pism_config:constants.standard_gravity = 9.81 ;
   +               pism_config:constants.standard_gravity = 3.728 ;
                   pism_config:start_year = 0. ;

.. _sec-saving-pism-config:

Saving PISM's configuration for post-processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to saving ``pism_config`` in the output file, PISM automatically adds this
variable to all files it writes (snap shots, time series of scalar and spatially-varying
diagnostic quantities, and backups). This may be useful for post-processing and analysis
of parameter sties as the user has easy access to all configuration options, model
choices, etc., without the need to keep run scripts around.

.. _sec-regridding:

Regridding
----------

It is common to want to interpolate a coarse grid model state onto a finer grid or vice
versa. For example, one might want to do the EISMINT II experiment on the default grid,
producing output ``foo.nc``, but then interpolate both the ice thickness and the
temperature onto a finer grid. The basic idea of "regridding" in PISM is that one starts
over from the beginning on the finer grid, but one extracts the desired variables stored
in the coarse grid file and interpolates these onto the finer grid before proceeding with
the actual computation.

The transfer from grid to grid is reasonably general --- one can go from coarse to fine or
vice versa in each dimension :math:`x,y,z` --- but the transfer must always be done by
*interpolation* and never *extrapolation*. (An attempt to do the latter will always
produce a PISM error.)

Such "regridding" is done using the :opt:`-regrid_file` and :opt:`-regrid_vars` commands
as in this example: }

.. code-block:: none

    pisms -eisII A -Mx 101 -My 101 -Mz 201 -y 1000 \
          -regrid_file foo.nc -regrid_vars thk,temp -o bar.nc

By specifying regridded variables "``thk,temp``", the ice thickness and temperature values
from the old grid are interpolated onto the new grid. Here one doesn't need to regrid the
bed elevation, which is set identically zero as part of the EISMINT II experiment A
description, nor the ice surface elevation, which is computed as the bed elevation plus
the ice thickness at each time step anyway.

A slightly different use of regridding occurs when "bootstrapping", as described in
section :ref:`sec-initboot` and illustrated by example in section :ref:`sec-start`.

See :numref:`tab-regridvar` for the regriddable variables using ``-regrid_file``.
Only model state variables are regriddable, while climate and boundary data generally are
not explicitly regriddable. (Bootstrapping, however, allows the same general interpolation
as this explicit regrid.)

.. list-table:: Regriddable variables.  Use ``-regrid_vars`` with these names.
   :header-rows: 1
   :name: tab-regridvar

   * - Name
     - Description
   * - :var:`age`
     - age of ice
   * - :var:`bwat`
     - effective thickness of subglacial melt water
   * - :var:`bmelt`
     - basal melt rate
   * - :var:`dbdt`
     - bedrock uplift rate
   * - :var:`litho_temp`
     - lithosphere (bedrock) temperature
   * - :var:`mask`
     - grounded/dragging/floating integer mask, see section :ref:`sec-floatmask`
   * - :var:`temp`
     - ice temperature
   * - :var:`thk`
     - land ice thickness
   * - :var:`topg`
     - bedrock surface elevation
   * - :var:`enthalpy`
     - ice enthalpy

Here is another example: suppose you have an output of a PISM run on a fairly coarse grid
(stored in ``foo.nc``) and you want to continue this run on a finer grid. This can be done
using ``-regrid_file`` along with ``-bootstrap``:

.. code-block:: none

   pismr -i foo.nc -bootstrap -Mx 201 -My 201 -Mz 21 -Lz 4000 \
         -regrid_file foo.nc -regrid_vars litho_temp,enthalpy -y 100 -o bar.nc \
         -surface constant

In this case all the model-state 2D variables present in ``foo.nc`` will be interpolated
onto the new grid during bootstrapping, which happens first, while three-dimensional
variables are filled using heuristics mentioned in section :ref:`sec-initboot`. Then
temperature in bedrock (``litho_temp``) and ice enthalpy (``enthalpy``) will be
interpolated from ``foo.nc`` onto the new grid during the regridding stage, overriding
values set at the bootstrapping stage. All of this, bootstrapping and regridding, occurs
before the first time step.

By default PISM checks the grid overlap and stops if the current computational domain is
not a subset of the one in a ``-regrid_file``. It is possible to disable this check and
allow constant extrapolation: use the option :opt:`-allow_extrapolation`.

For example, in a PISM run the ice thickness has to be lower than the vertical extent of
the computational domain. If the ice thickness exceeds ``Lz`` PISM saves the model state
and stops with an error message.

.. code-block:: none

   pismr -i input.nc -bootstrap -Mz 11 -Lz 1000 -z_spacing equal \
         -y 3e3 \
         -o too-short.nc
   PISM ERROR: Ice thickness exceeds the height of the computational box (1000.0000 m).
               The model state was saved to 'too-short_max_thickness.nc'.
               To continue this simulation, run with
               -i too-short_max_thickness.nc -bootstrap -regrid_file too-short_max_thickness.nc \
               -allow_extrapolation -Lz N [other options]
               where N > 1000.0000.

Regridding with extrapolation makes it possible to extend the vertical grid and continue a
simulation like this one --- just follow the instructions provided in the error message.

.. |pid| replace:: *PID*\s

.. _sec-signal:

Signals, to control a running PISM model
----------------------------------------

Ice sheet model runs sometimes take a long time, so the state of a run may need checking.
Sometimes the run needs to be stopped, but with the possibility of restarting. PISM
implements these behaviors using "signals" from the POSIX standard, included in Linux and
most flavors of Unix. :numref:`tab-signals` summarizes how PISM responds to signals.
A convenient form of ``kill``, for Linux users, is ``pkill`` which will find processes by
executable name. Thus "``pkill -USR1 pismr``" might be used to send all PISM processes the
same signal, avoiding an explicit list of |pid|.

.. list-table:: Signalling running PISM processes.  "|pid|" stands for list of all identifiers of the PISM processes.
   :name: tab-signals
   :header-rows: 1

   * - Command
     - Signal
     - PISM behavior
   * - ``kill -KILL`` |pid|
     - ``SIGKILL``
     - Terminate with extreme prejudice. PISM cannot catch it and no state is saved.
   * - ``kill -TERM`` |pid|
     - ``SIGTERM``
     - End process(es), but save the last model state in the output file, using ``-o``
       name or default name as normal. Note that the ``history`` string in the output file
       will contain an "``EARLY EXIT caused by signal SIGTERM``" indication.
   * - ``kill -USR1`` |pid|
     - ``SIGUSR1``
     - Process(es) will continue after saving the model state at the end of the current
       time step, using a file name including the current model year. Time-stepping is not
       altered. Also flushes output buffers of scalar time-series.
   * - ``kill -USR2`` |pid|
     - ``SIGUSR2``
     - Just flush time-series output buffers.
   
Here is an example. Suppose we start a long verification run in the background, with
standard out redirected into a file:

.. code-block:: none

   pismv -test G -Mz 101 -y 1e6 -o testGmillion.nc >> log.txt &

This run gets a Unix process id, which we assume is "8920". (Get it using ``ps`` or
``pgrep``.) If we want to observe the run without stopping it we send the ``USR1`` signal:


.. code-block:: none

   kill -USR1 8920

(With ``pkill`` one can usually type "``pkill -usr1 pismv``".) Suppose it happens that we
caught the run at year 31871.5. Then, for example, a NetCDF file ``pismv-31871.495.nc`` is
produced. Note also that in the standard out log file ``log.txt`` the line

.. code-block:: none

   caught signal SIGUSR1:  Writing intermediate file ... and flushing time series.

appears around that time step. Suppose, on the other hand, that the run needs to be
stopped. Then a graceful way is

.. code-block:: none

   kill -TERM 8920

because the model state is saved and can be inspected.

.. _sec-adapt:

Understanding adaptive time-stepping
------------------------------------

At each time step the PISM standard output includes "flags" and then a summary of the
model state using a few numbers. A typical example is

.. code-block:: none

   v$Eh  diffusivity (dt=0.83945 in 2 substeps; av dt_sub_mass_cont=0.41972)
   S -124791.571:  3.11640   2.25720      3.62041    18099.93737
   y  SSA:     3 outer iterations, ~17.0 KSP iterations each

The characters "``v$Eh``" at the beginning of the flags line, the first line in the above
example, give a very terse description of which physical processes were modeled in that
time step. Here "``v``" means that a stress balance was solved to compute the velocity.
Then the enthalpy was updated ("``E``") and the ice thickness and surface elevation were
updated ("``h``"). The rest of the flags line looks like

.. code-block:: none

   diffusivity (dt=0.83945 in 2 substeps; av dt_sub_mass_cont=0.41972)

Recall that the PISM time step is determined by an adaptive mechanism. Stable mass
conservation and conservation of energy solutions require such an adaptive time-stepping
scheme :cite:`BBL`. The first character we see here, namely "``diffusivity``", is the
adaptive-timestepping "reason" flag. See :numref:`tab-adaptiveflag`. We also see
that there was a major time step of :math:`0.83945` model years divided into :math:`2`
substeps of about :math:`0.42` years. The :opt:`-skip` option enables this mechanism,
while :opt:`-skip_max` sets the maximum number of such substeps. The adaptive mechanism
may choose to take fewer substeps than ``-skip_max`` so as to satisfy certain numerical
stability criteria, however.

The second line in the above, the line which starts with "``S``", is the summary. Its
format, and the units for these numbers, is simple and is given by a couple of lines
printed near the beginning of the standard output for the run:

.. code-block:: none

   P       YEAR:       ivol      iarea  max_diffusivity  max_hor_vel
   U      years   10^6_km^3  10^6_km^2         m^2 s^-1       m/year

That is, in each summary we have the total ice volume, total ice area, maximum diffusivity
(of the SIA mass conservation equation), and maximum horizontal velocity (i.e.
:math:`\max(\max(|u|), \max(|v|))`).

The third line of the above example shows that the SSA stress balance was solved.
Information on the number of nonlinear (outer) and linear (inner) iterations is provided
:cite:`BBssasliding`.

.. list-table:: Meaning of the adaptive time-stepping "reason" flag in the standard output
                flag line.
   :header-rows: 1
   :name: tab-adaptiveflag

   * - PISM output
     - Active adaptive constraint or PISM sub-system that limited time-step size

   * - ``3D CFL``
     - three-dimensional CFL for temperature/age advection :cite:`BBL`

   * - ``diffusivity``
     - diffusivity for SIA mass conservation :cite:`BBL`, :cite:`HindmarshPayne`

   * - ``end of the run``
     - end of prescribed run time

   * - ``max``
     - maximum allowed :math:`\Delta t` applies; set with ``-max_dt``

   * - ``internal (derived class)``
     - maximum :math:`\Delta t` was temporarily set by a derived class

   * - ``2D CFL``
     - 2D CFL for mass conservation in SSA regions (upwinded; :cite:`BBssasliding`)

   * - ``-ts_... reporting``
     - the ``-ts_times`` option and the configuration flag
       :config:`time_stepping.hit_ts_times`; see section :ref:`sec-saving-time-series`

   * - ``-extra_... reporting``
     - the ``-extra_times`` option; see section :ref:`sec-saving-diagnostics`

   * - ``surface``
     - a surface or an atmosphere model

   * - ``ocean``
     - an ocean model

   * - ``hydrology``
     - a hydrology model stability criterion, see section :ref:`sec-subhydro`

   * - ``BTU``
     - time-the bedrock thermal layer model, see section :ref:`sec-energy`

   * - ``eigencalving``
     - the eigen-calving model, see section :ref:`sec-calving`

.. list-table:: Options controlling time-stepping
   :header-rows: 1
   :name: tab-time-stepping

   * - Option
     - Description

   * - :opt:`-adapt_ratio`
     - Adaptive time stepping ratio for the explicit scheme for the mass balance equation.

   * - :opt:`-max_dt` (years)
     - The maximum time-step in years. The adaptive time-stepping scheme will make the
       time-step shorter than this as needed for stability, but not longer.

   * - :opt:`-skip`
     - Enables time-step skipping, see below.

   * - :opt:`-skip_max`
     - Number of mass-balance steps, including SIA diffusivity updates, to perform before
       temperature, age, and SSA stress balance computations are done. This is only
       effective if the time step is being limited by the diffusivity time step
       restriction associated to mass continuity using the SIA. The maximum recommended
       value for ``-skip_max`` is, unfortunately, dependent on the context. The
       temperature field should be updated when the surface changes significantly, and
       likewise the basal sliding velocity if it comes (as it should) from the SSA
       calculation.

   * - :opt:`-timestep_hit_multiples` (years)
     - Hit multiples of the number of model years specified. For example, if stability
       criteria require a time-step of 11 years and the ``-timestep_hit_multiples 3``
       option is set, PISM will take a 9 model year long time step. This can be useful to
       enforce consistent sampling of periodic climate data.

.. _sec-petscoptions:

PETSc options for PISM users
----------------------------

All PETSc programs including PISM accept command line options which control how PETSc
distributes jobs among parallel processors, how it solves linear systems, what additional
information it provides, and so on. The PETSc manual :cite:`petsc-user-ref` is the complete
reference on these options. We list some here that are useful to PISM users. They can be
mixed in any order with PISM options.

Both for PISM and PETSc options, there are ways of avoiding the inconvenience of long
commands with many runtime options. Obviously, and as illustrated by examples in the
previous sections, shell scripts can be set up to run PISM. But PETSc also provides two
mechanisms to give runtime options without retyping at each run command.

First, the environment variable ``PETSC_OPTIONS`` can be set. For example, a sequence of
runs might need the same refined grid, and you might want to know if other options are
read, ignored, or misspelled. Set (in Bash):

.. code-block:: none

   export PETSC_OPTIONS="-Mx 101 -My 101 -Mz 51 -options_left"

The runs

.. code-block:: none

   pismv -test F -y 100
   pismv -test G -y 100

then have the same refined grid in each run, and the runs report on which options were
read.

Alternatively, the file ``.petscrc`` is always read, if present, from the directory where
PISM (i.e. the PETSc program) is started. It can have a list of options, one per line. In
theory, these two PETSc mechanisms (``PETSC_OPTIONS`` and ``.petscrc``) can be used
together.

.. "-da_processors_x M -da_processors_y N" should not be documented here because they do not work. the reason is that IceModelVec2 and IceModelVec3 put the Mx, My dimensions in different arguments to the DACreate commands (FIXME: I don't think this is true.)

Now we address controls on how PETSc solves systems of linear equations, which uses the
PETSc "KSP" component (Krylov methods). Such linear solves are needed each time the
nonlinear SSA stress balance equations are used (e.g. with the option ``-stress_balance
ssa -ssa_method fd``).

Especially for solving the SSA equations with high resolution on multiple processors, it
is recommended that the option :opt:`-ssafd_ksp_rtol` be set lower than its default value
of :math:`10^{-5}`. For example,


.. code-block:: none

   mpiexec -n 8 ssa_testi -Mx 3 -My 769 -ssa_method fd

may fail to converge on a certain machine, but adding "``-ssafd_ksp_rtol 1e-10``" works
fine.

There is also the question of solver *type*, using option :opt:`-ssafd_ksp_type`. Based on
one processor evidence from ``ssa_testi``, the following are possible choices in the sense
that they work and allow convergence at some reasonable rate: ``cg``, ``bicg``, ``gmres``,
``bcgs``, ``cgs``, ``tfqmr``, ``tcqmr``, and ``cr``. It appears ``bicg``, ``gmres``,
``bcgs``, and ``tfqmr``, at least, are all among the best. The default is ``gmres``.

Actually the KSP uses preconditioning. This aspect of the solve is critical for parallel
scalability, but it gives results which are dependent on the number of processors. The
preconditioner type can be chosen with :opt:`-ssafd_pc_type`. Several choices are
possible, but for solving the ice stream and shelf equations we recommend only
``bjacobi``, ``ilu``, and ``asm``. Of these it is not currently clear which is fastest;
they are all about the same for ``ssa_testi`` with high tolerances (e.g. ``-ssa_rtol
1e-7`` ``-ssafd_ksp_rtol 1e-12``). The default (as set by PISM) is ``bjacobi``. To force
no preconditioning, which removes processor-number-dependence of results but may make the
solves fail, use ``-ssafd_pc_type none``.

For the full list of PETSc options controlling the SSAFD solver, run

.. code-block:: none

   ssa_testi -ssa_method fd -help | grep ssafd_ | less

.. _sec-scripts:

Utility and test scripts
------------------------

In the ``test/`` and ``util/`` subdirectories of the PISM directory the user will find
some python scripts and one Matlab script, listed in :numref:`tab-scripts-overview`.
The python scripts are all documented at the *Packages* tab on the `PISM Source Code
Browser`_. The Python scripts all take option ``--help``.

.. list-table:: Some scripts which help in using PISM
   :name: tab-scripts-overview
   :header-rows: 1

   * - Script
     - Function
   * - ``test/vfnow.py``
     - Organizes the process of verifying PISM. Specifies standard refinement paths for
       each of the tests (section :ref:`sec-verif`).
   * - ``test/vnreport.py``
     - Automates the creation of convergence graphs like figures :numref:`fig-thickerrsB`
       -- :numref:`fig-velerrsI`.
   * - ``util/fill_missing.py``
     - Uses an approximation to Laplace's equation :math:`\nabla^2 u = 0` to smoothly
       replace missing values in a two-dimensional NetCDF variable. The "hole" is filled
       with an average of the boundary non-missing values. Depends on ``netcdf4-python``
       and ``scipy`` Python packages.
   * - ``util/flowline.py``
     - See section :ref:`sec-flowline-modeling`.
   * - ``util/flowlineslab.py``
     - See section :ref:`sec-flowline-modeling`.
   * - ``util/check_stationarity.py``
     - Evaluate stationarity of a variable in a PISM ``-ts_file`` output.
   * - ``util/nc2cdo.py``
     - Makes a netCDF file ready for Climate Data Operators (CDO).
   * - ``util/nc2mat.py``
     - Reads specified variables from a NetCDF file and writes them to an output file in
       the MATLAB binary data file format ``.mat``, supported by MATLAB version 5 and
       later. Depends on ``netcdf4-python`` and ``scipy`` Python packages.
   * - ``util/nccmp.py``
     - A script comparing variables in a given pair of NetCDF files; used by PISM software
       tests.
   * - ``util/pism_config_editor.py``
     - Makes modifying or creating PISM configuration files easier.
   * - ``util/pism_matlab.m``
     - An example MATLAB script showing how to create a simple NetCDF file PISM can
       bootstrap from.
   * - ``util/PISMNC.py``
     - Used by many Python example scripts to generate a PISM-compatible file with the
       right dimensions and time-axis.

.. _sec-flowline-modeling:

Using PISM for flow-line modeling
---------------------------------

As described in sections :ref:`sec-coords` and :ref:`sec-grid`, PISM is a
three-dimensional model. Moreover, parameters ``Mx`` and ``My`` have to be greater than or
equal to three, so it is not possible to turn PISM into a 2D (flow-line) model by setting
``Mx`` or ``My`` to 1.

There is a way around this, though: by using the :opt:`-periodicity` option to tell PISM
to make the computational grid :math:`y`-periodic and providing initial and boundary
conditions that are functions of :math:`x` only one can ensure that there is no flow in
the :math:`y`\-direction. (Option :opt:`-periodicity` takes an argument specifying the
direction: ``none``, ``x``, ``y`` and ``xy`` --- for "periodic in both X- and
Y-directions".)

In this case ``Mx`` can be any number; we want to avoid unnecessary computations, though,
so "``-Mx 3``" is the obvious choice.

One remaining problem is that PISM still expects input files to contain both ``x`` and
``y`` dimensions. To help with this, PISM comes with a Python script ``flowline.py`` that
turns NetCDF files with :math:`N` grid points along a flow line into files with 2D fields
containing :math:`N\times3` grid points. [#]_

Here's an example which uses the script ``util/flowlineslab.py`` to create a minimal, and
obviously unrealistic, dataset. A file ``slab.nc`` is created by ``util/flowlineslab.py``,
but it is not ready to use with PISM. Proceed as follows, after checking that ``util/`` is
on your path:

.. code-block:: bash

   flowlineslab.py                         # creates slab.nc with only an x-direction
   flowline.py -o slab-in.nc --expand -d y slab.nc


produces a PISM-ready ``slab-in.nc``. Specifically, ``flowline.py`` "expands" its input
file in the y-direction. Now we can "bootstrap" from ``slab-in.nc``:

.. code-block:: none

   mpiexec -n 2 pismr -surface given -i slab-in.nc -bootstrap -periodicity y \
           -Mx 201 -My 3 -Lx 1000 -Ly 4 -Lz 2000 -Mz 11 -y 10000 -o pism-out.nc

To make it easier to visualize data in the file created by PISM, "collapse" it:

.. code-block:: none

   flowline.py -o slab-out.nc --collapse -d y pism-out.nc

.. _sec-code-modifications:

Managing source code modifications
----------------------------------

"Practical usage" may include editing the source code to extend, fix or replace parts of
PISM.

We provide both user-level (this manual) and developer-level documentation. Please see
source code browsers at |pism-docs| for the latter.

- To use your (modified) version of PISM, you will need to follow the compilation from
  sources instructions in the *Installation Manual*
- We find it very useful to be able to check if a recent source code change broke
  something. PISM comes with "regression tests", which check if certain parts of PISM
  perform the way it should. [#]_

  Run "``make test``" in the build directory to run PISM's regression tests.

  Note, though, that while a test failure usually means that the new code needs more work,
  passing all the tests does not guarantee that everything works as it should. We are
  constantly adding new tests, but so far only a subset of PISM's functionality can be
  tested automatically.
- We strongly recommend using a version control system to manage code changes. Not only is
  it safer than the alternative, it is also more efficient.

.. rubric:: Footnotes

.. [#] For example, resolutions of 2km and higher on the whole-Greenland scale.
.. [#] This is not likely to change.
.. [#] This will be fixed in a later release of PISM.
.. [#] For ``pismr``, grid parameters ``Mx``, ``My``, ``Mz``, ``Mbz``, ``Lz``, ``Lbz``,
       that must be set at bootstrapping, are exceptions.
.. [#] This script requires the ``numpy`` and ``netCDF4`` Python modules. Run
       ``flowline.py --help`` for a full list of options.
.. [#] This automates running verification tests described in section :ref:`sec-verif`,
       for example.
