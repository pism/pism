.. default-role:: math

.. _sec-modeling-computational:

Modeling choices: Grid and time
===============================

.. _sec-coords:

Computational box
-----------------

.. FIXME: This assumes that the sea (i.e. floatation) level is zero.

PISM does all simulations in a computational box which is rectangular in the PISM
coordinates. The coordinate system has horizontal coordinates `x,y` and a vertical
coordinate `z`. The `z` coordinate is measured positive upward from the base
of the ice. The vector of gravity is in the negative `z` direction. The surface
`z=0` is the base of the ice, however, and thus is usually not horizontal in the
sense of being parallel to the geoid. The surface `z=0` is the base of the ice both
when the ice is grounded and when the ice is floating.

When the ice is grounded, the true physical vertical coordinate `z'`, namely the
coordinate measure relative to a reference geoid, is given by `z'=z+b(x,y)` where
`b(x,y)` is the bed topography. The surface `z'=h(x,y)` is the surface of the
ice. In the grounded case the equation `h(x,y)=H(x,y)+b(x,y)` always applies if
`H(x,y)` is the thickness of the ice.

In the floating case, the physical vertical coordinate is

.. math::

   z'=z-\frac{\rho_i}{\rho_s} H(x,y)

where `\rho_i` is the density of ice and `\rho_s` the density of sea
water. Again `z=0` is the base of the ice, which is the surface

.. math::

   z' = -\frac{\rho_i}{\rho_s} H(x,y).

The surface of the ice is

.. math::

   h(x,y) = \left(1-\frac{\rho_i}{\rho_s}\right) H(x,y).

The *flotation criterion* `-\frac{\rho_i}{\rho_s} H(x,y) > b(x,y)` applies.

The computational box can extend downward into the bedrock. As `z=0` is the base of
the ice, the bedrock corresponds to negative `z` values regardless of its true (i.e.
`z'`) elevation.

The extent of the computational box, along with its bedrock extension downward, is
determined by four numbers ``Lx``, ``Ly``, ``Lz``, and ``Lbz`` (see Figure
:numref:`fig-rectilinearbox` and Table :numref:`tab-compbox`). The first two of these are
half-widths and have units of kilometers when set by options or displayed.

.. figure:: figures/rectilinearbox.png
   :name: fig-rectilinearbox

   PISM's computational box

.. list-table:: Options defining the extent of PISM's computational box
   :name: tab-compbox
   :header-rows: 1
   :widths: 20, 80

   * - Option
     - Description
   * - :opt:`-Lx` (km)
     - Half-width of the computational domain (in the `x`\-direction) 
   * - :opt:`-Ly` (km)
     - Half-width of the computational domain (in the `y`\-direction) 
   * - :opt:`-Lz` (meters)
     - Height of the computational domain; must exceed maximum ice thickness 
   * - :opt:`-Lbz` (meters)
     - Depth of the computational domain in the bedrock thermal layer

.. _sec-grid:
    
Spatial grid
------------

The PISM grid covering the computational box is equally spaced in horizontal (`x`
and `y`) directions. Vertical spacing in the ice is quadratic by default but
optionally equal spacing can be chosen; choose with options :opt:`-z_spacing`
[``quadratic``, ``equal``] at bootstrapping. The grid read from a "``-i``" input file is
used as is. The bedrock thermal layer model always uses equal vertical spacing.

The grid is described by four numbers, namely the number of grid points ``Mx`` in the
`x` direction, the number ``My`` in the `y` direction, the number ``Mz`` in
the `z` direction within the ice, and the number ``Mbz`` in the `z` direction
within the bedrock thermal layer. These are specified by options :opt:`-Mx`, :opt:`-My`,
:opt:`-Mz`, and :opt:`-Mbz`, respectively. The defaults are 61, 61, 31, and 1,
respectively. Note that ``Mx``, ``My``, ``Mz``, and ``Mbz`` all indicate the number of
grid *points* so the number of grid *spaces* are one less, thus 60, 60, 30, and 0 in the
default case.

The lowest grid point in a column of ice, at `z=0`, coincides with the highest grid
point in the bedrock, so ``Mbz`` must always be at least one. Choosing ``Mbz```>1`
is required to use the bedrock thermal model. When a thermal bedrock layer is used, the
distance ``Lbz`` is controlled by the ``-Lbz`` option. Note that ``Mbz`` is unrelated to
the bed deformation model (glacial isostasy model); see section :ref:`sec-beddef`.

In the quadratically-spaced case the spacing near the ice/bedrock interface is about four
times finer than it would be with equal spacing for the same value of ``Mz``, while the
spacing near the top of the computational box is correspondingly coarser. For a detailed
description of the spacing of the grid, see the documentation on
``IceGrid::compute_vertical_levels()`` in the PISM class browser.

The user should specify the grid when using ``-bootstrap`` or when initializing a
verification test (section :ref:`sec-verif`) or a simplified-geometry experiment (section
:ref:`sec-simp`). If one initializes PISM from a saved model state using ``-i`` then the
input file determines all grid parameters. For instance, the command

.. code-block:: none

   pismr -i foo.nc -y 100

should work fine if ``foo.nc`` is a PISM output file. Because ``-i`` input files take
precedence over options,

.. code-block:: none

   pismr -i foo.nc -Mz 201 -y 100

will give a warning that "``PISM WARNING: ignoring command-line option '-Mz'``".

.. _sec-domain-dstribution:

Parallel domain distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^


When running PISM in parallel with ``mpiexec -n N``, the horizontal grid is distributed
across `N` processes [#]_. PISM divides the grid into `N_x` parts in the
`x` direction and `N_y` parts in the `y` direction. By default this is
done automatically, with the goal that `N_x\times N_y = N` and `N_x` is as
close to `N_y` as possible. Note that `N` should, therefore, be a composite
(not prime) number.

Users seeking to override this default can specify `N_x` and `N_y` using the
:opt:`-Nx` and :opt:`-Ny` command-line options.

Once `N_x` and `N_y` are computed, PISM computes sizes of sub-domains
`M_{x,i}` so that `\sum_{i=1}^{N_x}M_{x,i} = \mathrm{Mx}` and `M_{x,i} -
\left\lfloor \mathrm{Mx} / N_x \right\rfloor < 1`. To specify strip widths `M_{x,i}`
and `M_{y,i}`, use command-line options :opt:`-procs_x` and :opt:`-procs_y`. Each
option takes a comma-separated list of numbers as its argument. For example,

.. code-block:: none

   mpiexec -n 3 pisms -Mx 101 -My 101 -Nx 1 -Ny 3 -procs_x 101 -procs_y 20,61,20

splits a `101 \times 101` grid into 3 strips along the `x` axis.

To see the parallel domain decomposition from a completed run, see the ``rank`` variable
in the output file, e.g. using ``-o_size big``. The same ``rank`` variable is available as
a spatial diagnostic field (subsection :ref:`sec-saving-spat-vari`).

.. _sec-time:

Model time
----------


Table :numref:`tab-timeoptions` gives the command-line options which control PISM time. If
option ``-ys`` is absent then the start year is read from the input file (if present) or
it defaults to zero. The default value for the end year is the start year plus the given
(``-y``) run length. If both ``-ys`` and ``-ye`` are used then the run length is set to
the difference. Using all three of ``-ys``, ``-y`` and ``-ys`` is not allowed; this
generates an error message.

.. list-table:: Command-line options controlling PISM time
   :name: tab-timeoptions
   :header-rows: 1
   :widths: 20, 80

   * - Option
     - Description
   * - :opt:`-y` (years)
     - Number of model years to run.
   * - :opt:`-ys` (years)
     - Model year at which to start the run. Also resets the model
       time, ignoring any time in the input file.
   * - :opt:`-ye` (years)
     - Model year at which to end the run.

.. _sec-calendars:

Calendars
---------


Most of PISM, and its ice dynamics core in particular, only needs to know the length of
the current time-step. Internally PISM stores time in "seconds since a specified moment"
and thus PISM generally does not use or need a calendar. [#]_ We refer to PISM internal
time as *model time*.

One can select a calendar for more precise control of the model time, however. A
"calendar" is a concept that is part of the `CF Conventions`_. Choosing a calendar is
appropriate for runs for specific temporal periods like "the 18th-century" or
"1989--2010". The calendar is generally needed because specific knowledge of lengths of
months and years is required to use climate data properly or to facilitate model
validation.

PISM uses CalCalcs_ by David W. Pierce to perform calendric computations. This lets us
support all the `calendars <CF-Calendars_>`_ defined by the CF Metadata Conventions
document except for the ``366_day`` (``all_leap``) calendar.

Time units in PISM's output files always contain a reference date because it is required
by the CF metadata conventions.

By default PISM does not use a calendar. This is appropriate for runs that do not require
precise application of forcing data or reporting on particular dates (paleo-climate runs,
for example). In this mode PISM ignores the reference date in time unit specifications
(such as "``days since 1969-7-20``"), though the value set using
:config:`time.reference_date` configuration parameter is saved in (is passed forward into)
output files.

.. list-table:: Calendars supported by PISM. Please see CalCalcs_ documentation for
                details
   :name: tab-calendars
   :header-rows: 1
   :widths: 20, 80

   * - Keyword
     - Meaning
   * - ``gregorian`` or ``standard``
     - Mixed Gregorian/Julian calendar used today.
   * - ``proleptic_gregorian``
     - Gregorian calendar extended to dates before 1582-10-15.
   * - ``noleap`` or ``365_day``
     - Calendar with fixed-length 365-day years
   * - ``360_day``
     - Calendar with fixed-length 360-day years divided into 30-day months
   * - ``julian``
     - Julian calendar 
   * - ``none``
     - no calendar

Selecting a calendar using the :config:`time.calendar` configuration parameter or the
:opt:`-calendar` command-line option enables calendar-based time management; see
:numref:`tab-calendars`. The implications of selecting a calendar are:

- PISM uses the ``units`` attribute of coordinate variables *literally* (including the
  reference date) in unit conversions. Please make sure that the :var:`time` variable in
  all forcing files has the units attribute such as "``days since 2012-1-1``". PISM will
  stop with an error message if a time variable does not have a reference date in its unit
  specification.

- It is important to use units that are a fixed multiple of "seconds", such as "``minutes
  since 1989-1-1``" or "``days since 1999-12-31``" and avoid "months" and "years". (PISM
  uses UDUNITS-2 to convert units, and in UDUNITS one month is always interpreted as
  `\frac{1}{12}\cdot 365.242198781` days.) Please see the `CF Conventions`_ document
  for details.

- PISM uses dates in standard output:

  .. code-block:: none

     ...
        time interval (length)   [2012-01-01, 2021-12-31]  (10.000 years)
     ...
     S 2012-05-26:  0.00011    0.6306   0.00000000           0.00000
     $v$Eh m (dt=0.10000)
     S 2012-07-01:  0.00014    0.6306   0.00000000           0.00000

  Just like in the no-calendar mode, run length, run start and run end times are specified
  using :opt:`-y`, :opt:`-ys` and :opt:`-ye` command-line options, respectively. Arguments
  of these options are interpreted in a slightly different manner, though:

- the run length option ``-y`` takes an *integer* argument, interpreted as the number of
  *calendar* years

- options ``-ys`` and ``-ye`` take *dates* as arguments.

For example, either of the following commands sets up a run covering the 21`^{st}`
century:

.. code-block:: none

   pismr -calendar gregorian -ys 2001-1-1 -y 100 ...
   pismr -calendar standard -ys 2001-1-1 -ye 2101-1-1 ...

(These option combinations are equivalent.)

It is also possible to run PISM for the duration of the available forcing data using the
:opt:`-time_file` option. The command

.. code-block:: none

   pismr -calendar gregorian -time_file forcing.nc

will extract the reference date and run length from ``forcing.nc``, respecting time
bounds.

When a non-trivial calendar is selected, spatial and scalar time-series can be saved
daily, monthly or yearly using these calendric computations. See sections
:ref:`sec-saving-time-series` and :ref:`sec-saving-spat-vari`.

.. _sec-time-file-restart:

Re-starting an interrupted run using ``-time_file``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a run using ``-time_file`` gets interrupted but manages to save a backup, re-starting
with ``-time_file`` will attempt to re-do the entire run because options ``-y``, ``-ys``,
and ``-ye`` are ignored:

.. code:: bash

   # This run gets killed but leaves backup.nc:
   pismr -i input.nc -time_file time.nc -o output.nc
   # This WILL NOT start from the time saved in backup.nc
   # and continue until the end time in time.nc
   pismr -i backup.nc -time_file time.nc -o output.nc

In this case we want to set the start time of the run from ``backup.nc``, but use the end
time from ``time.nc``. To achieve this, use the option :opt:`-time_file_continue_run`.

.. code-block:: bash

   # This run gets killed but leaves backup.nc:
   pismr -i input.nc -time_file time.nc -o output.nc
   # This WILL continue until the end time in time.nc, starting from backup.nc
   pismr -i backup.nc -time_file time.nc -o output.nc -time_file_continue_run

.. _sec-diagnostic-computations:

Diagnostic computations
-----------------------

A "diagnostic" computation can be defined as one where the internal state does not evolve.
The internal state of PISM is the set of variables read by "``-i``". You can ask PISM to
do a diagnostic computation by setting the run duration to a small number such as
`0.001` years (about `9` hours). The duration to use depends on the modeling
setup, but should be smaller than the maximum time-step allowed by PISM's stability
criteria. Such short runs can also be used to look at additional fields corresponding to
the current model state.

As an example, consider these two runs:

.. code-block:: none

   pisms -y 6000 -o foo.nc
   pismr -i foo.nc -y 0.001 -o bar.nc -o_size big

The result of the second (short) run is a NetCDF file ``bar.nc`` which contains the full
three-dimensional velocity field in the scalar NetCDF variables ``uvel``, ``vvel``, and
``wvel``, as well as many other variables. The file ``foo.nc`` does not contain many of
these fields because it was written with the default output size of ``medium``. The "``-y
0.001``" run has diagnostically "filled-in" all the fields which PISM can model at a time
step, but the run duration was chosen so as to avoid significant model state evolution
during the run.

This diagnostic mode is often associated to the modeling of ice shelves and ice streams.
Subsection :ref:`sec-ross` describes using a short "diagnostic" run to model the Ross ice
shelf [MacAyealetal]_. Verification tests I and J, section :ref:`sec-verif`, are
diagnostic calculations using the SSA.

The NetCDF model state saved by PISM at the end of an *evolution* run (i.e. with "``-y
Y``" for `Y>0`) does not, under the default ``-o_size medium`` output size, contain
the three-dimensional velocity field. Instead, it contains just a few more variables than
those which are needed to restart the run with ``-i``. One can force PISM to save all the
supported diagnostic quantities at the end of a time-stepping run using the option
``-o_size big``. Or one can go back and do a "``-y small_number``" diagnostic run using
``-o_size big``.

.. _sec-turning-off:

Disabling PISM components
-------------------------


Certain major model components, unlike more peripheral ones like bed deformation or
calving, are "on" by default. They do not need to be turned on explicitly. For example,
the SIA computation is so common that it would be a hassle to require an option to turn it
on every time you need it.

But sometimes one wants to disable particular components, during model spin-up, for
example. PISM has the following "off" switches:

- :opt:`-no_mass` disables the mass-continuity (conservation of mass) step
- :opt:`-energy none` disables the conservation of energy computation
- :opt:`-energy cold` makes PISM use temperature instead of enthalpy in the energy
  conservation code
- :opt:`-stress_balance none` disables the stress balance computation (useful for testing
  surface mass balance inputs)

.. _sec-hard-choices:

Dealing with more difficult modeling choices
--------------------------------------------

Most uses of an ice sheet model depend on careful modeling choices in situations where
there are considerable uncertainties *and* the model results depend strongly on those
choices. There may be, at the present state of knowledge, *no clear default values* that
PISM can provide. Furthermore, the available PISM options and sub-models are known to
*not* be sufficient for all users. Thus there are modelling situations for which we know
the user may have to do a great deal more hard work than just choose among PISM runtime
options.

Here are example cases where users have worked hard:

- User made use of available data in order to choose parameters for existing PISM models.
  These parameters then override PISM defaults.

  .. admonition:: Example
     :class: note

     Use regional atmosphere model output to identify PDD parameters suitable for modeling
     surface mass balance on a particular ice sheet. Then supply these parameters to PISM
     by a ``-config_override`` file.

     .. our UAF current situation with Greenland

- User wrote code, including code which modified current PISM internals, either to add
  additional processes or to "correct" PISM default process models.

  .. admonition:: Example
     :class: note

     Add a new sub-ice-shelf melt model by modifying C++ code in the ``src/coupler/``
     directory.

     .. PIK ocean models

- User simplified the model in use, instead of the default which was more elaborate.

  .. admonition:: Example
     :class: note

     Instead of using the PISM default mechanism connecting basal melt rate and basal
     strength, bypass this mechanism by generating a map of yield stress ``tauc`` directly
     and supplying it as input.

     .. Nick's -yield_stress constant choice

.. rubric:: Footnotes

.. [#] In most cases one process corresponds to one "core" of your computer.

.. [#] Note seconds are part of SI units.

