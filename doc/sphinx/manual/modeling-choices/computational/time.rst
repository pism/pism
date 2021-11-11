.. include:: ../../../global.txt

.. _sec-time:

Model time
----------

:numref:`tab-timeoptions` gives the command-line options which control PISM time. If
option ``-ys`` is absent then the start year is read from the input file (if present) or
it defaults to zero. The default value for the end of the run is the start time plus the
given (``-y``) run length. If both :opt:`-ys` and :opt:`-ye` are used then the run length
is set to the difference. If both :opt:`-y` and :opt:`-ye` are set PISM will use
:opt:`-ye` and ignore :opt:`-y`.

.. list-table:: Command-line options controlling PISM time
   :name: tab-timeoptions
   :header-rows: 1
   :widths: 20, 80

   * - Option
     - Description
   * - :opt:`-y` (years)
     - Number of model years to run (see :config:`time.run_length`).
   * - :opt:`-ys` (years)
     - Model year at which to start the run (see :config:`time.start`). Also resets the
       model time, ignoring any time in the input file.
   * - :opt:`-ye` (years)
     - Model year at which to end the run (see :config:`time.end`).

Times specified by :opt:`-ys` and :opt:`-ye` are interpreted as *years since the reference
date in the chosen calendar* (see :config:`time.reference_date` and
:config:`time.calendar`), so by default ``-ys 0`` will correspond to January 1 of year 1.

It is also possible to provide *units*: for example, ``-ys 1day -ye 2weeks`` will start the
run from January 2, 1 and run for 14 days (until January 15).

In addition to this, :opt:`-ys` and :opt:`-ye` (as well as :opt:`-extra_times`, etc) can
take *dates* of the form ``Y-M-D``, so ``-ys 2000-1-1 -ye 2100-1-1`` would tell PISM to
run for 100 years from January 1, 2000.

The run length given by :opt:`-y` is interpreted as the number of `365`\-day years. For
more precise control of the run duration, provide units (e.g. ``-y 1s`` for a very short
run or ``-y "1000 360day"`` to run for one thousand `360`\-day years) or use :opt:`-ye`
instead.

.. note::

   The fact that :opt:`-y` is not calendar-aware can lead to some confusion. Note, for
   example, that ``-ys 2000-1-1 -y 100 -calendar 360_day`` will end the run on
   ``2101-05-21``: this is the date in the `360`\-day calendar that is one hundred
   `365`\-day years from ``2000-1-1``.

It is also possible to run PISM for the duration of the available forcing data using the
:opt:`-time_file` option. The command

.. code-block:: none

   pismr -calendar gregorian -time_file forcing.nc

will extract the reference date and run length from ``forcing.nc``, respecting time
bounds.

.. note::

   * When preparing input files it is important to use time units that are a fixed multiple
     of "seconds", such as "``minutes since 1989-1-1``" or "``days since 1999-12-31``" and
     avoid "months" and "years". (PISM uses UDUNITS-2 to convert units, and in UDUNITS one
     month is always interpreted as `\frac{1}{12}\cdot 365.242198781` days.) Please see the
     `CF Conventions`_ document for details.

.. _sec-calendars:

Calendars
^^^^^^^^^

Most of PISM, and its ice dynamics core in particular, only needs to know the length of
the current time-step. Internally PISM stores time in "seconds since a specified moment"
and thus PISM generally does not use or need a calendar.\ [#f1]_ We refer to PISM internal
time as *model time*.

One can select a calendar for more precise control of the model time, however (see
:config:`time.calendar`). A "calendar" is a concept that is part of the `CF Conventions`_.
Choosing a calendar is appropriate for runs for specific temporal periods like "the
18th-century" or "1989--2010". The calendar is generally needed because specific knowledge
of lengths of months and years is required to use climate data properly or to facilitate
model validation.

PISM uses CalCalcs_ by David W. Pierce to perform calendric computations. This lets us
support all the `calendars <CF-Calendars_>`_ defined by the CF Metadata Conventions
document except for the ``366_day`` (``all_leap``) calendar.

By default PISM uses the ``365_day`` calendar. This is appropriate for runs that do not
require precise application of forcing data or reporting on particular dates
(paleo-climate runs, for example).

In addition to setting run start and end times, the calendar setting also affects the
interpretation of "monthly" and "yearly" reporting with :opt:`-extra_times`,
:opt:`-ts_times`, and :opt:`-save_times` (see :ref:`sec-saving-diagnostics`,
:ref:`sec-saving-time-series`, and :ref:`sec-snapshots`).

.. note::

   This does **not** affect unit conversion: the factor used to convert `m/s` to
   `m/\text{year}` does not depend on the calendar choice.

.. list-table:: Calendars supported by PISM. Please see CalCalcs_ documentation for
                details
   :name: tab-calendars
   :header-rows: 1
   :widths: 1, 3

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

.. _sec-time-file-restart:

Re-starting an interrupted run using ``-time_file``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a run using ``-time_file`` gets interrupted but manages to save a backup, re-starting
with ``-time_file`` will attempt to re-do the entire run because options ``-y``, ``-ys``,
and ``-ye`` are ignored:

.. code-block:: bash

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
^^^^^^^^^^^^^^^^^^^^^^^

A "diagnostic" computation can be defined as one where the internal state does not evolve.
The internal state of PISM is the set of variables read by "``-i``". You can ask PISM to
do a diagnostic computation by setting the run duration to a small number such as
`1` hours (``-y 10hours``). The duration to use depends on the modeling
setup, but should be smaller than the maximum time-step allowed by PISM's stability
criteria. Such short runs can also be used to look at additional fields corresponding to
the current model state.

As an example, consider these two runs:

.. code-block:: none

   pismr -eisII A -y 6000 -o foo.nc
   pismr -i foo.nc -y 10hours -o bar.nc -o_size big

The result of the second (short) run is a NetCDF file ``bar.nc`` which contains the full
three-dimensional velocity field in the scalar NetCDF variables ``uvel``, ``vvel``, and
``wvel``, as well as many other variables. The file ``foo.nc`` does not contain many of
these fields because it was written with the default output size of ``medium``. The "``-y
10hours``" run has diagnostically "filled-in" all the fields which PISM can model at a time
step, but the run duration was chosen so as to avoid significant model state evolution
during the run.

This diagnostic mode is often associated to the modeling of ice shelves and ice streams.
section :ref:`sec-ross` describes using a short "diagnostic" run to model the Ross ice
shelf :cite:`MacAyealetal`. Verification tests I and J, section :ref:`sec-verif`, are
diagnostic calculations using the SSA.

The NetCDF model state saved by PISM at the end of an *evolution* run (i.e. with "``-y
Y``" for `Y>0`) does not, under the default ``-o_size medium`` output size, contain
the three-dimensional velocity field. Instead, it contains just a few more variables than
those which are needed to restart the run with ``-i``. One can force PISM to save all the
supported diagnostic quantities at the end of a time-stepping run using the option
``-o_size big``. Or one can go back and do a "``-y small_number``" diagnostic run using
``-o_size big``.


.. rubric:: Footnotes

.. [#f1] Note seconds are part of SI units.
