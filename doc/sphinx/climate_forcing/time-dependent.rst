.. include:: shortcuts.txt

.. _sec-forcing-time-dependent:

Using time-dependent forcing
----------------------------

.. contents::

Introduction
++++++++++++

PISM can use time-dependent *scalar* and *spatially-variable* forcing inputs.

The ``ncdump`` output for a typical forcing file would look similar to :ref:`this example
<src-forcing-time-dependent>`.

.. literalinclude:: delta_T.cdl
   :language: none
   :caption: Example time-dependent forcing
   :name: src-forcing-time-dependent

A data set like this one could be used to model a scenario in which the temperature at
the top surface of the ice drops by 30 degrees over the course of 100 years, remains
constant for 800 years, then increases by 30 degrees over 100 years. See
:ref:`sec-surface-delta-t`.

In addition to *times*, all forcing files have to contain *time bounds* defining time
intervals corresponding to individual records.

.. note::

   PISM will check if the modeled time interval (set using :config:`time.start` and
   :config:`time.end` or :config:`time.run_length`) is a subset of the time covered by a
   forcing file and **stop** if it is not.

   Set :config:`input.forcing.time_extrapolation` to ``true`` to tell PISM to use constant
   extrapolation instead.

.. _sec-forcing-scalar:

Scalar time-dependent inputs
++++++++++++++++++++++++++++

All scalar time-dependent inputs are interpreted as *piecewise-linear*.

Time bounds are used to compute period length for periodic forcing and the time interval
covered by provided data otherwise. Only the left end point of the first interval and the
right end point of the last interval are used.


.. _sec-forcing-spatially-variable:

Spatially-variable time-dependent inputs
++++++++++++++++++++++++++++++++++++++++

To make :ref:`balancing the books <sec-mass-conservation>` possible PISM interprets
*fluxes* such as the top surface mass balance, precipitation, and the sub shelf mass flux
as *piecewise-constant in time* over intervals defined by time bounds.

In this case *times* corresponding to individual forcing records are irrelevant (and are
ignored) and only time bounds are used.

Other 2D input fields (examples: ice surface temperature, near-surface air temperature,
sea level elevation) are interpreted as *piecewise-linear in time*.

In this case times are read from the file and time bounds are used to compute period
length for periodic forcing and the time interval covered by provided data otherwise.

.. _sec-periodic-forcing:

Periodic forcing
++++++++++++++++

A PISM forcing file with a periodic forcing **has to contain exactly one period**.

The left end point of the first time interval defines the *start of the period*.

The total duration of forcing in the file defines the length of the period. Specifically,
the length of the period is the difference of the right end point of the last interval and
the left end point of the first interval.

When used as periodic forcing, :ref:`src-forcing-time-dependent` would be interpreted as
having the period of one thousand `365`\-day years, with the period starting on January
`1` of year `1`.

.. note::

   - A real life (Gregorian, Julian, etc) calendar does not usually make sense in
     simulations using periodic forcing.

   - It is usually a good idea to use time units that are an integer multiple of one
     second, for example "day" or "365 days" as in the example above. This makes forcing
     files easier to interpret. (The units "years" corresponds to the mean tropical year.
     This is appropriate when converting from ``m/s`` to ``m/year``, for example, but
     not for keeping track of time.)

.. _sec-adding-time-bounds:

Adding time bounds
++++++++++++++++++

NCO_\'s ``ncap2`` makes it easy to add time bounds to a data set.

Save one of the scripts below to ``add_time_bounds.txt`` and then run

.. code-block:: bash

   ncap2 -O -S add_time_bounds.txt forcing.nc forcing-with-bounds.nc

to add time bounds to ``forcing.nc``.

Times as mid-points of intervals
================================

Use this script to interpret each time as a mid-point of the corresponding interval.

.. code-block:: none

   defdim("nv",2);
   time_bnds=make_bounds(time,$nv,"time_bnds");

or just use this one command:

.. code-block:: bash
   :caption: Adding time bounds, interpreting times as mid-points of intervals

   ncap2 -O -s 'defdim("nv",2);time_bnds=make_bounds(time,$nv,"time_bnds");' \
         forcing.nc forcing-with-bounds.nc

Times as left end points of intervals
=====================================

Use this script to interpret each time as a the left end point of the corresponding
interval.

.. code-block:: none
   :caption: Adding time bounds, interpreting times as left end points of intervals

   time@bounds="time_bnds";
   defdim("nv",2);
   time_bnds=array(0,0,/$time,$nv/);
   time_bnds(:,0)=time;
   time_bnds(:-2,1)=time(1:);
   time_bnds(-1,1)=2*time(-1)-time(-2);

Times as right end points of intervals
======================================

Use this script to interpret each time as a the right end point of the corresponding
interval.

.. code-block:: none
   :caption: Adding time bounds, interpreting times as right end points of intervals

   time@bounds="time_bnds";
   defdim("nv",2);
   time_bnds=array(0,0,/$time,$nv/);
   time_bnds(:,1)=time;
   time_bnds(1:,0)=time(:-2);
   time_bnds(0,0)=2*time(0)-time(1);
