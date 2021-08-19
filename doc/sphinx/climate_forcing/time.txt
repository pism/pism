.. include:: shortcuts.txt

.. _sec-model-time:

Managing model time
-------------------

Most of PISM only needs to know how long the current time step is. The climate forcing
(reporting) code, on the other hand, uses time in a precise manner to provide (and report)
the correct values at the right time. For example: the February mass balance should be
used for 28 days (except during leap years) and not `365/12 = 30.4167` days.

.. _sec-periodic-forcing:

Periodic forcing
++++++++++++++++

All components reading time-dependent forcing data from files can interpret it as
"periodic". The length of the period (in years) is specified using a :opt:`-..._period`
option. For example, to prescribe a periodic climate which has the same values each year
but which includes inter-annual variations, using the :opt:`-surface given` option, set:

.. code-block:: none

    -surface given -surface_given_period 1 -surface_given_file forcing.nc

Each component has a unique command-line option prefix for a :opt:`-..._period` option.
Please refer to corresponding sections for allowed prefixes.

If forcing data has the period other than one year it is also necessary to specify the
"starting time" using the :opt:`-..._reference_year` option.

For example, to use a 20 year long climate record as periodic climate starting at the
beginning of the model year 10, do

.. code-block:: none

    -surface given -surface_given_period 20 -surface_given_file forcing.nc \
    -surface_given_reference_year 10

Note that the reference year is given in *model years*, not calendar years.

The :var:`time` variable in a forcing file that is to be used as periodic should start at
`0`. (In other words, time in a file with periodic forcing data is *time since the
beginning of a period*.) Please see the :ref:`sec-time` for a discussion of time units
appropriate in forcing files.

.. _sec-time-bounds:

Using time bounds in forcing data
+++++++++++++++++++++++++++++++++

To make :ref:`balancing the books <sec-mass-conservation>` possible, PISM interprets *fluxes*
(examples: top surface mass balance, precipitation, sub shelf mass flux) as
piecewise-constant in time. A forcing file is required to contain time bounds
corresponding to each record.

Other 2D input fields (examples: ice surface temperature, near-surface air temperature)
are interpreted as piecewise-linear in time.

*Scalar* time-dependent inputs are interpreted as piecewise-constant in time if an input
file contains time bounds and piecewise-linear otherwise.

PISM supports time bounds specified according to `CF Conventions`_. The ``ncdump -h``
output from a conforming file would look similar to:

.. code-block:: none

    netcdf forcing {
    dimensions:
            time = UNLIMITED ; // (214 currently)
            nv = 2 ;
    variables:
            double time(time) ;
                    time:units = "seconds since 2000-1-1" ;
                    time:axis = "T" ;
                    time:bounds = "time_bounds" ;
                    time:calendar = "gregorian" ;
                    time:long_name = "time" ;
            double nv(nv) ;
            double time_bounds(time, nv) ;

The :var:`time_bounds` variable stores the starting and the ending time for each interval
in the forcing. This variable is assumed to have the same units as the :var:`time`
variable it is associated with, which is why these attributes are not set in this example.
