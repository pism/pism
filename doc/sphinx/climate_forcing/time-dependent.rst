.. include:: shortcuts.txt

.. _sec-forcing-time-dependent:

Using time-dependent forcing
----------------------------

This section describes the way PISM interprets time-dependent forcing.

All time-dependent forcing

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

.. _sec-periodic-forcing:

Periodic forcing
++++++++++++++++
