.. include:: shortcuts.txt

.. _sec-forcing-examples:

Examples and corresponding options
----------------------------------

This section gives a very brief overview of some coupling options. Please see sections
referenced below for more information.

One way coupling to a climate model
+++++++++++++++++++++++++++++++++++

One-way coupling of PISM to a climate model can be achieved by reading a NetCDF file with
time- and space-dependent climate data produced by a climate model.

There are two cases:

- coupling to a climate model that includes surface (firn, snow) processes
- coupling to a climate model providing near-surface air temperature and precipitation

.. _sec-example-surface-given:

Reading ice surface temperature and mass balance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the simplest case. It is often the preferred case, for example when the climate
model in use has high quality surface mass and energy sub-models which are then preferred
to the highly simplified (e.g. temperature index) surface models in PISM.

:|variables|: :var:`climatic_mass_balance`, :var:`ice_surface_temp`
:|options|: :opt:`-surface given -surface_given_file forcing.nc`
:|seealso|: :ref:`sec-surface-given`

.. _sec-example-atmosphere-given:

Reading air temperature and precipitation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned above, if a climate model provides near-surface air temperature and
precipitation, these data need to be converted into top-of-the-ice temperature and
climatic mass balance.

One way to do that is by using a temperature index (PDD) model component included in PISM.
This component has adjustable parameters; default values come from :cite:`RitzEISMINT`.

:|variables|: :var:`precipitation`, :var:`air_temp`
:|options|: :opt:`-atmosphere given -atmosphere_given_file forcing.nc -surface pdd`
:|seealso|: :ref:`sec-atmosphere-given`, :ref:`sec-surface-pdd`

If melt is negligible :opt:`-surface pdd` should be replaced with :opt:`-surface simple`
(see section :ref:`sec-surface-simple`).

.. _sec-example-atmosphere-anomalies:

Using climate anomalies
+++++++++++++++++++++++

Prognostic modeling experiments frequently use time- and space-dependent air temperature
and precipitation anomalies.

:|variables|: :var:`precipitation`,
              :var:`air_temp`,
              :var:`precipitation_anomaly`,
              :var:`air_temp_anomaly`
:|options|: :opt:`-atmosphere given,anomaly`,
            :opt:`-atmosphere_given_file forcing.nc`,
            :opt:`-atmosphere_anomaly_file anomalies.nc`,
            :opt:`-surface simple`
:|seealso|: :ref:`sec-atmosphere-given`,
            :ref:`sec-atmosphere-anomaly`,
            :ref:`sec-surface-simple`

The ``simple`` surface model component re-interprets precipitation as climatic mass
balance, which is useful in cases when there is no melt (Antarctic simulations is an
example).

Simulations of the Greenland ice sheet typically use :opt:`-surface pdd` instead of
:opt:`-surface simple`.

.. _sec-example-searise-greenland:

SeaRISE-Greenland
+++++++++++++++++

The SeaRISE-Greenland setup uses a parameterized near-surface air temperature
:cite:`Faustoetal2009` and a constant-in-time precipitation field read from an input
(:opt:`-i`) file. A temperature-index (PDD) scheme is used to compute the climatic mass
balance.


:|variables|: :var:`precipitation`,
              :var:`lat`,
              :var:`lon`
:|options|:  :opt:`-atmosphere searise_greenland -surface pdd`
:|seealso|: :ref:`sec-atmosphere-searise-greenland`,
            :ref:`sec-surface-pdd`

The air temperature parameterization is a function of latitude (:var:`lat`), longitude
(:var:`lon`) and surface elevation (dynamically updated by PISM).

.. _sec-example-searise-greenland-paleo:

SeaRISE-Greenland paleo-climate run
+++++++++++++++++++++++++++++++++++

The air temperature parameterization in the previous section is appropriate for present
day modeling. PISM includes some mechanisms allowing for corrections taking into account
differences between present and past climates. In particular, one can use ice-core derived
scalar air temperature offsets :cite:`JohnsenetalGRIP`, precipitation adjustments
:cite:`Huybrechts02`, and sea level offsets from SPECMAP :cite:`Imbrieetal1984`.

:|variables|: :var:`precipitation`,
              :var:`delta_T`,
              :var:`delta_SL`,
              :var:`lat`,
              :var:`lon`
:|options|: :opt:`-atmosphere searise_greenland,delta_T -atmosphere_delta_T_file
            delta_T.nc -surface pdd -sea_level constant,delta_sl -ocean_delta_sl_file
            delta_SL.nc`
:|seealso|: :ref:`sec-atmosphere-searise-greenland`,
            :ref:`sec-atmosphere-delta-t`,
            :ref:`sec-surface-pdd`,
            :ref:`sec-ocean-constant`,
            :ref:`sec-ocean-delta-sl`
    
Note that the temperature offsets are applied to *air* temperatures at the *atmosphere
level*. This ensures that `\Delta T` influences the PDD computation.

.. _sec-example-antarctica-paleo:

Antarctic paleo-climate runs
++++++++++++++++++++++++++++

:|variables|: :var:`climatic_mass_balance`,
              :var:`air_temp`,
              :var:`delta_T`,
              :var:`delta_SL`
:|options|: :opt:`-surface given,delta_T -surface_delta_T_file delta_T.nc -sea_level
            constant,delta_sl -ocean_delta_sl_file delta_SL.nc`
:|seealso|: :ref:`sec-surface-given`,
            :ref:`sec-surface-delta-t`,
            :ref:`sec-ocean-constant`,
            :ref:`sec-ocean-delta-sl`
