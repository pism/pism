.. include:: shortcuts.txt

Atmosphere model components
---------------------------

.. contents::

.. _sec-atmosphere-given:

Reading boundary conditions from a file
+++++++++++++++++++++++++++++++++++++++

:|options|: ``-atmosphere given``
:|variables|: :var:`air_temp`, :var:`precipitation` |flux|
:|implementation|: ``pism::atmosphere::Given``
:|seealso|: :ref:`sec-surface-given`

.. note:: This is the default choice.

Command-line options:

- :opt:`-atmosphere_given_file` prescribes an input file
- :opt:`-atmosphere_given_period` (*years*) makes PISM interpret data in
  ``-atmosphere_given_file`` as periodic. See section :ref:`sec-periodic-forcing`.
- :opt:`-atmosphere_given_reference_year` sets the reference model year; see section
  :ref:`sec-periodic-forcing`.

A file ``foo.nc`` used with ``-atmosphere given -atmosphere_given_file foo.nc`` should
contain several records; the :var:`time` variable should describe what model time these
records correspond to.

This model component was created to force PISM with sampled (possibly periodic) climate
data, e.g. using monthly records of :var:`air_temp` and :var:`precipitation`.

It can also used to drive a temperature-index (PDD) climatic mass balance computation
(section :ref:`sec-surface-pdd`).

.. _sec-atmosphere-yearly-cycle:

Cosine yearly cycle
+++++++++++++++++++

:|options|: :opt:`-atmosphere yearly_cycle`
:|variables|: :var:`air_temp_mean_annual`, 
              :var:`air_temp_mean_july`,
              :var:`precipitation` |flux|
              :var:`amplitude_scaling`
:|implementation|: ``pism::atmosphere::CosineYearlyCycle``

This atmosphere model component computes the near-surface air temperature using the
following formula:

.. math::

   T(\mathrm{time}) = T_{\text{mean annual}}
   + A(\mathrm{time})\cdot(T_{\text{mean July}} - T_{\text{mean annual}}) \cdot \cos(2\pi t),

where `t` is the year fraction "since last July"; the summer peak of the cycle is on
:config:`atmosphere.fausto_air_temp.summer_peak_day`, which is set to day `196` by
default (approximately July 15).

Here `T_{\text{mean annual}}` (variable :var:`air_temp_mean_annual`) and
`T_{\text{mean July}}` (variable :var:`air_temp_mean_july`) are read from a file
selected using the :opt:`-atmosphere_yearly_cycle_file` command-line option. A
time-independent precipitation field (variable :var:`precipitation`) is read from the same
file.

Optionally a time-dependent scalar amplitude scaling `A(t)` can be used. Specify a
file to read it from using the :opt:`-atmosphere_yearly_cycle_scaling_file` command-line
option. Without this option `A(\mathrm{time}) \equiv 1`.

.. _sec-atmosphere-searise-greenland:

SeaRISE-Greenland
+++++++++++++++++
    
:|options|: ``-atmosphere searise_greenland``
:|variables|: :var:`lon`,
              :var:`lat`,
              :var:`precipitation` |flux|
:|implementation|: ``pism::atmosphere::SeaRISEGreenland``
:|seealso|: :ref:`sec-atmosphere-paleo-precip`

This atmosphere model component implements a longitude, latitude, and elevation dependent
near-surface air temperature parameterization and a cosine yearly cycle described in
:cite:`Faustoetal2009` and uses a constant in time ice-equivalent precipitation field (in units
of thickness per time, variable :var:`precipitation`) that is read from an input (``-i``)
file. To read time-independent precipitation from a different file, use the option
:opt:`-atmosphere_searise_greenland_file`.

The air temperature parameterization is controlled by configuration parameters with the
``snow_temp_fausto`` prefix.

See also the ``-atmosphere ...,paleo_precip`` modifier, section
:ref:`sec-atmosphere-paleo-precip`, for an implementation of the SeaRISE-Greenland formula
for paleo-precipitation correction from present; a 7.3\% change of precipitation rate for
every one degree Celsius of temperature change :cite:`Huybrechts02`.

.. _sec-atmosphere-pik:

PIK
+++
    
:|options|: :opt:`-atmosphere pik`
:|variables|: :var:`lat`,
              :var:`precipitation`
:|implementation|: ``pism::atmosphere::PIK``

This model component reads a time-independent precipitation field from an input
(:opt:`-i`) file and computes near-surface air temperature using a latitude and surface
elevation-dependent formula.

The parameterization is the same as in the :opt:`-surface pik` model, section
:ref:`sec-surface-pik`.

.. _sec-atmosphere-one-station:

One weather station
+++++++++++++++++++

:|options|: :opt:`-atmosphere one_station`
            :opt:`-atmosphere_one_station_file`
:|variables|: :var:`air_temp` [Kelvin],
              :var:`precipitation` |flux|
:|implementation|: ``pism::atmosphere::WeatherStation``

This model component reads scalar time-series of the near-surface air temperature and
precipitation from a file specified using the :opt:`-atmosphere_one_station_file` option
and uses them at *all* grid points in the domain. In other words, resulting climate fields
are constant in space but not necessarily in time.

The :opt:`-atmosphere one_station` model should be used with a modifier such as
``lapse_rate`` (see section :ref:`sec-atmosphere-lapse-rate`) to create spatial
variablitity.

.. _sec-atmosphere-delta-t:

Scalar temperature offsets
++++++++++++++++++++++++++

:|options|: ``-atmosphere ...,delta_T``
:|variables|: :var:`delta_T`
:|implementation|: ``pism::atmosphere::Delta_T``

This modifier applies scalar time-dependent air temperature offsets to the output of an
atmosphere model. It takes the following command-line options.

- :opt:`-atmosphere_delta_T_file` sets the name of the file PISM will read :var:`delta_T`
  from.
- :opt:`-atmosphere_delta_T_period` (*years*) sets the period of the forcing data (section
  :ref:`sec-periodic-forcing`).
- :opt:`-atmosphere_delta_T_reference_year` sets the reference year (section
  :ref:`sec-periodic-forcing`).

Please make sure that :var:`delta_T` has the units of "``Kelvin``".

.. _sec-atmosphere-delta-p:

Scalar precipitation offsets
++++++++++++++++++++++++++++

:|options|: :opt:`-atmosphere ...,delta_P`
:|variables|: :var:`delta_P` |flux|
:|implementation|: ``pism::atmosphere::Delta_P``

This modifier applies scalar time-dependent precipitation offsets to the output of an
atmosphere model. It takes the following command-line options.

- :opt:`-atmosphere_delta_P_file` sets the name of the file PISM will read :var:`delta_P`
  from.
- :opt:`-atmosphere_delta_P_period` (*years*) sets the period of the forcing data (section
  :ref:`sec-periodic-forcing`).
- :opt:`-atmosphere_delta_P_reference_year` sets the reference year (section
  :ref:`sec-periodic-forcing`).

.. _sec-atmosphere-frac-p:

Scalar precipitation scaling
++++++++++++++++++++++++++++

:|options|: ``-atmosphere ...,frac_P``
:|variables|: :var:`frac_P` [no unit]
:|implementation|: ``pism::atmosphere::Frac_P``

This modifier scales precipitation output of an atmosphere model using a scalar
time-dependent precipitation fraction, with a value of one corresponding to no change in
precipitation. It takes the following command-line options:

- :opt:`-atmosphere_frac_P_file` sets the name of the file PISM will read :var:`frac_P`
  from.
- :opt:`-atmosphere_frac_P_period` (*years*) sets the period of the forcing data (section
  :ref:`sec-periodic-forcing`).
- :opt:`-atmosphere_frac_P_reference_year` sets the reference year (section
  :ref:`sec-periodic-forcing`).

.. _sec-atmosphere-paleo-precip:

Precipitation correction using scalar temperature offsets
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-atmosphere ...,paleo_precip``
:|variables|: :var:`delta_T` [degrees Kelvin]
:|implementation|: ``pism::atmosphere::PaleoPrecip``

This modifier implements the SeaRISE-Greenland formula for a precipitation correction from
present; a 7.3\% change of precipitation rate for every one degree Celsius of air
temperature change :cite:`Huybrechts02`. See `SeaRISE Greenland model initialization
<SeaRISE-Greenland_>`_ for details. The input file should contain air temperature offsets
in the format used by ``-atmosphere ...,delta_T`` modifier, see section :ref:`sec-atmosphere-delta-t`.

It takes the following command-line options.

- :opt:`-atmosphere_paleo_precip_file` sets the name of the file PISM will read
  :var:`delta_T` from.
- :opt:`-atmosphere_paleo_precip_period` (*years*) sets the period of the forcing data
  (section :ref:`sec-periodic-forcing`).
- :opt:`-atmosphere_paleo_precip_reference_year` sets the reference year (section
  :ref:`sec-periodic-forcing`).

.. _sec-atmosphere-lapse-rate:

Temperature and precipitation lapse rate corrections
++++++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: :opt:`-atmosphere ...,lapse_rate`
:|variables|: :var:`surface_altitude` (CF standard name)
:|implementation|: ``pism::atmosphere::LapseRates``

The ``lapse_rate`` modifier allows for correcting air temperature and precipitation using
elevation lapse rates. It uses the following options.

- :opt:`-temp_lapse_rate` gives the temperature lapse rate, in `K/km`. Note that we
  use the following definition of the temperature lapse rate:

  .. math::

    \gamma = -\frac{dT}{dz}.

- :opt:`-precip_lapse_rate` gives the precipitation lapse rate, in `(m/year)/km`.
  Here `\gamma = -\frac{dM}{dz}`.
- :opt:`-atmosphere_lapse_rate_file` specifies a file containing the reference surface
  elevation field (standard name: :var:`surface_altitude`). This file may contain several
  surface elevation records to use lapse rate corrections relative to a time-dependent
  surface. If one record is provided, the reference surface elevation is assumed to be
  time-independent.
- :opt:`-atmosphere_lapse_rate_period` gives the period, in model years; see section
  :ref:`sec-periodic-forcing`.
- :opt:`-atmosphere_lapse_rate_reference_year` specifies the reference date; see section
  :ref:`sec-periodic-forcing`.

.. _sec-atmosphere-anomaly:

Using climate data anomalies
++++++++++++++++++++++++++++

:|options|: :opt:`-atmosphere ...,anomaly`
:|variables|: :var:`air_temp_anomaly`,
              :var:`precipitation_anomaly` |flux|
:|implementation|: ``pism::atmosphere::Anomaly``

This modifier implements a spatially-variable version of ``-atmosphere
...,delta_T,delta_P``.

It takes the following options:

- :opt:`-atmosphere_anomaly_file` specifies a file containing variables
  :var:`air_temp_anomaly` and :var:`precipitation_anomaly`.
- :opt:`-atmosphere_anomaly_period` (years) specifies the period of the forcing data, in
  model years; section :ref:`sec-periodic-forcing`.
- :opt:`-atmosphere_anomaly_reference_year` specifies the reference year; section
  :ref:`sec-periodic-forcing`.

See also to ``-surface ...,anomaly`` (section :ref:`sec-surface-anomaly`), which is
similar, but applies anomalies at the surface level.
