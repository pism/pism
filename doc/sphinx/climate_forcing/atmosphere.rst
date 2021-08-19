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

This model is controlled by the following parameters with the prefix ``atmosphere.given.``.

.. pism-parameters::
   :prefix: atmosphere.given.

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

Here `T_{\text{mean annual}}` (variable :var:`air_temp_mean_annual`) and `T_{\text{mean
July}}` (variable :var:`air_temp_mean_july`) are read from a file selected using the
command-line option :opt:`-atmosphere_yearly_cycle_file`. A time-independent precipitation
field (variable :var:`precipitation`) is read from the same file.

Optionally a time-dependent scalar amplitude scaling `A(t)` can be used. Specify a
file to read it from using the :opt:`-atmosphere_yearly_cycle_scaling_file` command-line
option. Without this option `A(\mathrm{time}) \equiv 1`.

This model is controlled by the following parameters with the prefix
``atmosphere.yearly_cycle.``.

.. pism-parameters::
   :prefix: atmosphere.yearly_cycle.

.. _sec-atmosphere-searise-greenland:

SeaRISE-Greenland
+++++++++++++++++
    
:|options|: ``-atmosphere searise_greenland``
:|variables|: :var:`lon`,
              :var:`lat`,
              :var:`precipitation` |flux|
:|implementation|: ``pism::atmosphere::SeaRISEGreenland``
:|seealso|: :ref:`sec-atmosphere-precip-scaling`

This atmosphere model component implements a longitude, latitude, and elevation dependent
near-surface air temperature parameterization and a cosine yearly cycle described in
:cite:`Faustoetal2009` and uses a constant in time ice-equivalent precipitation field (in units
of thickness per time, variable :var:`precipitation`) that is read from an input (``-i``)
file. To read time-independent precipitation from a different file, use the option
:opt:`-atmosphere_searise_greenland_file`.

The air temperature parameterization is controlled by configuration parameters with the
prefix ``atmosphere.fausto_air_temp``:

.. pism-parameters::
   :prefix: atmosphere.fausto_air_temp.

See :ref:`sec-atmosphere-precip-scaling` for an implementation of the SeaRISE-Greenland
formula for precipitation adjustment using air temperature offsets relative to present; a
7.3\% change of precipitation rate for every one degree Celsius of temperature change
:cite:`Huybrechts02`.

.. _sec-atmosphere-pik:

PIK
+++

:|options|: :opt:`-atmosphere pik`
:|variables|: :var:`lat`, :var:`lon`, :var:`precipitation`
:|implementation|: ``pism::atmosphere::PIK``

This model component reads a time-independent precipitation field from an input file
specified by :config:`atmosphere.pik.file` and computes near-surface air temperature using
a parameterization selected using :config:`atmosphere.pik.parameterization`.

.. note::

   * Parameterizations implemented in this model are appropriate for the **Antarctic** ice
     sheet.

   * All parameterizations except for the first one implement a cosine annual cycle of the
     near-surface air temperature.

.. list-table:: Near-surface air temperature parameterizations
   :header-rows: 1
   :widths: 1,2

   * - Keyword
     - Description

   * - ``martin`` (default)
     - Uses equation (1) from :cite:`Martinetal2011` to parameterize mean annual
       temperature with *no seasonal variation in temperature.*

   * - ``huybrechts_dewolde``
     - Mean annual and mean summer temperatures are computed using parameterizations for
       the Antarctic ice sheet described in :cite:`HuybrechtsdeWolde` (equations C1 and C2).

   * - ``martin_huybrechts_dewolde``
     - Hybrid of the two above: mean annual temperature as in :cite:`Martinetal2011` with
       the amplitude of the annual cycle from :cite:`HuybrechtsdeWolde`.

   * - ``era_interim``
     - Mean annual and mean summer temperatures use parameterizations based on multiple
       regression analysis of ERA INTERIM data.

   * - ``era_interim_sin``
     - Mean annual and mean summer temperatures use parameterizations based on multiple
       regression analysis of ERA INTERIM data with a `\sin(\text{latitude})` dependence.

   * - ``era_interim_lon``
     - Mean annual and mean summer temperatures use parameterizations based on multiple
       regression analysis of ERA INTERIM data with a `\cos(\text{longitude})` dependence.

See :ref:`sec-surface-pik` for a surface model that implements the ``martin`` choice as a
parameterization of the ice temperature at its top surface.

.. _sec-atmosphere-one-station:

One weather station
+++++++++++++++++++

:|options|: :opt:`-atmosphere one_station`
            :opt:`-atmosphere_one_station_file`
:|variables|: :var:`air_temp` [Kelvin],
              :var:`precipitation` |flux|
:|implementation|: ``pism::atmosphere::WeatherStation``

This model component reads scalar time-series of the near-surface air temperature and
precipitation from a file specified using :config:`atmosphere.one_station.file`
and uses them at *all* grid points in the domain. In other words, resulting climate fields
are constant in space but not necessarily in time.

The :opt:`-atmosphere one_station` model should be used with a modifier such as
``elevation_change`` (see section :ref:`sec-atmosphere-elevation-change`) to create spatial
variablitity.

.. _sec-atmosphere-delta-t:

Scalar temperature offsets
++++++++++++++++++++++++++

:|options|: ``-atmosphere ...,delta_T``
:|variables|: :var:`delta_T`
:|implementation|: ``pism::atmosphere::Delta_T``

This modifier applies scalar time-dependent air temperature offsets to the output of an
atmosphere model. It uses the following parameters:

.. pism-parameters::
   :prefix: atmosphere.delta_T.

Please make sure that :var:`delta_T` has the units of "``Kelvin``".

.. _sec-atmosphere-delta-p:

Scalar precipitation offsets
++++++++++++++++++++++++++++

:|options|: :opt:`-atmosphere ...,delta_P`
:|variables|: :var:`delta_P` |flux|
:|implementation|: ``pism::atmosphere::Delta_P``

This modifier applies scalar time-dependent precipitation offsets to the output of an
atmosphere model. It uses the following parameters:

.. pism-parameters::
   :prefix: atmosphere.delta_P.

.. _sec-atmosphere-frac-p:

Precipitation scaling
+++++++++++++++++++++

:|options|: ``-atmosphere ...,frac_P``
:|variables|: :var:`frac_P` [1]
:|implementation|: ``pism::atmosphere::Frac_P``

This modifier scales precipitation output of an atmosphere model using a time-dependent
precipitation fraction, with a value of *one* corresponding to no change in precipitation.
It supports both 1D (scalar) and 2D (spatially-variable) factors.

If the variable :var:`frac_P` in the input file (see :config:`atmosphere.frac_P.file`) depends
on time only it is used as a time-dependent constant-in-space scaling factor.

If the variable :var:`frac_P` has more than one dimension PISM tries to use it as a
time-and-space-dependent scaling factor.

See configuration parameters with the prefix ``atmosphere.frac_P.``

.. pism-parameters::
   :prefix: atmosphere.frac_P.

.. _sec-atmosphere-precip-scaling:

Precipitation correction using scalar temperature offsets
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-atmosphere ...,precip_scaling``
:|variables|: :var:`delta_T` [degrees Kelvin]
:|implementation|: ``pism::atmosphere::PrecipitationScaling``

This modifier implements the SeaRISE-Greenland formula for a precipitation correction from
present; a 7.3\% change of precipitation rate for every one degree Celsius of air
temperature change :cite:`Huybrechts02`. See `SeaRISE Greenland model initialization
<SeaRISE-Greenland_>`_ for details. The input file should contain air temperature offsets
in the format used by ``-atmosphere ...,delta_T`` modifier, see section :ref:`sec-atmosphere-delta-t`.

This mechanisms increases precipitation by `100(\exp(C) - 1)` percent for each degree of
temperature increase, where

`C =` :config:`atmosphere.precip_exponential_factor_for_temperature`.

It uses the following parameters:

.. pism-parameters::
   :prefix: atmosphere.precip_scaling.

.. _sec-atmosphere-elevation-change:

Adjustments using modeled change in surface elevation
+++++++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: :opt:`-atmosphere ...,elevation_change`
:|variables|: :var:`surface_altitude` (CF standard name)
:|implementation|: ``pism::atmosphere::ElevationChange``

The ``elevation_change`` modifier adjusts air temperature and precipitation using modeled
changes in surface elevation relative to a reference elevation read from a file.

The near-surface air temperature is modified using an elevation lapse rate `\gamma_T =`
:config:`atmosphere.elevation_change.temperature_lapse_rate`. Here

.. math::
   \gamma_T = -\frac{dT}{dz}.

Two methods of adjusting precipitation are available:

- Scaling using an exponential factor

  .. math::

     \mathrm{P} = \mathrm{P_{input}} \cdot \exp(C \cdot \Delta T),

  where `C =` :config:`atmosphere.precip_exponential_factor_for_temperature` and `\Delta
  T` is the temperature difference produced by applying
  :config:`atmosphere.elevation_change.temperature_lapse_rate`.

  This mechanisms increases precipitation by `100(\exp(C) - 1)` percent for each degree of
  temperature increase.

  To use this method, set :opt:`-precip_adjustment scale`.

- Elevation lapse rate for precipitation

  .. math::

     \mathrm{P} = \mathrm{P_{input}} - \Delta h \cdot \gamma_P,

  where `\gamma_P =` :config:`atmosphere.elevation_change.precipitation.lapse_rate` and
  `\Delta h` is the difference between modeled and reference surface elevations.

  To use this method, set :opt:`-smb_adjustment shift`.

It uses the following parameters:

.. pism-parameters::
   :prefix: atmosphere.elevation_change.

The file specified using :config:`atmosphere.elevation_change.file` may contain several
surface elevation records to use lapse rate corrections relative to a time-dependent
surface. If one record is provided, the reference surface elevation is assumed to be
time-independent.

.. _sec-atmosphere-anomaly:

Using climate data anomalies
++++++++++++++++++++++++++++

:|options|: :opt:`-atmosphere ...,anomaly`
:|variables|: :var:`air_temp_anomaly`,
              :var:`precipitation_anomaly` |flux|
:|implementation|: ``pism::atmosphere::Anomaly``

This modifier implements a spatially-variable version of ``-atmosphere
...,delta_T,delta_P``.

It uses the following parameters:

.. pism-parameters::
   :prefix: atmosphere.anomaly.

See also to ``-surface ...,anomaly`` (section :ref:`sec-surface-anomaly`), which is
similar but applies anomalies at the surface level.

.. _sec-orographic-precipitation:

Orographic precipitation
++++++++++++++++++++++++

:|options|: :opt:`-atmosphere ...,orographic_precipitation`
:|variables|: None
:|implementation|: ``pism::atmosphere::OrographicPrecipitation``

This modifier implements the linear orographic precipitation model described in
:cite:`SmithBarstad2004` with a modification incorporating the influence of the Coriolis
force from :cite:`SmithBarstadBonneau2005`.

We compute the Fourier transform of the precipitation field using the formula below (see
equation 49 in :cite:`SmithBarstad2004` or equation 3 in :cite:`SmithBarstadBonneau2005`).

.. math::
   :label: eq-orographic-precipitation

   \hat P_{\text{LT}}(k, l) = \frac{C_w i \sigma \hat h(k, l)}
   {(1 - i m H_w) (1 + i \sigma \tau_c) (1 + i \sigma \tau_f)},

where `h` is the surface elevation, `C_w = \rho_{S_{\text{ref}}} \Gamma_m / \gamma`
relates the condensation rate to vertical motion (see the appendix of
:cite:`SmithBarstad2004`), `m` is the vertical wavenumber (see equation 6 in
:cite:`SmithBarstadBonneau2005`), and `\sigma` is the intrinsic frequency. The rest of the
constants are defined below.

The spatial pattern of precipitation is recovered using an inverse Fourier transform
followed by post-processing:

.. math::
   :label: eq-orographic-post-processing

   P = \max(P_{\text{pre}} + P_{\text{LT}}, 0) \cdot S + P_{\text{post}}.

.. note::

   - Discontinuities in the surface gradient (e.g. at ice margins) may cause oscillations
     in the computed precipitation field (probably due to the Gibbs phenomenon). To address
     this our implementation includes the ability to smooth the surface topography using a
     Gaussian filter. Set
     :config:`atmosphere.orographic_precipitation.smoothing_standard_deviation` to a
     positive number to enable smoothing. Values around `\dx` appear to be effective.

   - The spectral method used to implement this model requires that the input (i.e.
     surface elevation `h`) is periodic in `x` and `y`. To simulate periodic `h` the
     implementation uses an extended grid (see
     :config:`atmosphere.orographic_precipitation.grid_size_factor`) and pads modeled
     surface elevation with zeros.

     It is worth noting that the resulting precipitation field `P_{\text{LT}}` *is also
     periodic* in `x` and `y` on the *extended* grid. The appropriate size of this
     extended grid may differ depending on the spatial scale of the domain and values of
     model parameters.

It is implemented as a "modifier" that overrides the precipitation field provided by an
input model. Use it with a model providing air temperatures to get a complete model. For
example, ``-atmosphere yearly_cycle,orographic_precipitation ...`` would use the annual
temperature cycle from ``yearly_cycle`` combined with precipitation computed using this
model.

The only spatially-variable input of this model is the surface elevation (`h` above)
modeled by PISM. It is controlled by a number of configuration parameters. See parameters
with the prefix ``atmosphere.orographic_precipitation.`` listed below.

.. pism-parameters::
   :prefix: atmosphere.orographic_precipitation.
