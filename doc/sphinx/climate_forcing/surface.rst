.. include:: shortcuts.txt

Surface mass and energy process model components
------------------------------------------------

.. contents::

.. _sec-surface-simple:

The "invisible" model
+++++++++++++++++++++

:|options|: ``-surface simple``
:|variables|: none
:|implementation|: ``pism::surface::Simple``

This is the simplest "surface model" available in PISM, enabled using ``-surface simple``.
Its job is to re-interpret precipitation as climatic mass balance, and to re-interpret
mean annual near-surface (2m) air temperature as the temperature of the ice at the depth
at which firn processes cease to change the temperature of the ice. (I.e. the temperature
*below* the firn.) This implies that there is no melt. Though primitive, this model
component may be desired in cold environments (e.g. East Antarctic ice sheet) in which
melt is negligible and heat from firn processes is ignored.

.. _sec-surface-given:

Reading top-surface boundary conditions from a file
+++++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-surface given``
:|variables|: :var:`ice_surface_temp`, :var:`climatic_mass_balance` |flux|
:|implementation|: ``pism::surface::Given``

.. note::

   This is the default choice.

This model component was created to force PISM with sampled (possibly periodic) climate
data by reading ice upper surface boundary conditions from a file. These fields are
provided directly to the ice dynamics code (see :ref:`sec-climate-inputs` for details).

PISM will stop if variables :var:`ice_surface_temp` (ice temperature at the ice surface
but below firn) and :var:`climatic_mass_balance` (top surface mass flux into the ice) are
not present in the input file.

A file ``foo.nc`` used with ``-surface given -surface_given_file foo.nc`` may contain
several records. If this file contains one record (i.e. fields corresponding to one time
value only), provided forcing data is interpreted as time-independent. Variables
:var:`time` and :var:`time_bnds` should specify model times corresponding to individual
records.

For example, to use monthly periodic forcing with a period of 1 year starting at the
beginning of `1980` (let's use the `360`\-day calendar for simplicity), create a file
(say, "``foo.nc``") with 12 records. The :var:`time` variable may contain `15, 45, 75,
\dots, 345` (mid-month for all `12` months) and have the units of "``days since
1980-1-1``". (It is best to avoid units of "months" and "years" because their meanings
depend on the calendar.) Next, add the :var:`time_bounds` variable for the time dimension
with the values `0, 30, 30, 60, \dots` specifying times corresponding to beginnings and
ends of records for each month and set the ``time:bounds`` attribute accordingly. Now run

.. code-block:: none

    pismr -surface given -surface_given_file foo.nc -surface.given.periodic

See :ref:`sec-forcing-time-dependent` for more information.

.. note::

   - This surface model *ignores* the atmosphere model selection made using the option
     :opt:`-atmosphere`.
   - PISM can handle files with virtually any number of records: it will read and store in
     memory at most :config:`input.forcing.buffer_size` records at any given time
     (default: 60, or 5 years' worth of monthly fields).
   - when preparing a file for use with this model, it is best to use the ``t,y,x``
     variable storage order: files using this order can be read in faster than ones using
     the ``t,x,y`` order, for reasons :ref:`explained in the User's Manual
     <sec-pism-io-performance>`.
   
     To change the storage order in a NetCDF file, use ``ncpdq``:
   
     .. code-block:: none
   
       ncpdq -a t,y,x input.nc output.nc
   
     will copy data from ``input.nc`` into ``output.nc``, changing the storage order to
     ``t,y,x`` at the same time.

.. rubric:: Parameters

Prefix: ``surface.given.``

.. pism-parameters::
   :prefix: surface.given.


.. _sec-surface-elevation:

Elevation-dependent temperature and mass balance
++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-surface elevation``
:|variables|: none
:|implementation|: ``pism::surface::Elevation``

This surface model component parameterizes the ice surface temperature `T_{h}` =
:var:`ice_surface_temp` and the mass balance `m` = :var:`climatic_mass_balance` as
*piecewise-linear* functions of surface elevation `h`.

The option :opt:`-ice_surface_temp` (*list of 4 numbers*) determines the surface
temperature using the 4 parameters `\T{min}`, `\T{max}`, `\h{min}`,
`\h{max}`. Let

.. math::

  \diff{T}{h} = (\T{max} - \T{min}) / (\h{max} - \h{min})

be the temperature gradient. Then

.. math::

  T(x,y) =
  \begin{cases}
    \T{min}, & h(x,y) \le \h{min}, \\
    \T{min} + \diff{T}{h} \, (h(x,y) - \h{min}), & \h{min} < h(x,y) < \h{max}, \\
    \T{max}, & \h{max} \le h(x,y).
  \end{cases}

The option :opt:`-climatic_mass_balance` (*list of 5 numbers*) determines the surface mass
balance using the 5 parameters `\m{min}`, `\m{max}`, `\h{min}`,
`\h{ELA}`, `\h{max}`. Let

.. math::

   \diff{\m{abl}}{h} = -\m{min} / (\ELA - \h{min})

and

.. math::

   \diff{\m{acl}}{h} = \m{max} / (\h{max} - \ELA)

be the mass balance gradient in the ablation and in the accumulation area, respectively.
Then

.. math::

   m(x,y) =
   \begin{cases}
    \m{min}, & h(x,y) \le \h{min}, \\
    \diff{\m{abl}}{h} \, (h(x,y) - \ELA), & \h{min} < h(x,y) \le \ELA, \\
    \diff{\m{acl}}{h} \, (h(x,y) - \ELA), & \ELA < h(x,y) \le \h{max}, \\
    \m{max}, & \h{max} < h(x,y).
   \end{cases}

The option :opt:`-climatic_mass_balance_limits` (*list of 2 numbers*) limits the mass
balance below `\h{min}` to `\ms{min}` and above `\h{max}` to
`\ms{max}`, thus

.. math::

   m(x,y) =
   \begin{cases}
    m^{*}_{\text{min}}, & h(x,y) \le \h{min}, \\
    \diff{\m{abl}}{h} \, (h(x,y) - \ELA), & \h{min} < h(x,y) \le \ELA, \\
    \diff{\m{acl}}{h} \, (h(x,y) - \ELA), & \ELA < h(x,y) \le \h{max}, \\
    m^{*}_{\text{max}}, & \h{max} < h(x,y).
   \end{cases}

Note: this surface model *ignores* the atmosphere model selection made using the
:opt:`-atmosphere` option.

.. _sec-surface-pdd:

Temperature-index scheme
++++++++++++++++++++++++

:|options|: ``-surface pdd``
:|variables|: :var:`air_temp_sd`, :var:`snow_depth`
:|implementation|: ``pism::surface::TemperatureIndex``
                   
The default PDD model used by PISM, turned on by option :opt:`-surface pdd`, is based on
:cite:`CalovGreve05` and EISMINT-Greenland intercomparison (see :cite:`RitzEISMINT`).

Our model computes the solid (snow) precipitation rate using the air temperature threshold
with a linear transition. All precipitation during periods with air temperatures above
:config:`surface.pdd.air_temp_all_precip_as_rain` (default of `2^\circ C`) is interpreted as
rain; all precipitation during periods with air temperatures below
:config:`surface.pdd.air_temp_all_precip_as_snow` (default of `0^\circ C`) is interpreted as
snow.

For long-term simulations, a PDD model generally uses an idealized seasonal temperature
cycle. "White noise" is added to this cycle to simulate additional daily variability
associated to the vagaries of weather. This additional random variation is quite
significant, as the seasonal cycle may never reach the melting point but that point may be
reached with some probability, in the presence of the daily variability, and thus melt may
occur. Concretely, a normally-distributed, mean zero random temperature increment is added
to the seasonal cycle. There is no assumed spatial correlation of daily variability. The
standard deviation of the daily variability is controlled by configuration parameters with
the prefix ``surface.pdd.std_dev.``:

.. pism-parameters::
   :prefix: surface.pdd.std_dev.

A file ``foo.nc`` used with ``-surface pdd -pdd_sd_file foo.nc`` should contain standard
deviation of near-surface air temperature in variable :var:`air_temp_sd`, and the
corresponding time coordinate in variable :var:`time`. If ``-pdd_sd_file`` is not set,
PISM uses a constant value for standard deviation, which is set by the
configuration parameter :config:`surface.pdd.std_dev.value`. The default value is `5.0` degrees
:cite:`RitzEISMINT`. However, this approach is not recommended as it induces significant
errors in modeled surface mass balance in both ice-covered and ice-free regions
:cite:`RogozhinaRau2014`, :cite:`Seguinot2013`.

Over ice-covered grid cells, daily variability can also be parameterized as a linear
function of near-surface air temperature `\sigma = a \cdot T + b` using the
:config:`surface.pdd.std_dev.use_param` configuration flag, and the corresponding
parameters :config:`surface.pdd.std_dev.param_a` and
:config:`surface.pdd.std_dev.param_b`. This parametrization replaces prescribed standard
deviation values over glacierized grid cells as defined by the :var:`mask` variable (see
:config:`geometry.ice_free_thickness_standard`). Default values for the slope `a` and
intercept `b` were derived from the ERA-40 reanalysis over the Greenland ice sheet
:cite:`SeguinotRogozhina2014`.

The number of positive degree days is computed as the magnitude of the temperature
excursion above `0\!\phantom{|}^\circ \text{C}` multiplied by the duration (in days)
when it is above zero.

In PISM there are two methods for computing the number of positive degree days. The first
computes only the expected value, by the method described in :cite:`CalovGreve05`. This is
the default when a PDD is chosen (i.e. option :opt:`-surface pdd`). The second is a Monte
Carlo simulation of the white noise itself, chosen by adding the option :opt:`-pdd_method
random_process`. This Monte Carlo simulation adds the same daily variation at every point,
though the seasonal cycle is (generally) location dependent. If repeatable randomness is
desired use :opt:`-pdd_method repeatable_random_process` instead.

.. figure:: figures/pdd-model-flowchart.png
   :name: fig-pdd-model

   PISM's positive degree day model. `F_s` and `F_i` are PDD factors for snow
   and ice, respectively; `\theta_{\text{refreeze}}` is the refreeze fraction.

By default, the computation summarized in :numref:`fig-pdd-model` is performed every week.
(This frequency is controlled by the parameter :config:`surface.pdd.max_evals_per_year`.)
To compute mass balance during each week-long time-step, PISM keeps track of the current
snow depth (using units of ice-equivalent thickness). This is necessary to determine if
melt should be computed using the degree day factor for snow
(:config:`surface.pdd.factor_snow`) or the corresponding factor for ice
(:config:`surface.pdd.factor_ice`).

A fraction of the melt controlled by the configuration parameter :config:`surface.pdd.refreeze`
(`\theta_{\text{refreeze}}` in :numref:`fig-pdd-model`, default: `0.6`)
refreezes. The user can select whether melted ice should be allowed to refreeze using the
configuration flag :config:`surface.pdd.refreeze_ice_melt`.

Since PISM does not have a principled firn model, the snow depth is set to zero at the
beginning of the balance year. See :config:`surface.mass_balance_year_start_day`. Default is
`274`, corresponding to October 1\ :superscript:`st`.

Our PDD implementation is meant to be used with an atmosphere model implementing a cosine
yearly cycle such as ``searise_greenland`` (section
:ref:`sec-atmosphere-searise-greenland`), but it is not restricted to parameterizations
like these.

This code also implements latitude- and mean July temperature dependent ice and snow
factors using formulas (6) and (7) in :cite:`Faustoetal2009`; set :opt:`-pdd_fausto` to enable.
The default standard deviation of the daily variability (option :opt:`-pdd_std_dev`) is
2.53 degrees when :opt:`-pdd_fausto` is set :cite:`Faustoetal2009`. See also configuration
parameters with the prefix ``surface.pdd.fausto.``:

.. pism-parameters::
   :prefix: surface.pdd.fausto.

Note that when used with periodic climate data (air temperature and precipitation) that is
read from a file (see section :ref:`sec-atmosphere-given`), use of
:config:`time_stepping.hit_multiples` is recommended: set it to the length of the climate
data period in years.

This model provides the following scalar:

- :var:`surface_accumulation_rate`
- :var:`surface_melt_rate`
- :var:`surface_runoff_rate`

and these 2D diagnostic quantities (averaged over reporting intervals; positive flux
corresponds to ice gain):

- :var:`surface_accumulation_flux`
- :var:`surface_melt_flux`
- :var:`surface_runoff_flux`

This makes it easy to compare the surface mass balance computed by the model to its
individual components:

.. code::

   SMB = surface_accumulation_flux - surface_runoff_flux

.. rubric:: Parameters

Prefix: ``surface.pdd.``

.. pism-parameters::
   :prefix: surface.pdd.

.. _sec-surface-debm-simple:

Diurnal Energy Balance Model "dEBM-simple"
++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-surface debm_simple``
:|variables|: :var:`surface_albedo`
:|implementation|: ``pism::surface::DEBMSimple``

This PISM module implements the "simple" version of the diurnal energy balance model
developed by :cite:`Zeitzetal2021`. It follows :cite:`KrebsKanzow2018` and includes
parameterizations of the surface albedo and the atmospheric transmissivity that make it
possible to run the model in a standalone, prognostic mode.

It is designed to use time-dependent forcing by near-surface air temperature and total
(i.e. liquid *and* solid) precipitation provided by one of PISM's atmosphere models. The
temperature forcing should resolve the annual cycle, i.e. it should use *monthly or more
frequent* temperature records if forced using :ref:`sec-atmosphere-given`. In cases when
only annual temperature records are available we recommend using the
:ref:`sec-atmosphere-yearly-cycle` approximation.

.. note::

   We recommend setting :config:`time_stepping.hit_multiples` to the length of the climate
   data period in years when forcing dEBM-simple with periodic climate data that are read
   from a file.

Similarly to other surface models, the outputs are

- ice temperature at its top surface
- climatic mass balance (SMB)

dEBM-simple re-interprets near-surface air temperature to produce the ice temperature at
its top surface.

SMB is defined as

.. math::

   \text{SMB} &= \text{Solid accumulation} - \text{Runoff},

   \text{SMB} &= \text{Solid accumulation} - \left(\text{Melt} - \text{Refreeze} \right).

Solid accumulation is approximated using provided total precipitation and a linear
transition from interpreting it as "all snow" when the air temperature is below
:config:`surface.debm_simple.air_temp_all_precip_as_snow` (default of `0^\circ C`) to "all
rain" when the air temperature is above
:config:`surface.debm_simple.air_temp_all_precip_as_rain` (default of `2^\circ C`).
Alternatively, all precipitation is interpreted as snow if
:config:`surface.debm_simple.interpret_precip_as_snow` is set.

.. note::

   Part of the precipitation that is interpreted as rain is assumed to run off
   instantaneously and *does not* contribute to reported modeled runoff.

A fraction `\theta` (:config:`surface.debm_simple.refreeze`) of computed melt amount is
assumed to re-freeze:

.. math::
   :label: eq-debm-refreeze

   \text{Refreeze} = \theta\, \text{Melt}.

By default only snow melt is allowed to refreeze; set
:config:`surface.debm_simple.refreeze_ice_melt` to refreeze both snow and ice melt. To
distinguish between melted snow and ice dEBM-simple keeps track of the evolving snow depth
during a balance year and resets it to zero once a year on the day set using
:config:`surface.mass_balance_year_start_day`.

Let `T` be the near-surface air temperature and define

.. math::
   :label: eq-debm-melt-components

   M_I &= \frac{\Delta t_{\Phi}}{\Delta t \rho_{\text{w}} L_{\text{m}}} \tau_\text{A}
   ( 1 - \alpha ) \bar S_{\Phi},

   M_T &= \frac{\Delta t_{\Phi}}{\Delta t \rho_{\text{w}} L_{\text{m}}}
   c_1 T_\text{eff},

   M_O &= \frac{\Delta t_{\Phi}}{\Delta t \rho_{\text{w}} L_{\text{m}}}
   c_2,

then the average daily melt rate is approximated by

.. math::
   :label: eq-debm-melt

   M &=
   \begin{cases}
   M_I + M_T + M_O,&
   T \ge T_{\text{min}},\\
   0,& T < T_{\text{min}}.
   \end{cases},

Here

- `M_I` is the *insolation-driven melt contribution* representing the net uptake of
  incoming solar shortwave radiation of the surface during the diurnal melt period,
- `M_T` is the *temperature-driven melt contribution* representing the
  air-temperature-dependent part of the incoming longwave radiation as well as turbulent
  sensible heat fluxes,
- `M_O` is the negative *melt offset* representing the outgoing longwave radiation and the
  air temperature-independent part of the incoming longwave radiation :cite:`Garbe2023`.

See :numref:`tab-debm-simple-notation` for details and note that :eq:`eq-debm-melt` is
equation 1 in :cite:`Zeitzetal2021`.

The following two sections describe implementations of insolation-driven and
temperature-driven melt contributions.

.. list-table:: Notation used in :eq:`eq-debm-melt-components` and :eq:`eq-debm-melt`
   :header-rows: 1
   :widths: 2,4
   :name: tab-debm-simple-notation

   * - Quantity
     - Description

   * - `\Phi`
     - Threshold for the solar elevation angle (:config:`surface.debm_simple.phi`). It is
       assumed that melt can occur only when the sun is above this angle. `\Phi` should be
       treated as a tuning parameter since its value is not well constrained.

   * - `\Delta t_{\Phi} / \Delta t`
     - Fraction of the day during which the sun is above the elevation angle `\Phi` and
       melt can occur

   * - `\tau_A`
     - Transmissivity of the atmosphere

   * - `\alpha`
     - Surface albedo

   * - `\bar S_{\Phi}`
     - Mean top of the atmosphere insolation during the part of the day when the sun
       is above the elevation angle `\Phi`.

   * - `T_{\text{eff}}`
     - "Effective air temperature" computed using provided air temperature forcing and
       additional stochastic variations used to model the effect of daily temperature
       variations (see :eq:`eq-debm-t-eff`, :cite:`Zeitzetal2021` and
       :cite:`CalovGreve05`)

   * - `c_1`
     - Tuning parameter that controls the temperature influence on melt
       (:config:`surface.debm_simple.c1`)

   * - `c_2`
     - Tuning parameter that controls the (negative) melt offset
       (:config:`surface.debm_simple.c2`)

   * - `\rho_w`
     - Fresh water density (:config:`constants.fresh_water.density`)

   * - `L_m`
     - Latent heat of fusion (:config:`constants.fresh_water.latent_heat_of_fusion`)

   * - `T_{\text{min}}`
     - Threshold temperature (:config:`surface.debm_simple.melting_threshold_temp`). Melt
       is prohibited if the air temperature is below `T_{\text{min}}` to avoid melt rates
       from high insolation values and low albedo values when it is too cold to actually
       melt.

.. _sec-debm-insolation-driven-melt:

Insolation-driven melt contribution
===================================

.. math::
   :label: eq-debm-insolation-melt

   M_I = \frac{\Delta t_{\Phi}}{\Delta t \rho_{\text{w}} L_{\text{m}}} \left(
   \tau_\text{A} \left( 1 - \alpha \right) \bar S_{\Phi} \right),

This term models the influence of the *mean insolation during the melt period* `\bar
S_{\Phi}`, the atmosphere transmissivity `\tau_{A}` and the surface albedo `\alpha`.

.. _sec-debm-simple-insolation:

Mean top of the atmosphere insolation
#####################################

The mean top of the atmosphere insolation during the part of the day when the sun is above
`\Phi` degrees is approximated by

.. math::
   :label: eq-debm-toa-insolation

   \bar S_{\Phi} = \frac{S_0}{h_{\Phi}}  \left( \frac{\bar d}{d} \right)^2
   \left( h_{\Phi}  \sin(\phi)  \sin(\delta) + \cos(\phi)  \cos(\delta)  \sin(h_{\Phi}) \right)

where

- `S_0` is the solar constant :config:`surface.debm_simple.solar_constant` :cite:`Kopp2011`,
- `\bar d / d` is the ratio of the length of the semimajor axis of the Earth's orbit to
  the Earth-Sun distance,
- `h_{\Phi}` is the hour angle when the sun has an elevation angle of at least `\Phi`,
- `\phi` is the latitude
- `\delta` is the solar declination angle.

In short, `\bar S_{\Phi}` is a function of latitude, the factor `\bar d / d`, and the
solar declination angle `\delta`. In the "present day" case both `\bar d / d` and `\delta`
have the period of one year and are approximated using trigonometric expansions (see
:cite:`Liou2002`).

.. rubric:: Paleo simulations

Trigonometric expansions for `\bar d / d` and `\delta` mentioned above are not applicable
when modeling times far from present; in this case we use more general (and more
computationally expensive) formulas (:cite:`Liou2002`, chapter 2). Set
:config:`surface.debm_simple.paleo.enabled` to switch to using the "paleo" mode.

In this case

- the ratio `\bar d / d` is a function of the eccentricity of the Earth's orbit and the
  perihelion longitude,
- the solar declination `\delta` is a function of the eccentricity of the Earth's orbit,
  the perihelion longitude, and the Earth's obliquity.

The values of these are set using the following configuration parameters (prefix:
``surface.debm_simple.paleo.``):

.. pism-parameters::
   :prefix: surface.debm_simple.paleo.
   :exclude: .+(enabled|file|periodic)$

Alternatively, PISM can read in scalar time series of variables :var:`eccentricity`,
:var:`obliquity`, and :var:`perihelion_longitude` from a file specified using
:config:`surface.debm_simple.paleo.file`.

.. note::

   We provide a script (``examples/debm_simple/orbital_parameters.py``) that can be used
   to generate time series of these parameters using trigonometric expansions due to
   :cite:`Berger1978` with corrections made by the authors of the GISS GCM ModelE
   (expansion coefficients used in ``orbital_parameters.py`` come from the GISS ModelE
   version 2.1.2). See https://data.giss.nasa.gov/modelE/ar5plots/srorbpar.html for details.

   These expansions are considered to be valid for about 1 million years.


.. _sec-debm-simple-albedo:

Surface albedo
##############

To capture melt processes driven by changes in albedo without requiring a more
sophisticated surface process model (including the firn layer, for example), dEBM-simple
assumes that the surface albedo is a piecewise linear function of the modeled melt rate.

.. math::
   :label: eq-debm-surface-albedo

   \alpha = \max\left( \alpha_{\text{max}} + \alpha_{s} \cdot M, \alpha_{\text{min}}\right).

Here `M` is the estimated melt rate from the previous time step (meters liquid water
equivalent per second) and `\alpha_s` is a *negative* tuning parameter
(:config:`surface.debm_simple.albedo_slope`).

In this approach, the albedo decreases linearly with increasing melt from the maximum
value `\alpha_{\text{max}}` (the "fresh snow" albedo
:config:`surface.debm_simple.albedo_max`) for regions with no melting to the minimum value
`\alpha_{\text{min}}` (the "bare ice" albedo :config:`surface.debm_simple.albedo_min`).

Alternatively, albedo (variable :var:`surface_albedo`; no units) can be read from a file
specified using :config:`surface.debm_simple.albedo_input.file`.

.. note::

   - It is recommended to use monthly records of albedo in
     :config:`surface.debm_simple.albedo_input.file`.

   - The fresh snow albedo is also treated as a tuning parameter. Default values of
     `\alpha_{\text{max}}` and `\alpha_s` were obtained by fitting this approximation to
     the output of a regional climate model.

.. _sec-debm-simple-transmissivity:

Atmosphere transmissivity
#########################

dEBM-simple assumes that the transmissivity of the atmosphere `\tau_A` is a linear
function of the local surface altitude. Similarly to the albedo parameterization, the
default values of coefficients `a` and `b` below were obtained using linear regression
of an RCM output. This parameterization also relies on the assumption that no other
processes (e.g. changing mean cloud cover in a changing climate) affect `\tau_A`.

.. math::
   :label: eq-debm-transmissivity

   \tau_A = a + b\cdot z,

where `a` is set by :config:`surface.debm_simple.tau_a_intercept`, `b` is set by
:config:`surface.debm_simple.tau_a_slope`, and `z` is the ice surface altitude in meters.

.. _sec-debm-temperature-driven-melt:

Temperature-driven melt contribution
====================================

.. math::
   :label: eq-debm-temperature

   M_T = \frac{\Delta t_{\Phi}}{\Delta t \rho_{\text{w}} L_{\text{m}}}
   c_1 T_\text{eff},

where

.. math::
   :label: eq-debm-t-eff

   T_{\text{eff}}(T, \sigma) =
   \frac{1}{\sigma \sqrt{2 \pi}}
   \int_{T_{\text{pos}}}^{\infty} \xi\,
   \exp\left( -\frac{(\xi - T)^2}{2 \sigma^2}  \right)\, d\xi.

The "effective temperature" `T_{\text{eff}}` is the expected value of "positive"
excursions, i.e. excursions above the positivity threshold `T_{\text{pos}}`
(:config:`surface.debm_simple.positive_threshold_temp`, usually `0\!\phantom{|}^\circ
\text{C}`) of stochastic temperature variations added to the provided temperature
forcing.

Similarly to the PDD :ref:`sec-surface-pdd`, these stochastic variations are assumed to
follow the normal distribution with the mean of zero and the standard deviation `\sigma`
and are used to model the effect of daily temperature variations *not resolved* by this
model either because of the temporal resolution of the provided forcing or the chosen time
step length.

.. note::

   The standard deviation `\sigma` of added daily variations should be treated as a tuning
   parameter. The appropriate value may change depending on the application domain (for
   example: Greenland vs Antarctica), the temporal resolution of the air temperature
   forcing and lengths of time steps taken by dEBM-simple; see
   :config:`surface.debm_simple.max_evals_per_year`.

Here `\sigma` can be

- constant in time and space (the default; set using :config:`surface.debm_simple.std_dev`),
- read from a file containing the two dimensional variable :var:`air_temp_sd` that can be
  constant in time or time-dependent (units: *Kelvin*; specify the file name using
  :config:`surface.debm_simple.std_dev.file`), or
- parameterized as a function of air temperature `T`: `\sigma = \max(a\, (T -
  T_{\text{melting}}) + b, 0)` with `T_{\text{melting}} = 273.15` Kelvin.

These mechanisms are controlled by parameters with the prefix ``surface.debm_simple.std_dev.``:

.. pism-parameters::
   :prefix: surface.debm_simple.std_dev.

.. _sec-debm-simple-tuning:

Tuning parameters
=================

Default values of many parameters come from :cite:`Zeitzetal2021` and are appropriate for
Greenland; their values will need to change to use this model in other contexts. See Table
1 in :cite:`Garbe2023` for parameter values more appropriate in an Antarctic setting and
for the description of a calibration procedure that can be used to obtain some of these
values.

.. list-table:: Notable tuning parameters in the order of decreasing importance
   :header-rows: 1
   :widths: 1,1,3
   :name: tab-debm-simple-tuning-parameters

   * - Parameter
     - Equation
     - Configuration parameters

   * - `c_1`
     - :eq:`eq-debm-temperature`
     - :config:`surface.debm_simple.c1`

   * - `c_2`
     - :eq:`eq-debm-melt-components`
     - :config:`surface.debm_simple.c2`

   * - `T_{\text{min}}`
     - :eq:`eq-debm-melt`
     - :config:`surface.debm_simple.melting_threshold_temp`

   * - `\sigma`
     - :eq:`eq-debm-t-eff`
     - :config:`surface.debm_simple.std_dev` and others at the end of
       :ref:`sec-debm-temperature-driven-melt`

   * - `\theta`
     - :eq:`eq-debm-refreeze`
     - :config:`surface.debm_simple.refreeze`, :config:`surface.debm_simple.refreeze_ice_melt`

   * - `\alpha_{\text{max}}`
     - :eq:`eq-debm-surface-albedo`
     - :config:`surface.debm_simple.albedo_max`

   * - `\alpha_s`
     - :eq:`eq-debm-surface-albedo`
     - :config:`surface.debm_simple.albedo_slope`

   * - `a`
     - :eq:`eq-debm-transmissivity`
     - :config:`surface.debm_simple.tau_a_intercept`

   * - `b`
     - :eq:`eq-debm-transmissivity`
     - :config:`surface.debm_simple.tau_a_slope`

   * - `\Phi`
     - :eq:`eq-debm-insolation-melt`, :eq:`eq-debm-toa-insolation`
     - :config:`surface.debm_simple.phi`

.. _sec-surface-pik:

PIK
+++

:|options|: ``-surface pik``
:|variables|: :var:`climatic_mass_balance` |flux|,
              :var:`lat` (latitude), (degrees north)
:|implementation|: ``pism::surface::PIK``

This surface model component implements the setup used in :cite:`Martinetal2011`. The
:var:`climatic_mass_balance` is read from an input (``-i``) file; the ice surface
temperature is computed as a function of latitude (variable :var:`lat`) and surface
elevation (dynamically updated by PISM). See equation (1) in :cite:`Martinetal2011`.

.. _sec-surface-delta-t:

Scalar temperature offsets
++++++++++++++++++++++++++

:|options|: ``-surface ...,delta_T``
:|variables|: :var:`delta_T`
:|implementation|: ``pism::surface::Delta_T``

The time-dependent scalar offsets :var:`delta_T` are added to :var:`ice_surface_temp`
computed by a surface model.

Please make sure that :var:`delta_T` has the units of "``Kelvin``".

This modifier is identical to the corresponding atmosphere modifier, but applies offsets
at a different stage in the computation of top-surface boundary conditions needed by the
ice dynamics core.

.. rubric:: Parameters

Prefix: ``surface.delta_T.``

.. pism-parameters::
   :prefix: surface.delta_T.

.. _sec-surface-elevation-change:

Adjustments using modeled change in surface elevation
+++++++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-surface ...,elevation_change``
:|variables|: :var:`surface_altitude` (CF standard name),
:|implementation|: ``pism::surface::LapseRates``

The ``elevation_change`` modifier adjusts ice-surface temperature and surface mass balance
using modeled changes in surface elevation relative to a reference elevation read from a
file.

The surface temperature is modified using an elevation lapse rate
`\gamma_T =` :config:`surface.elevation_change.temperature_lapse_rate`. Here

.. math::
   \gamma_T = -\frac{dT}{dz}.

Two methods of adjusting the SMB are available:

- Scaling using an exponential factor

  .. math::

     \mathrm{SMB} = \mathrm{SMB_{input}} \cdot \exp(C \cdot \Delta T),

  where `C =` :config:`surface.elevation_change.smb.exp_factor` and `\Delta T` is the
  temperature difference produced by applying
  :config:`surface.elevation_change.temperature_lapse_rate`.

  This mechanisms increases the SMB by `100(\exp(C) - 1)` percent for each degree of
  temperature increase.

  To use this method, set :opt:`-smb_adjustment scale`.

- Elevation lapse rate for the SMB

  .. math::

     \mathrm{SMB} = \mathrm{SMB_{input}} - \Delta h \cdot \gamma_M,

  where `\gamma_M =` :config:`surface.elevation_change.smb.lapse_rate` and `\Delta h` is the
  difference between modeled and reference surface elevations.

  To use this method, set :opt:`-smb_adjustment shift`.

.. rubric:: Parameters

Prefix: ``surface.elevation_change.``.

.. pism-parameters::
   :prefix: surface.elevation_change.

.. _sec-surface-forcing:

Mass flux adjustment
++++++++++++++++++++
    
:|options|: ``-surface ...,forcing``
:|variables|: :var:`thk` (ice thickness), :var:`ftt_mask` (mask of zeros and ones; 1 where
              surface mass flux is adjusted and 0 elsewhere)
:|implementation|: ``pism::surface::ForceThickness``

The ``forcing`` modifier implements a surface mass balance adjustment mechanism which
forces the thickness of grounded ice to a target thickness distribution at the end of the
run. The idea behind this mechanism is that spinup of ice sheet models frequently requires
the surface elevation to come close to measured values at the end of a run. A simpler
alternative to accomplish this, namely option ``-no_mass``, represents an unmodeled,
frequently large, violation of the mass continuity equation.

In more detail, let `H_{\text{tar}}` be the target thickness. Let `H` be the
time-dependent model thickness. The surface model component described here produces the
term `M` in the mass continuity equation:

.. math::

   \frac{\partial H}{\partial t} = M - S - \nabla\cdot \mathbf{q}.

(Other details of this equation do not concern us here.) The ``forcing`` modifier causes
`M` to be adjusted by a multiple of the difference between the target thickness and
the current thickness,

.. math::

   \Delta M = \alpha (H_{\text{tar}} - H)

where `\alpha>0`. We are adding mass (`\Delta M>0`) where
`H_{\text{tar}} > H` and ablating where `H_{\text{tar}} < H`.

Option :opt:`-force_to_thickness_file` identifies the file containing the target ice
thickness field ``thk`` and the mask ``ftt_mask``. A basic run modifying surface model
``given`` would look like

.. code-block:: none

    pismr -i foo.nc -surface given,forcing -force_to_thickness_file bar.nc

In this case ``foo.nc`` contains fields :var:`climatic_mass_balance` and
:var:`ice_surface_temp`, as normal for ``-surface given``, and ``bar.nc`` contains fields
:var:`thk` which will serve as the target thickness and :var:`ftt_mask` which defines the
map plane area where this adjustment is applied. Option :opt:`-force_to_thickness_alpha`
adjusts the value of `\alpha`, which has a default value specified in the
:ref:`sec-parameter-list`.

In addition to this one can specify a multiplicative factor `C` used in areas where
the target thickness field has less than
:opt:`-force_to_thickness_ice_free_thickness_threshold` meters of ice;
`\alpha_{\text{ice free}} = C \times \alpha`. Use the
:opt:`-force_to_thickness_ice_free_alpha_factor` option to set `C`.

.. _sec-surface-anomaly:

Using climate data anomalies
++++++++++++++++++++++++++++
    
:|options|: :opt:`-surface ...,anomaly`
:|variables|: :var:`ice_surface_temp_anomaly`,
              :var:`climatic_mass_balance_anomaly` |flux|
:|implementation|: ``pism::surface::Anomaly``

This modifier implements a spatially-variable version of ``-surface ...,delta_T`` which
also applies time-dependent climatic mass balance anomalies.

See also ``-atmosphere ...,anomaly`` (section :ref:`sec-atmosphere-anomaly`), which is
similar but applies anomalies at the atmosphere level.

.. rubric:: Parameters

Prefix: ``surface.anomaly.``

.. pism-parameters::
   :prefix: surface.anomaly.

.. _sec-surface-cache:

The caching modifier
++++++++++++++++++++

:|options|: ``-surface ...,cache``
:|implementation|: ``pism::surface::Cache``
:|seealso|: :ref:`sec-ocean-cache`
    
This modifier skips surface model updates, so that a surface model is called no more than
every :config:`surface.cache.update_interval` 365-day "years". A time-step of `1` year is
used every time a surface model is updated.

This is useful in cases when inter-annual climate variability is important, but one year
differs little from the next. (Coarse-grid paleo-climate runs, for example.)

.. rubric:: Parameters

Prefix: ``surface.cache.``

.. pism-parameters::
   :prefix: surface.cache.

.. _sec-surface-no-gl-retreat:

Preventing grounding line retreat
+++++++++++++++++++++++++++++++++

:|options|: ``-surface ...,no_gl_retreat``
:|implementation|: ``pism::surface::NoGLRetreat``

This modifier adjust the surface mass balance to prevent the retreat of the grounding
line. See :ref:`sec-tillphi-optimization` for an application.

.. note::

   - This modifier *adds mass* in violation of mass conservation. Save the diagnostic
     :var:`no_gl_retreat_smb_adjustment` to get an idea about the amount added. Note,
     though, that this is an imperfect measure: it includes mass added to maintain
     non-negativity of ice thickness.

   - We assume that the sea level and the bed elevation remain constant throughout the
     simulation.

   - This does *not* prevent grounding line retreat caused by the thinning of the ice due
     to the melt at the base. Set :config:`geometry.update.use_basal_melt_rate` to "false"
     to ensure that basal melt has no effect on the position of the grounding line
