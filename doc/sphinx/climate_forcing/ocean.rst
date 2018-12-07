.. include:: shortcuts.txt

Ocean model components
----------------------

PISM Ocean model components provide sub-shelf ice temperature (:var:`shelfbtemp`) and
sub-shelf mass flux (:var:`shelfbmassflux`) to the ice dynamics core.

The sub-shelf ice temperature is used as a Dirichlet boundary condition in the energy
conservation code. The sub-shelf mass flux is used as a source in the mass-continuity
(transport) equation. Positive flux corresponds to ice loss; in other words, this
sub-shelf mass flux is a "melt rate".

.. contents::

.. _sec-ocean-constant:

Constant in time and space
++++++++++++++++++++++++++

:|options|: ``-ocean constant``
:|variables|: none
:|implementation|: ``pism::ocean::Constant``

.. note:: This is the default choice.

This ocean model component implements boundary conditions at the ice/ocean interface that
are constant *both* in space and time.

The sub-shelf ice temperature is set to pressure melting and the sub-shelf melt rate is
assumed to be proportional to the heat flux from the ocean into the ice (configuration
parameter :config:`ocean.sub_shelf_heat_flux_into_ice`).

Alternatively, the sub-shelf melt rate in meters per year can be set using the
:opt:`-shelf_base_melt_rate` command-line option.

.. _sec-ocean-given:

Reading forcing data from a file
++++++++++++++++++++++++++++++++

:|options|: ``-ocean given``
:|variables|: :var:`shelfbtemp` Kelvin,
              :var:`shelfbmassflux`  |flux|
:|implementation|: ``pism::ocean::Given``

This ocean model component reads sub-shelf ice temperature :var:`shelfbtemp` and the
sub-shelf mass flux :var:`shelfbmassflux` from a file. It takes the following command-line
options.

- :opt:`-ocean_given_file`: sets the name of the file to read forcing data from. The file
  may contain several records. If only one record is provided it is interpreted as
  time-independent.
- :opt:`-ocean_given_reference_year` specifies the reference date; see section :ref:`sec-periodic-forcing`.
- :opt:`-ocean_given_period` specifies the length of the period of the forcing data, in
  model years; see section :ref:`sec-periodic-forcing`.

Variables :var:`shelfbtemp` and :var:`shelfbmassflux` may be time-dependent. (The ``-ocean
given`` component is very similar to ``-surface given`` and ``-atmosphere given``.)

.. _sec-ocean-pik:

PIK
+++

:|options|: ``-ocean pik``
:|variables|: none
:|implementation|: ``pism::ocean::PIK``

This ocean model component implements the ocean forcing setup used in :cite:`Martinetal2011`.
The sub-shelf ice temperature is set to pressure-melting; the sub-shelf mass flux
computation follows :cite:`BeckmannGoosse2003`.

It takes one command-line option:

- :opt:`-meltfactor_pik`: a melt factor `F_{\mathrm{melt}}` in sub-shelf-melting
  parameterization, see equation (5) in :cite:`Martinetal2011`.

.. _sec-ocean-th:

Basal melt rate and temperature from thermodynamics in boundary layer
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-ocean th``
:|variables|: :var:`theta_ocean` (absolute potential ocean temperature), [Kelvin],
              :var:`salinity_ocean` (salinity of the adjacent ocean), [g/kg]
:|implementation|: ``pism::ocean::GivenTH``

This ocean model component derives basal melt rate and basal temperature from
thermodynamics in a boundary layer at the base of the ice shelf. It uses a set of three
equations describing

#. the energy flux balance,
#. the salt flux balance,
#. the pressure and salinity dependent freezing point in the boundary layer.

This model is described in :cite:`HollandJenkins1999` and :cite:`Hellmeretal1998`.

Inputs are potential temperature (variable :var:`theta_ocean`) and salinity (variable
:var:`salinity_ocean`) read from a file.

No ocean circulation is modeled, so melt water computed by this model is not fed back into
the surrounding ocean.

This implementation uses different approximations of the temperature gradient at the base
of an ice shelf column depending on whether there is sub-shelf melt, sub-shelf freeze-on,
or neither (see :cite:`HollandJenkins1999` for details).

It takes two command-line option:

- :opt:`-ocean_th_file`: specifies the NetCDF file providing potential temperature and
  salinity fields.
- :opt:`-clip_shelf_base_salinity`: if this is set (which is the default), the sub-shelf
  salinity is clipped so that it stays in the `[4, 40]` psu range. This is done to
  ensure that we stay in the range of applicability of the melting point temperature
  parameterization; see :cite:`HollandJenkins1999`. To disable salinity clipping, use the
  :opt:`-no_clip_shelf_base_salinity` option or set the configuration parameter
  :config:`ocean.three_equation_model_clip_salinity`  to "no".

.. _sec-pico:

PICO
++++

:|options|: ``-ocean pico``
:|variables|: :var:`theta_ocean` (potential ocean temperature), [Kelvin],

              :var:`salinity_ocean` (salinity of the adjacent ocean), [g/kg],

              :var:`basins` (mask of large-scale ocean basins that ocean input is averaged over), [integer]
:|implementation|: ``pism::ocean::Pico``

The PICO model provides sub-shelf melt rates and temperatures consistent with the vertical
overturning circulation in ice shelf cavities that drives the exchange with open ocean
water masses. It is based on the ocean box model of :cite:`OlbersHellmer2010` and includes
a geometric approach which makes it applicable to ice shelves that evolve in two
horizontal dimensions. For each ice shelf, PICO solves the box model equations describing
the transport between coarse ocean boxes. It applies a boundary layer melt formulation
:cite:`HellmerOlbers1989`, :cite:`HollandJenkins1999`. The overturning circulation is
driven by the ice-pump :cite:`LewisPerkin1986`: melting at the ice-shelf base reduces the
density of ambient water masses. Buoyant water rising along the shelf base draws in ocean
water at depth, which flows across the continental shelf towards the deep grounding lines.
The model captures this circulation by defining consecutive boxes following the flow
within the ice shelf cavity, with the first box adjacent to the grounding line. The
extents of the ocean boxes are computed adjusting to the evolving grounding lines and
calving fronts. Open ocean properties in front of the shelf as well as the geometry of the
shelf determine basal melt rate and basal temperature at each grid point.

The main equations reflect the

#. heat and salt balance for each ocean box in contact with the ice shelf base,
#. overturning flux driven by the density difference between open-ocean and grounding-line box,
#. boundary layer melt formulation.

The PICO model is described in detail in :cite:`ReeseAlbrecht2018`.

Inputs are potential temperature (variable :var:`theta_ocean`), salinity (variable
:var:`salinity_ocean`) and ocean basin mask (variable :var:`basins`). Variables
:var:`theta_ocean` and :var:`salinity_ocean` may be time-dependent.

Forcing ocean temperature and salinity are taken from the water masses that occupy the sea
floor in front of the ice shelves, which extends down to a specified continental shelf
depth (see :config:`ocean.pico.continental_shelf_depth`). These water masses are
transported by the overturning circulation into the ice shelf cavity and towards the
grounding line. The basin mask defines regions of similar, large-scale ocean conditions;
each region is marked with a distinct positive integer. In PICO, ocean input temperature
and salinity are averaged on the continental shelf within each basins. For each ice shelf,
the input values of the overturning circulation are calculated as an area-weighted average
over all basins that intersect the ice shelf. If ocean input parameters cannot be
identified, standard values are used (**Warning:** this could strongly influence melt
rates computed by PICO). In regions where the PICO geometry cannot be identified,
:cite:`BeckmannGoosse2003` is applied.

PICO has one command-line option and 7 configuration parameters:

- :opt:`-ocean_pico_file`: specifies the NetCDF file containing potential temperature
  (:var:`theta_ocean`), salinity (:var:`salinity_ocean`) and ocean basins (:var:`basins`).
- :config:`ocean.pico.heat_exchange_coefficent` sets the coefficient for turbulent heat
  exchange from the ambient ocean across the boundary layer beneath the ice shelf base.
- :config:`ocean.pico.overturning_coefficent`: sets the coefficient in the overturning
  parameterization.
- :config:`ocean.pico.number_of_boxes`: For each ice shelf the number of ocean boxes is
  determined by interpolating between 1 and number_of_boxes depending on its size and
  geometry such that larger ice shelves are resolved with more boxes; a value of 5 is
  suitable for the Antarctic setup.
- :config:`ocean.pico.number_of_basins`
- :config:`ocean.pico.exclude_ice_rises`: If set to true, grounding lines of ice rises are
  excluded in the geometrical routines that determine the ocean boxes; using this option
  is recommended.
- :config:`ocean.pico.continental_shelf_depth`: specifies the depth up to which oceanic
  input temperatures and salinities are averaged over the continental shelf areas in front
  of the ice shelf cavities.
- :config:`ocean.pico.maximum_ice_rise_area`: specifies an area threshold that separates
  ice rises from continental regions.

.. _sec-ocean-delta-sl:

Scalar sea level offsets
++++++++++++++++++++++++

:|options|: :opt:`-sea_level ...,delta_sl`
:|variables|: :var:`delta_SL` (meters)
:|implementation|: ``pism::ocean::sea_level::Delta_SL``

The ``delta_sl`` modifier implements sea level forcing using scalar offsets.

It takes the following command-line options:

- :opt:`-ocean_delta_sl_file`: specifies the name of the file containing forcing data.
  This file has to contain the :var:`delta_SL` variable using units "meters" or
  equivalent.
- :opt:`-ocean_delta_sl_period` specifies the length of the period of the forcing data, in
  model years; see section :ref:`sec-periodic-forcing`.
- :opt:`-ocean_delta_sl_reference_year` specifies the reference date; see section
  :ref:`sec-periodic-forcing`.

.. _sec-ocean-delta-sl-2d:

Two-dimensional sea level offsets
+++++++++++++++++++++++++++++++++

:|options|: :opt:`-sea_level ...,delta_sl_2d`
:|variables|: :var:`delta_SL` (meters)
:|implementation|: ``pism::ocean::sea_level::Delta_SL_2D``

The ``delta_sl`` modifier implements sea level forcing using time-dependent and
spatially-variable offsets.

It uses the following configuration parameters:

- :config:`ocean.delta_sl_2d.file`: specifies the name of the file containing forcing
  data. This file has to contain the :var:`delta_SL` variable using units "meters" or
  equivalent.
- :config:`ocean.delta_sl_2d.period` specifies the length of the period of the forcing
  data, in model years; see section :ref:`sec-periodic-forcing`.
- :config:`ocean.delta_sl_2d.reference_year` specifies the reference date; see section
  :ref:`sec-periodic-forcing`.

.. _sec-ocean-delta-t:

Scalar sub-shelf temperature offsets
++++++++++++++++++++++++++++++++++++


:|options|: :opt:`-ocean ...,delta_T`
:|variables|: :var:`delta_T` (Kelvin)
:|implementation|: ``pism::ocean::Delta_T``

This modifier implements forcing using sub-shelf ice temperature offsets.

It takes the following command-line options:

- :opt:`-ocean_delta_T_file`: specifies the name of the file containing forcing data. This
  file has to contain the :var:`delta_T` variable using units of "Kelvin" or equivalent.
- :opt:`-ocean_delta_T_period` specifies the length of the period of the forcing data, in
  model years; see section :ref:`sec-periodic-forcing`.
- :opt:`-ocean_delta_T_reference_year` specifies the reference date; see section
  :ref:`sec-periodic-forcing`.

.. _sec-ocean-delta-smb:

Scalar sub-shelf mass flux offsets
++++++++++++++++++++++++++++++++++

:|options|: ``-ocean ...,delta_SMB``
:|variables|: :var:`delta_SMB` |flux|
:|implementation|: ``pism::ocean::Delta_SMB``

This modifier implements forcing using sub-shelf mass flux (melt rate) offsets.

It takes the following command-line options:

- :opt:`-ocean_delta_SMB_file`: specifies the name of the file containing forcing data.
  This file has to contain the :var:`delta_SMB` variable using units |flux| or equivalent.
- :opt:`-ocean_delta_SMB_period` specifies the length of the period of the forcing data,
  in model years; see section :ref:`sec-periodic-forcing`.
- :opt:`-ocean_delta_SMB_reference_year` specifies the reference date; see section
  :ref:`sec-periodic-forcing`.

.. _sec-ocean-frac-smb:

Scalar sub-shelf mass flux fraction offsets
+++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-ocean ...,frac_SMB``
:|variables|: :var:`frac_SMB` [1]
:|implementation|: ``pism::ocean::Frac_SMB``

This modifier implements forcing using sub-shelf mass flux (melt rate) fraction offsets.

It takes the following command-line options:

- :opt:`-ocean_frac_SMB_file`: specifies the name of the file containing forcing data.
  This file has to contain the :var:`frac_SMB` variable.
- :opt:`-ocean_frac_SMB_period` specifies the length of the period of the forcing data, in
  model years; see section :ref:`sec-periodic-forcing`.
- :opt:`-ocean_frac_SMB_reference_year` specifies the reference date; see section
  :ref:`sec-periodic-forcing`.

.. _sec-ocean-anomaly:

Two-dimensional sub-shelf mass flux offsets
+++++++++++++++++++++++++++++++++++++++++++

:|options|: :opt:`-ocean ...,anomaly`
:|variables|: :var:`shelf_base_mass_flux_anomaly` |flux|
:|implementation|: ``pism::ocean::Anomaly``

This modifier implements a spatially-variable version of ``-ocean ...,delta_SMB`` which
applies time-dependent shelf base mass flux anomalies, as used for initMIP or LARMIP
model intercomparisons.

It takes the following command-line options:

- :opt:`-ocean_anomaly_file` specifies a file containing the variable
  :var:`shelf_base_mass_flux_anomaly`.
- :opt:`-ocean_anomaly_period` (years) specifies the period of the forcing data, in
  model years; see :ref:`sec-periodic-forcing`
- :opt:`-ocean_anomaly_reference_year` specifies the reference year; see
  :ref:`sec-periodic-forcing`

  See also to ``-atmosphere ...,anomaly`` or
  ``-surface ...,anomaly`` (section :ref:`sec-surface-anomaly`)
  which is similar, but applies anomalies at the atmosphere or surface level,
  respectively.

.. _sec-ocean-frac-mbp:

Scalar melange back pressure fraction
+++++++++++++++++++++++++++++++++++++

:|options|: :opt:`-ocean ...,frac_MBP`
:|variables|: :var:`frac_MBP`
:|implementation|: ``pism::ocean::Frac_MBP``

This modifier implements forcing using melange back pressure fraction offsets. The
variable :var:`frac_MBP` should take on values from 0 to 1; it is understood as the
fraction of the maximum melange back pressure possible at a given location. (We assume
that melange back pressure cannot exceed the pressure of the ice column at a calving
front.)

Please see :ref:`sec-model-melange-pressure` for details.

This modifier takes the following command-line options:

- :opt:`-ocean_frac_MBP_file`: specifies the name of the file containing forcing data.
  This file has to contain the :var:`frac_MBP` variable using units of "1" (a
  dimensionless parameter)
- :opt:`-ocean_frac_MBP_period` specifies the length of the period of the forcing data, in
  model years; see section :ref:`sec-periodic-forcing`.
- :opt:`-ocean_frac_MBP_reference_year` specifies the reference date; see section
  :ref:`sec-periodic-forcing`.

.. _sec-ocean-cache:

The caching modifier
++++++++++++++++++++

:|options|: :opt:`-ocean ...,cache`
:|implementation|: ``pism::ocean::Cache``
:|seealso|: :ref:`sec-surface-cache`

This modifier skips ocean model updates, so that a ocean model is called no more than
every :opt:`-ocean_cache_update_interval` years. A time-step of `1` year is used every
time a ocean model is updated.

This is useful in cases when inter-annual climate variability is important, but one year
differs little from the next. (Coarse-grid paleo-climate runs, for example.)

It takes the following options:

- :opt:`-ocean_cache_update_interval` (*years*) Specifies the minimum interval between
  updates. PISM may take longer time-steps if the adaptive scheme allows it, though.
