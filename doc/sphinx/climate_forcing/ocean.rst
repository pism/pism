.. include:: shortcuts.txt

Ocean model components
----------------------

PISM's ocean model components provide sub-shelf ice temperature (:var:`shelfbtemp`) and
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
controlled by :config:`ocean.constant.melt_rate`.

.. _sec-ocean-given:

Reading forcing data from a file
++++++++++++++++++++++++++++++++

:|options|: ``-ocean given``
:|variables|: :var:`shelfbtemp` Kelvin,
              :var:`shelfbmassflux`  |flux|
:|implementation|: ``pism::ocean::Given``

This ocean model component reads sub-shelf ice temperature :var:`shelfbtemp` and the
sub-shelf mass flux :var:`shelfbmassflux` from a file.

It uses the following parameters:

.. pism-parameters::
   :prefix: ocean.given.

Variables :var:`shelfbtemp` and :var:`shelfbmassflux` may be time-dependent. (The ``-ocean
given`` component is very similar to ``-surface given`` and ``-atmosphere given``.)

.. _sec-ocean-pik:

PIK
+++

:|options|: ``-ocean pik``
:|variables|: none
:|implementation|: ``pism::ocean::PIK``

This ocean model component implements the ocean forcing setup used in
:cite:`Martinetal2011`. The sub-shelf ice temperature is set to pressure-melting; the
sub-shelf mass flux computation follows :cite:`BeckmannGoosse2003`.

It uses the following parameters:

.. pism-parameters::
   :prefix: ocean.pik_

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
:var:`salinity_ocean`) read from a file. A constant salinity (see
:config:`constants.sea_water.salinity`) is used if the input file does not contain
:var:`salinity_ocean`.

No ocean circulation is modeled, so melt water computed by this model is not fed back into
the surrounding ocean.

This implementation uses different approximations of the temperature gradient at the base
of an ice shelf column depending on whether there is sub-shelf melt, sub-shelf freeze-on,
or neither (see :cite:`HollandJenkins1999` for details).

It uses the following parameters:

.. pism-parameters::
   :prefix: ocean.th.

.. note::

   If :config:`ocean.th.clip_salinity` is set (the default), the sub-shelf salinity is
   clipped so that it stays in the `[4, 40]` psu range. This is done to ensure that we
   stay in the range of applicability of the melting point temperature parameterization;
   see :cite:`HollandJenkins1999`.

   Set :config:`ocean.th.clip_salinity` to ``false`` if restricting salinity is not
   appropriate.

See :ref:`sec-ocean-th-details` for implementation details.

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
over all basins that intersect the ice shelf. Only those basins are considered in the average, 
in which the ice shelf has in fact a connection to the ocean. Large ice shelves, that cover 
across two basins, that do not share an ocean boundary, are considered as two separate ice 
shelves with individual ocean inputs. If ocean input parameters cannot be
identified, standard values are used (**Warning:** this could strongly influence melt
rates computed by PICO). In regions where the PICO geometry cannot be identified,
:cite:`BeckmannGoosse2003` is applied.

PICO uses the following configuration parameters (prefix: ``ocean.pico.``):

.. pism-parameters::
   :prefix: ocean.pico.

.. _sec-ocean-delta-sl:

Scalar sea level offsets
++++++++++++++++++++++++

:|options|: :opt:`-sea_level ...,delta_sl`
:|variables|: :var:`delta_SL` (meters)
:|implementation|: ``pism::ocean::sea_level::Delta_SL``

The ``delta_sl`` modifier implements sea level forcing using scalar offsets.

It uses the following parameters:

.. pism-parameters::
   :prefix: ocean.delta_sl.

.. _sec-ocean-delta-sl-2d:

Two-dimensional sea level offsets
+++++++++++++++++++++++++++++++++

:|options|: :opt:`-sea_level ...,delta_sl_2d`
:|variables|: :var:`delta_SL` (meters)
:|implementation|: ``pism::ocean::sea_level::Delta_SL_2D``

The ``delta_sl`` modifier implements sea level forcing using time-dependent and
spatially-variable offsets.

It uses the following configuration parameters:

.. pism-parameters::
   :prefix: ocean.delta_sl_2d.

.. _sec-ocean-delta-t:

Scalar sub-shelf temperature offsets
++++++++++++++++++++++++++++++++++++


:|options|: :opt:`-ocean ...,delta_T`
:|variables|: :var:`delta_T` (Kelvin)
:|implementation|: ``pism::ocean::Delta_T``

This modifier implements forcing using sub-shelf ice temperature offsets.

It uses the following parameters:

.. pism-parameters::
   :prefix: ocean.delta_T.

.. _sec-ocean-delta-smb:

Scalar sub-shelf mass flux offsets
++++++++++++++++++++++++++++++++++

:|options|: ``-ocean ...,delta_SMB``
:|variables|: :var:`delta_SMB` |flux|
:|implementation|: ``pism::ocean::Delta_SMB``

This modifier implements forcing using sub-shelf mass flux (melt rate) offsets.

It uses the following parameters:

.. pism-parameters::
   :prefix: ocean.delta_mass_flux.

.. _sec-ocean-frac-smb:

Scalar sub-shelf mass flux fraction offsets
+++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-ocean ...,frac_SMB``
:|variables|: :var:`frac_SMB` [1]
:|implementation|: ``pism::ocean::Frac_SMB``

This modifier implements forcing using sub-shelf mass flux (melt rate) fraction offsets.

It uses the following parameters:

.. pism-parameters::
   :prefix: ocean.frac_mass_flux.

.. _sec-ocean-anomaly:

Two-dimensional sub-shelf mass flux offsets
+++++++++++++++++++++++++++++++++++++++++++

:|options|: :opt:`-ocean ...,anomaly`
:|variables|: :var:`shelf_base_mass_flux_anomaly` |flux|
:|implementation|: ``pism::ocean::Anomaly``

This modifier implements a spatially-variable version of ``-ocean ...,delta_SMB`` which
applies time-dependent shelf base mass flux anomalies, as used for initMIP or LARMIP
model intercomparisons.

It uses the following parameters:

.. pism-parameters::
   :prefix: ocean.anomaly.

See also to ``-atmosphere ...,anomaly`` or ``-surface ...,anomaly`` (section
:ref:`sec-surface-anomaly`) which is similar, but applies anomalies at the atmosphere or
surface level, respectively.

.. _sec-ocean-delta-mbp:

Scalar melange back pressure offsets
++++++++++++++++++++++++++++++++++++

:|options|: :opt:`-ocean ...,delta_MBP`
:|variables|: :var:`delta_MBP` [Pascal]
:|implementation|: ``pism::ocean::Delta_MBP``

The scalar time-dependent variable :var:`delta_MBP` (units: Pascal) has the meaning of the
melange back pressure `\sigma_b` in :cite:`Krug2015`. It is assumed that `\sigma_b` is
applied over the thickness of melange `h` specified using
:config:`ocean.delta_MBP.melange_thickness`.

To convert to the average pressure over the ice front thickness, we compute

.. math::
   :label: eq-melange-pressure

   \bar p_{\text{melange}} = \frac{\sigma_b\cdot h}{H},

where `H` is ice thickness.

See :ref:`sec-model-melange-pressure` for details.

This modifier uses the following configuration parameters:

.. pism-parameters::
   :prefix: ocean.delta_MBP.

.. _sec-ocean-frac-mbp:

Melange back pressure as a fraction of pressure difference
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: :opt:`-ocean ...,frac_MBP`
:|variables|: :var:`frac_MBP`
:|implementation|: ``pism::ocean::Frac_MBP``

This modifier implements forcing using melange back pressure fraction (scaling).

Here we assume that the total vertically-averaged back pressure at an ice margin cannot
exceed the vertically-averaged ice pressure at the same location:

.. math::

   \bar p_{\text{water}} + \bar p_{\text{melange}} &\le \bar p_{\text{ice}},\, \text{or}

   \bar p_{\text{melange}} &\le \bar p_{\text{ice}} - \bar p_{\text{water}}.

We introduce `\lambda \in [0, 1]` such that

.. math::

   \bar p_{\text{melange}} = \lambda (\bar p_{\text{ice}} - \bar p_{\text{water}}).


The scalar time-dependent variable :var:`frac_MBP` should take on values between 0 and 1
and has the meaning of `\lambda` above.

Please see :ref:`sec-model-melange-pressure` for details.

This modifier uses the following configuration parameters:

.. pism-parameters::
   :prefix: ocean.frac_MBP.

.. _sec-ocean-cache:

The caching modifier
++++++++++++++++++++

:|options|: :opt:`-ocean ...,cache`
:|implementation|: ``pism::ocean::Cache``
:|seealso|: :ref:`sec-surface-cache`

This modifier skips ocean model updates, so that a ocean model is called no more than
every :config:`ocean.cache.update_interval` 365-day "years". A time-step of `1` year
(respecting the chosen calendar) is used every time a ocean model is updated.

This is useful in cases when inter-annual climate variability is important, but one year
differs little from the next. (Coarse-grid paleo-climate runs, for example.)

It uses the following parameters:

.. pism-parameters::
   :prefix: ocean.cache.
