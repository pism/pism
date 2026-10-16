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
:|variables|: :var:`shelfbtemp` kelvin,
              :var:`shelfbmassflux`  |flux|
:|implementation|: ``pism::ocean::Given``

This ocean model component reads sub-shelf ice temperature :var:`shelfbtemp` and the
sub-shelf mass flux :var:`shelfbmassflux` from a file.

Variables :var:`shelfbtemp` and :var:`shelfbmassflux` may be time-dependent. (The ``-ocean
given`` component is very similar to ``-surface given`` and ``-atmosphere given``.)

.. rubric:: Parameters

Prefix: ``ocean.given.``

.. pism-parameters::
   :prefix: ocean.given.

.. _sec-ocean-pik:

PIK
+++

:|options|: ``-ocean pik``
:|variables|: none
:|implementation|: ``pism::ocean::PIK``

This ocean model component implements the ocean forcing setup used in
:cite:`Martinetal2011`. The sub-shelf ice temperature is set to pressure-melting; the
sub-shelf mass flux computation follows :cite:`BeckmannGoosse2003`.

.. rubric:: Parameters

Prefix: ``ocean.pik_``

.. pism-parameters::
   :prefix: ocean.pik_

.. _sec-ocean-th:

Basal melt rate and temperature from thermodynamics in boundary layer
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-ocean th``
:|variables|: :var:`theta_ocean` (absolute potential ocean temperature), [kelvin],
              :var:`salinity_ocean` (salinity of the adjacent ocean), [g/kg]
:|implementation|: ``pism::ocean::GivenTH``

This ocean model component derives basal melt rate and basal temperature from
thermodynamics in a boundary layer at the base of the ice shelf. It uses a set of three
equations describing

#. the energy flux balance,
#. the salt flux balance,
#. the pressure and salinity dependent freezing point in the boundary layer.

Sub-shelf circulation is not modeled. This model is described in
:cite:`HollandJenkins1999` and :cite:`Hellmeretal1998`.

Inputs are two-dimensional, possibly time-dependent potential temperature (variable
:var:`theta_ocean`) and salinity (variable :var:`salinity_ocean`) read from a file
:config:`ocean.th.file`. A constant salinity (see :config:`constants.sea_water.salinity`)
is used if the input file does not contain :var:`salinity_ocean`.

This implementation uses different approximations of the temperature gradient at the base
of an ice shelf column depending on whether there is sub-shelf melt, sub-shelf freeze-on,
or neither (see :cite:`HollandJenkins1999` and :ref:`sec-ocean-th-details` for details).

.. rubric:: Parameters

Prefix: ``ocean.th.``

.. pism-parameters::
   :prefix: ocean.th.

.. note::

   If :config:`ocean.th.clip_salinity` is set (the default), the sub-shelf salinity is
   clipped so that it stays in the `[4, 40]` psu range. This is done to ensure that we
   stay in the range of applicability of the melting point temperature parameterization;
   see :cite:`HollandJenkins1999`.

   Set :config:`ocean.th.clip_salinity` to ``false`` if restricting salinity is not
   appropriate.

.. _sec-pico:

PICO
++++

:|options|: ``-ocean pico``
:|variables|: :var:`theta_ocean` (potential ocean temperature), [kelvin],

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

Inputs are two-dimensional, possibly time-dependent potential temperature (variable
:var:`theta_ocean`), salinity (variable :var:`salinity_ocean`) and a constant in time
ocean basin mask (variable :var:`basins`) read from a file :config:`ocean.pico.file`.

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

.. rubric:: Parameters

Prefix: ``ocean.pico.``

.. pism-parameters::
   :prefix: ocean.pico.

.. _sec-plume:

Plume
+++++

:|options|: ``-ocean plume``
:|variables|: :var:`theta_ocean` (potential ocean temperature, or thermal forcing; see below), [kelvin],

              :var:`salinity_ocean` (salinity of the adjacent ocean), [g/kg]
:|implementation|: ``pism::ocean::Plume``

The plume model computes sub-shelf melt rates from the theory of buoyant meltwater plumes
that rise along the base of an ice shelf, starting at the grounding line. It combines

- the *ambient* melt parameterization of :cite:`Lazeroms2018`, as adapted by
  :cite:`Pelle2019`, which is based on the one-dimensional plume model of
  :cite:`Jenkins1991`: the plume is driven by the thermal forcing of the ambient ocean and
  by the melt water it entrains along the way, and
- the *subglacial discharge* extension of :cite:`Pelle2023`: fresh water that leaves the
  glacier bed at the grounding line makes the plume more buoyant and increases melt near
  the outflow.

The ambient temperature `T_a` and salinity `S_a` are read from a file at every floating
grid cell. This is the model to use when the forcing already resolves the near-glacier
ocean state -- for example thermal forcing extrapolated into the fjords, as in the ISMIP6
and ISMIP7 Greenland protocols. (:ref:`sec-picop` uses the same plume physics, but obtains
`T_a` and `S_a` from the PICO box model.)

.. note::

   Positive melt rates correspond to ice loss. The plume is defined on floating ice only.
   The model needs ice velocities from the stress balance model. It also needs a
   hydrology model that provides the subglacial water flux (``-hydrology routing``, see
   :ref:`sec-hydrology-routing`) unless the discharge contribution is turned off using
   :config:`ocean.plume.add_fresh_water_melt`.

.. rubric:: Ambient conditions

If :config:`ocean.plume.temperature_as_thermal_forcing` is set, :var:`theta_ocean` is
interpreted as thermal forcing `\mathrm{TF}` (temperature above the freezing point, in
degrees Celsius) and the ambient temperature is

.. math::
   :label: eq-plume-ta

   T_a = T_f(S_a, z_{\mathrm{gl}}) + \mathrm{TF},

so that the plume is driven by exactly the thermal forcing in the file. Otherwise
:var:`theta_ocean` is used as `T_a` directly. In both cases `T_a` is bounded below by the
surface freezing point `T_f(S_a, 0)`.

The freezing point of sea water with salinity `S` at elevation `z` (negative below sea
level) is

.. math::
   :label: eq-plume-freezing-point

   T_f(S, z) = \lambda_1 S + \lambda_2 + \lambda_3 z,

with the coefficients `\lambda_{1,2,3}` set by ``ocean.plume.freezing_point_*``. The
sub-shelf ice temperature (:var:`shelfbtemp`) is set to `T_f(S_a, z_{\mathrm{b}})`, where
`z_{\mathrm{b}}` is the elevation of the shelf base.

.. rubric:: Geometry: grounding line elevation and basal slope

A plume originates at the grounding line, so the elevation `z_{\mathrm{gl}}` of the point
where the water that reaches a given floating cell left the grounding line is needed
everywhere under the shelf. It is obtained by transporting the bed elevation at the
grounding line along the ice flow (the modeled depth-averaged ice velocity `\mathbf{v}`):

.. math::
   :label: eq-grounding-line-depth

    \mathbf{v} \cdot \nabla z_{\mathrm{gl}} = 0 \text{ on floating ice},\qquad
    z_{\mathrm{gl}} = z_{\mathrm{bed}} \text{ at the grounding line}.

PISM solves this with a semi-Lagrangian (backward-characteristics) scheme, iterated until
the maximum relative change falls below a tolerance. By default each iteration starts
from the previous solution (see :config:`ocean.plume.transport_warm_start`). Because a
plume cannot originate above the local shelf base, `z_{\mathrm{gl}}` is limited to
`\min(z_{\mathrm{gl}}, z_{\mathrm{b}})`.

The local basal slope is `\alpha = \arctan |\nabla z_{\mathrm{b}}|`. Finite differences
are not taken across grounding lines or calving fronts.

.. rubric:: Ambient melt

With `\Delta T = T_a - T_f(S_a, z_{\mathrm{gl}})` define the effective heat exchange
coefficient

.. math::
   :label: eq-plume-gamma

   \Gamma_{TS} = \frac{C_d^T}{\sqrt{C_d}}
   \left[\gamma_1 + \gamma_2\,
   \frac{\Delta T}{\lambda_3}\,
   \frac{E_0 \sin\alpha}{C_{d,TS}^0 + E_0\sin\alpha}\right],

where `C_d` is the drag coefficient, `C_d^T` the turbulent heat exchange coefficient,
`E_0` the entrainment coefficient, and `\gamma_1`, `\gamma_2`, `C_{d,TS}^0` are
parameters (see the table below). Write `C_{d,TS} = \sqrt{C_d}\,\Gamma_{TS}`. The geometric
factor is

.. math::
   :label: eq-plume-geometry

   g(\alpha) &= G_1 G_2 G_3,

   G_1 &= \sqrt{\frac{\sin\alpha}{C_d + E_0\sin\alpha}},

   G_2 &= \sqrt{\frac{C_{d,TS}}{C_{d,TS} + E_0\sin\alpha}},

   G_3 &= \frac{E_0\sin\alpha}{C_{d,TS} + E_0\sin\alpha}.

The melt scale is

.. math::
   :label: eq-plume-melt-scale

   M = M_0\, g(\alpha)\, \Delta T^{\,\beta},

with `M_0` (:config:`ocean.plume.melt_rate_parameter`) and the thermal forcing exponent
`\beta` (:config:`ocean.plume.power_beta`, default 2). The position along the plume is
described by the dimensionless coordinate

.. math::
   :label: eq-dimensionless-coordinate

   \hat X = \mathrm{clip}\!\left(\frac{z_{\mathrm{b}} - z_{\mathrm{gl}}}{l},\, 0,\, 1\right),
   \qquad
   l = f\,\frac{\Delta T}{\lambda_3}\,
   \frac{x_0 C_{d,TS} + E_0\sin\alpha}{x_0\,(C_{d,TS} + E_0\sin\alpha)},

where `x_0` is :config:`ocean.plume.dimensionless_scaling_factor` and `f` is
:config:`ocean.plume.length_scale_factor` (default 1, which reproduces the published
parameterization). The ambient melt rate is then

.. math::
   :label: eq-melt-rate

   \dot m_a = \hat M(\hat X)\, M,

where `\hat M` is the dimensionless melt curve of :cite:`Lazeroms2018` (a degree-11
polynomial, using the coefficients from the corrigendum). It peaks at `\hat X \approx
0.18` and vanishes at `\hat X = x_0 = 0.56`, so melt is largest some distance downstream
of the grounding line.

.. note::

   The melt curve was fitted to Antarctic cavities where the draft rises by about 1 km.
   On shelves that rise by a few hundred meters (e.g. in Greenland fjords) `l` may be
   too large, so that melt increases toward the ice front instead of peaking near the
   grounding line. Use :config:`ocean.plume.length_scale_factor` (a value near 0.15 gives
   `l` of order 0.5 km) to correct this.

.. rubric:: Subglacial discharge

If :config:`ocean.plume.add_fresh_water_melt` is set, the melt caused by fresh water that
enters the ocean at the grounding line is added following :cite:`Pelle2023`. For a
discharge flux `q_{sg}` (volume of water per unit width of the grounding line and unit
time, `\mathrm{m^2\,s^{-1}}`) the plume is buoyant if

.. math::
   :label: eq-plume-buoyancy

   \Delta\rho_i = \beta_S S_a - \beta_T \Delta T > 0

(`\beta_S` and `\beta_T` are the haline contraction and thermal expansion coefficients;
the discharge has zero salinity). In this case the discharge-driven melt rate is

.. math::
   :label: eq-plume-fresh-water-melt

   \dot m_{fw} = \frac{c_p}{L_f}\, C_{d,TS}^0\, G_1\,
   \left(G_2\, g\, q_{sg}\, \Delta\rho_i\right)^{\nu} \Delta T,

where `c_p` and `L_f` are the specific heat capacity and latent heat of fusion of fresh
water, `g` is the acceleration due to gravity, and `\nu` is
:config:`ocean.plume.power_alpha` (default `1/3`, the published value). `\dot m_{fw} = 0`
where `q_{sg} = 0` or the plume is not buoyant.

The ambient and discharge contributions are combined (Eqn. 16 in :cite:`Pelle2023`) as

.. math::
   :label: eq-plume-melt-total

   \dot m = \frac{\hat M(\hat X)\, M^2}{M + \dot m_{fw}} + \dot m_{fw},

which reduces to `\dot m_a` if there is no discharge and to `\dot m_{fw}` if `\hat M M
\ll \dot m_{fw}`. The sub-shelf mass flux passed to the ice dynamics core is
`\rho_i \dot m`.

.. rubric:: Where the discharge comes from: hydrology routing

PISM models subglacial water under grounded ice only, so the discharge enters the cavity
at the grounding line. The plume model uses the magnitude of the subglacial water flux
`\mathbf{q}` provided by the hydrology model (``-hydrology routing``; see
:ref:`sec-hydrology-routing`, Eq. :eq:`eq-flux`), which is already a flux per unit width in
`\mathrm{m^2\,s^{-1}}`, so no conversion using a channel width is needed.

Every grounded cell that is next to floating ice and has a non-zero flux is an *outflow*
with `q_{sg,0} = |\mathbf{q}|`. For each outflow `\alpha`, `z_{\mathrm{gl}}`, `T_a` and
`S_a` are averaged over floating cells within 5 km, and the *governing length scale*

.. math::
   :label: eq-plume-governing-length

   R = 5 L' = \frac{5\, q_{sg,0}}{\dot m_{fw,0}}

is computed from the discharge melt rate `\dot m_{fw,0}` of Eq.
:eq:`eq-plume-fresh-water-melt` (Eqn. 15 in :cite:`Pelle2023`). The discharge flux on the
floating cells (the field `q_{sg}(x,y)` used above) decays quadratically with the distance
`d` from the outflow:

.. math::
   :label: eq-plume-discharge-decay

   q_{sg} = q_{sg,0}\left(1 - \frac{d}{R}\right)^2,\quad d < R,

and is zero beyond `R`. Cells closer than `1.5\,\max(\Delta x, \Delta y)` to an outflow
receive the full `q_{sg,0}`, so a plume that is smaller than a grid cell is not lost.
:config:`ocean.plume.discharge_method` selects how `d` is measured and how overlapping
plumes are treated:

``along_flow`` (default)
   `q_{sg,0}`, `R` and the along-flow path length `d` from the outflow are transported
   downstream with the ice flow using the same semi-Lagrangian scheme as `z_{\mathrm{gl}}`.
   Each floating cell therefore sees a single upstream outflow.

``isotropic``
   `d` is the Euclidean distance from the outflow; the maximum is taken where plumes
   overlap.

``downstream_gate``
   As ``isotropic``, but only cells downstream of the outflow (positive dot product of
   the outflow-to-cell vector with the ice flow direction) are affected.

.. rubric:: Parameters

Prefix: ``ocean.plume.``

.. pism-parameters::
   :prefix: ocean.plume.

.. rubric:: Diagnostics

The fields `\dot m` (``plume_basal_melt_rate``), `\dot m_{fw}`
(``plume_fresh_water_melt_rate``), `q_{sg}` (``plume_discharge_flux``), `z_{\mathrm{gl}}`
(``plume_grounding_line_elevation``), `\alpha` (``plume_local_slope``), `z_{\mathrm{b}}`
(``plume_shelf_base_elevation``), `T_a` (``plume_temperature``) and `S_a`
(``plume_salinity``) are available as spatial diagnostics. They are reported on floating
ice only; elsewhere they contain :config:`output.fill_value`.

.. _sec-picop:

PICOP
+++++

:|options|: ``-ocean picop``
:|variables|: :var:`theta_ocean` (potential ocean temperature), [kelvin],

              :var:`salinity_ocean` (salinity of the adjacent ocean), [g/kg],

              :var:`basins` (mask of large-scale ocean basins that ocean input is averaged over), [integer]
:|implementation|: ``pism::ocean::Picop``
:|seealso|: :ref:`sec-pico`, :ref:`sec-plume`

PICOP :cite:`Pelle2019` combines the PICO box model (:ref:`sec-pico`) with the buoyant
plume melt rate parameterization (:ref:`sec-plume`). It addresses the main limitation of
PICO: the box model produces a smooth melt pattern that depends on the box a cell
belongs to, while observations show melt rates that are largest near the grounding line
and strongly controlled by the local basal slope. The plume parameterization captures
this, but needs ambient conditions at every floating cell.

PICOP provides them by running PICO first:

#. PICO averages the forcing (:var:`theta_ocean`, :var:`salinity_ocean`) on the
   continental shelf in each basin and solves the box model for each ice shelf, as
   described in :ref:`sec-pico`. This yields, in every floating cell, the temperature
   `T_a` and salinity `S_a` of the water in the box that contains the cell (the PICO
   fields :var:`pico_temperature` and :var:`pico_salinity`). Water that has been cooled
   and freshened by melt upstream is therefore colder and fresher downstream.
#. These fields are used as the ambient conditions `T_a` and `S_a` of the plume model,
   which computes the melt rate `\dot m` using the equations in :ref:`sec-plume`:
   transport of the grounding line elevation, the basal slope, ambient melt, and
   (if enabled) subglacial discharge.
#. The sub-shelf ice temperature is the one computed by PICO.

The PICO melt rate is *replaced* by the plume melt rate on floating ice; PICO is
used for the ambient ocean state only. In regions where PICO's geometry cannot be
identified, PICO falls back to :cite:`BeckmannGoosse2003` and these values are used as
ambient conditions.

Unlike ``-ocean plume``, PICOP interprets :var:`theta_ocean` as potential temperature (as
PICO does) and
:config:`ocean.plume.temperature_as_thermal_forcing` has no effect. The forcing is
read from :config:`ocean.pico.file`.

Like the plume model, PICOP requires the stress balance model. If it is not available
(for example during bootstrapping), PICOP falls back to the PICO melt rates. The
subglacial discharge needs a hydrology model that provides the water flux (see
:ref:`sec-hydrology-routing`).

PICOP has no parameters of its own: the box model is configured with ``ocean.pico.*``
(see :ref:`sec-pico`) and the plume with ``ocean.plume.*`` (see :ref:`sec-plume`). The
``plume_*`` diagnostics listed in :ref:`sec-plume` are available; `T_a` and `S_a` are the
PICO fields.

.. _sec-ocean-delta-sl:

Scalar sea level offsets
++++++++++++++++++++++++

:|options|: :opt:`-sea_level ...,delta_sl`
:|variables|: :var:`delta_SL` (meters)
:|implementation|: ``pism::ocean::sea_level::Delta_SL``

The ``delta_sl`` modifier implements sea level forcing using scalar offsets.

.. rubric:: Parameters

Prefix: ``ocean.delta_sl.``

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

.. rubric:: Parameters

Prefix: ``ocean.delta_sl_2d.``

.. pism-parameters::
   :prefix: ocean.delta_sl_2d.

.. _sec-ocean-delta-t:

Scalar sub-shelf temperature offsets
++++++++++++++++++++++++++++++++++++


:|options|: :opt:`-ocean ...,delta_T`
:|variables|: :var:`delta_T` (kelvin)
:|implementation|: ``pism::ocean::Delta_T``

This modifier implements forcing using sub-shelf ice temperature offsets.

.. rubric:: Parameters

Prefix: ``ocean.delta_T.``

.. pism-parameters::
   :prefix: ocean.delta_T.

.. _sec-ocean-delta-smb:

Scalar sub-shelf mass flux offsets
++++++++++++++++++++++++++++++++++

:|options|: ``-ocean ...,delta_SMB``
:|variables|: :var:`delta_SMB` |flux|
:|implementation|: ``pism::ocean::Delta_SMB``

This modifier implements forcing using sub-shelf mass flux (melt rate) offsets.

.. rubric:: Parameters

Prefix: ``ocean.delta_mass_flux.``

.. pism-parameters::
   :prefix: ocean.delta_mass_flux.

.. _sec-ocean-frac-smb:

Scalar sub-shelf mass flux fraction offsets
+++++++++++++++++++++++++++++++++++++++++++

:|options|: ``-ocean ...,frac_SMB``
:|variables|: :var:`frac_SMB` [1]
:|implementation|: ``pism::ocean::Frac_SMB``

This modifier implements forcing using sub-shelf mass flux (melt rate) fraction offsets.

.. rubric:: Parameters

Prefix: ``ocean.frac_mass_flux.``

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

See also to ``-atmosphere ...,anomaly`` or ``-surface ...,anomaly`` (section
:ref:`sec-surface-anomaly`) which is similar, but applies anomalies at the atmosphere or
surface level, respectively.

.. rubric:: Parameters

Prefix: ``ocean.anomaly.``

.. pism-parameters::
   :prefix: ocean.anomaly.

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

.. rubric:: Parameters

Prefix: ``ocean.delta_MBP.``

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

.. rubric:: Parameters

Prefix: ``ocean.frac_MBP.``

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

.. rubric:: Parameters

Prefix: ``ocean.cache.``

.. pism-parameters::
   :prefix: ocean.cache.
