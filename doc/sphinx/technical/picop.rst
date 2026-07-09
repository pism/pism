.. include:: ../global.txt

.. math::

.. _sec-technical-picop:

PICOP with subglacial discharge: technical details
==================================================

.. contents::

.. _sec-picop-intro:

Introduction
------------

PICOP (:cite:`Pelle2019`) combines the box model PICO (:cite:`ReeseAlbrecht2018`)
with the buoyant-meltwater-plume parameterization of :cite:`Lazeroms2018` to compute
sub–ice-shelf basal melt rates. PICO supplies the ambient ocean temperature and salinity
in the sub-shelf cavity; the plume model turns those, together with the local ice-shelf
draft geometry, into a melt rate.

PISM additionally implements the **subglacial-discharge** extension of
:cite:`Pelle2023`, in which fresh subglacial water crossing the grounding line provides
an extra source of buoyancy that locally enhances melt within a distance `5L'` of the
outflow. This note documents the equations and their implementation.

Select the model with ``-ocean.models picop``. The implementation lives in
``src/coupler/ocean/Picop.{cc,hh}`` (coupling, geometry, discharge field) and
``src/coupler/ocean/PicopPhysics.{cc,hh}`` (the plume equations). Equation numbers below
refer to :cite:`Pelle2023`.

.. note::

   Unlike the reference ISSM implementation, which works in m/year, PISM works in SI
   units throughout. Melt rates are therefore in `\text{m}\,\text{s}^{-1}` and the
   paper's `31536000` (m/s → m/year) conversion factor is omitted.

.. _sec-picop-plume:

Meltwater plume melt rate
-------------------------

The plume originates at the grounding line and travels up the ice-shelf base, driven by
the buoyancy of the meltwater it produces. Its behavior is controlled by the ambient
temperature and salinity `T_a`, `S_a`, the local basal slope angle `\alpha`, and the
depth `z_{gl}` of the grounding line that feeds a given shelf point.

The depth-dependent freezing point enters through

.. math::

   \Delta T_{f,gl} = \frac{T_a - \lambda_1 S_a + \lambda_2 + \lambda_3 z_{gl}}{\lambda_3},

and the *thermal forcing* that drives melt is `\lambda_3 \Delta T_{f,gl}`, which in the
code is formed as `T_a - T_{f,gl}` where `T_{f,gl} = \lambda_1 S_a + \lambda_2 +
\lambda_3 z_{gl}` (``characteristic_freezing_point``). Dimensionless geometric factors
`G_1, G_2, G_3` (``geometric_scaling``), the melt scale `M` (``melt_function``,
`M = M_0\,G(\alpha)\,(\lambda_3\Delta T_{f,gl})^2`), the dimensionless along-plume
coordinate `\hat X \in [0,1]` (``dimensionless_coordinate``), and the universal melt
curve `M(\hat X)` (``dimensionless_melt_curve``) combine into the plume melt rate
`m_o = M\,M(\hat X)` (their Eq. 12). See :cite:`Lazeroms2018` for the derivation.

.. _sec-picop-discharge:

Subglacial-discharge enhancement
--------------------------------

Because the plume melt rate vanishes at the grounding line (no buoyancy source there),
:cite:`Pelle2023` add a discharge-driven melt rate derived from :cite:`Jenkins2011`.
With the subglacial discharge flux `q_{sg}` (units `\text{m}^2\,\text{s}^{-1}`) and the
density contrast `\Delta\rho_i = \beta_S S_a - \beta_T\,(\lambda_3\Delta T_{f,gl})`
(assuming fresh discharge, `S = 0`), the discharge melt rate is (their Eq. 13, in SI):

.. math::
   :label: eq-picop-mfw

   m_{fw} = \frac{c_p}{L_{fw}}\; C_d^{1/2}\Gamma_{Ts0}\; G_1\; G_2^{1/3}\;
            \bigl(g\, q_{sg}\, \Delta\rho_i\bigr)^{1/3}\;
            (\lambda_3\Delta T_{f,gl}).

Two subtleties in :eq:`eq-picop-mfw` that are easy to get wrong (``fresh_water_melt_rate``):

* the prefactor is `c_p / L_{fw} = 1/(L_{fw}/c_p) \approx 0.01\,\text{K}^{-1}` --- **not**
  `1/(L_{fw}\,c_p)`; and
* both `\Delta\rho_i` and the trailing factor use the **thermal forcing**
  `\lambda_3\Delta T_{f,gl} = T_a - T_{f,gl}`, not the bare freezing temperature.

The `(g\,q_{sg})^{1/3}` term makes `m_{fw}` scale with the one-third power of discharge,
as expected for a horizontal plume.

.. _sec-picop-lengthscale:

Governing length scale
----------------------

The discharge enhancement applies only within a radius `5L'` of each outflow, where
plume buoyancy is discharge-dominated. Rather than evaluate the paper's explicit Eq. 15,
PISM uses the fact that Eq. 13 and Eq. 15 share the same geometric/thermal block and are
reciprocal in it, so they collapse to (``governing_length_scale``)

.. math::
   :label: eq-picop-length

   5L' = \frac{5\, q_{sg}}{m_{fw}}.

This is algebraically identical to Eq. 15: the paper's leading constant `500 = 5\,
L_{fw}/c_p` is carried by `1/m_{fw}` (whose prefactor is `c_p/L_{fw}`), so no magic
number appears. Because `G_1` and `G_2` (and hence the slope dependence) cancel in the
product `m_{fw}\cdot L'`, :eq:`eq-picop-length` retains the **full** slope dependence and
does **not** rely on the `\sin\alpha \approx 0.01` simplification used to derive the
closed form of Eq. 15 for typical Antarctic shelves. Physically, `5L'` is the distance
over which the accumulated melt `m_{fw}\times\text{distance}` balances the input flux.

.. _sec-picop-field:

Constructing the discharge field `q_{sg}(x,y)`
----------------------------------------------

PISM routes subglacial water under grounded ice only, so discharge enters the cavity at
the grounding line and must be spread onto floating cells. The field is rebuilt every
ocean time step (the grounding line moves) in ``Picop::build_discharge_field``:

#. **Outflows.** A grounding-line outflow is a grounded-ice cell adjacent to floating ice
   (``next_to_floating_ice``) with nonzero water flux. Its `q_{sg,0}` is the magnitude of
   ``hydrology->flux()`` there. PISM's hydrology flux is already a flux per unit width
   (`\text{m}^2\,\text{s}^{-1}`), so --- unlike ISSM --- no channel-width conversion is
   needed.

#. **5-km averaging.** Following :cite:`Pelle2023`, `\alpha`, `z_{gl}`, `T_a`, and `S_a`
   are averaged over floating cells within 5 km of each outflow before `m_{fw,0}` and
   `5L'` are computed. The accumulation is summed across ranks with ``GlobalSum`` so the
   average spans MPI subdomain boundaries.

#. **Painting.** For each floating cell within `5L'` of an outflow,

   .. math::

      q_{sg}(x,y) = q_{sg,0}\left(1 - \frac{d}{5L'}\right)^2,

   with `d` the distance to the outflow; where radii overlap the maximum is taken
   (the paper's ``maskdist`` `\times` ``maskdis``). Cells beyond every `5L'` get
   `q_{sg} = 0`.

The final melt rate blends the plume and discharge contributions (their Eq. 16):

.. math::

   m = \frac{M(\hat X)\, M^2}{M + m_{fw}} + m_{fw},

which reduces to the plume melt `m_o` where `m_{fw} = 0`.

.. _sec-picop-parallel:

Parallel implementation
-----------------------

The search radius `5L'` (order 10–20 km for Denman-scale fluxes) spans several grid
cells and can cross MPI subdomain boundaries. To avoid wide-stencil communication, the
outflow list --- only `O(\text{grounding-line length}/\Delta x)` points --- is gathered
onto every rank with ``MPI_Allgatherv`` (helper ``gather_outflows``). Each rank then
paints its own floating cells against the global outflow list. The per-outflow 5-km
averages are formed with a single ``GlobalSum`` over a packed accumulator buffer.

.. _sec-picop-io:

Configuration and diagnostics
-----------------------------

The plume constants are configuration parameters under ``ocean.picop.*`` (e.g.
``ocean.picop.entrainment_coefficient``, ``ocean.picop.heat_exchange_parameter``,
``ocean.picop.drag_coefficient``); PICO parameters remain under ``ocean.pico.*``.

Spatial diagnostics:

* ``picop_basal_melt_rate`` --- total melt rate `m` (`\text{m}\,\text{s}^{-1}`);
* ``picop_fresh_water_melt_rate`` --- discharge contribution `m_{fw}`;
* ``picop_discharge_flux`` --- the constructed `q_{sg}(x,y)`
  (`\text{m}^2\,\text{s}^{-1}`);
* ``picop_grounding_line_elevation``, ``picop_local_slope``.

As a validation check, `q_{sg}` should be nonzero in a band `\sim 5L'` wide adjacent to
the grounding line with a peak near the grounding-line water flux, and `m_{fw}` should
reach tens of m/year near the grounding line (:cite:`Pelle2023`, Fig. 2).
