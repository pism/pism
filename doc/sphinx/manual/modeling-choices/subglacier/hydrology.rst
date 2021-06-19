.. include:: ../../../global.txt

.. _sec-subhydro:

Subglacial hydrology
--------------------

At the present time, two simple subglacial hydrology models are implemented *and
documented* in PISM, namely ``-hydrology null`` (and its modification ``steady``) and
``-hydrology routing``; see :numref:`tab-hydrologychoice` and :cite:`BuelervanPelt2015`.
In both models, some of the water in the subglacial layer is stored locally in a layer of
subglacial till by the hydrology model. In the ``routing`` model water is conserved by
horizontally-transporting the excess water (namely ``bwat``) according to the gradient of
the modeled hydraulic potential. In both hydrology models a state variable ``tillwat`` is
the effective thickness of the layer of liquid water in the till; it is used to compute
the effective pressure on the till (see :ref:`sec-basestrength`). The pressure of the
transportable water ``bwat`` in the ``routing`` model does not relate directly to the
effective pressure on the till.

.. note::

   Both ``null`` and ``routing`` described here provide all diagnostic quantities needed
   for mass accounting, even though ``null`` is not mass-conserving. See
   :ref:`sec-mass-conservation-hydrology` for details.

.. list-table:: Command-line options to choose the hydrology model
   :name: tab-hydrologychoice
   :header-rows: 1
   :widths: 2,5

   * - Option
     - Description
   * - :opt:`-hydrology null`
     - The default model with only a layer of water stored in till. Not mass conserving in
       the map-plane but much faster than ``-hydrology routing``. Based on "undrained
       plastic bed" model of :cite:`Tulaczyketal2000b`. The only state variable is
       ``tillwat``.

   * - :opt:`-hydrology steady`
     - A version of the ``null`` model that computes an approximation of the steady state
       subglacial water flux that can be used as an input for a frontal melt
       parameterization such as the one described in :ref:`sec-ismip6-frontal-melt`.

   * - :opt:`-hydrology routing`
     - A mass-conserving horizontal transport model in which the pressure of transportable
       water is equal to overburden pressure. The till layer remains in the model, so this
       is a "drained and conserved plastic bed" model. The state variables are ``bwat``
       and ``tillwat``.

See :numref:`tab-hydrology` for options which apply to all hydrology models. Note
that the primary water source for these models is the energy conservation model which
computes the basal melt rate ``basal_melt_rate_grounded``. There is, however, also option
:opt:`-hydrology_input_to_bed_file` which allows the user to *add* water directly into the
subglacial layer, in addition to the computed ``basal_melt_rate_grounded`` values. Thus
``-hydrology_input_to_bed_file`` allows the user to model drainage directly to the bed
from surface runoff, for example. Also option :opt:`-hydrology_bmelt_file` allows the user
to replace the computed ``basal_melt_rate_grounded`` rate by values read from a file,
thereby effectively decoupling the hydrology model from the ice dynamics
(esp. conservation of energy).

To use water runoff computed by a PISM surface model, set
:config:`hydrology.surface_input_from_runoff`. (The :ref:`sec-surface-pdd` computes runoff
using computed melt and the assumption that a fraction of this melt re-freezes, all other
models assume that all melt turns into runoff.)

.. list-table:: Subglacial hydrology command-line options which apply to all hydrology models
   :name: tab-hydrology
   :header-rows: 1
   :widths: 1,1

   * - Option
     - Description
   * - :opt:`-hydrology.surface_input.file`
     - Specifies a NetCDF file which contains a time-dependent field ``water_input_rate``
       which has units of water thickness per time. This rate is *added to* the basal melt
       rate computed by the energy conservation code.
   * - :opt:`-hydrology_tillwat_max` (m)
     - Maximum effective thickness for water stored in till.
   * - :opt:`-hydrology_tillwat_decay_rate` (m/a)
     - Water accumulates in the till at the basal melt rate ``basal_melt_rate_grounded``,
       minus this rate.
   * - :opt:`-hydrology_use_const_bmelt`
     - Replace the conservation-of-energy basal melt rate ``basal_melt_rate_grounded``
       with a constant.

.. _sec-hydrology-null:

The default model: ``-hydrology null``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this model the water is *not* conserved but it is stored locally in the till up to a
specified amount; the configuration parameter :config:`hydrology.tillwat_max` sets that
amount. The water is not conserved in the sense that water above the
:config:`hydrology.tillwat_max` level is lost permanently. This model is based on the
"undrained plastic bed" concept of :cite:`Tulaczyketal2000b`; see also
:cite:`BBssasliding`.

In particular, denoting ``tillwat`` by `W_{till}`, the till-stored water layer effective
thickness evolves by the simple equation

.. math::
   :label: eq-tillwatevolve

   \frac{\partial W_{till}}{\partial t} = \frac{m}{\rho_w} - C

where `m=` :var:`basal_melt_rate_grounded` (kg `\text{m}^{-2}\,\text{s}^{-1}`), `\rho_w`
is the density of fresh water, and `C` :var:`hydrology_tillwat_decay_rate`. At all times
bounds `0 \le W_{till} \le W_{till}^{max}` are satisfied.

This ``-hydrology null`` model has been extensively tested in combination with the
Mohr-Coulomb till (section :ref:`sec-basestrength`) for modelling ice streaming (see
:cite:`AschwandenAdalgeirsdottirKhroulev` and :cite:`BBssasliding`, among others).

.. _sec-hydrology-steady:

"Steady state flow" model
~~~~~~~~~~~~~~~~~~~~~~~~~

This model (``-hydrology steady``) adds an approximation of the subglacial water flux to
the ``null`` model. It *does not* model the steady state distribution of the subglacial
water thickness.

Here we assume that the water input from the surface read from the file specified using
:config:`hydrology.surface_input.file` instantaneously percolates to the base of the ice
and enters the subglacial water system.

We also assume that the subglacial drainage system instantaneously reaches its steady
state. Under this assumption the water flux can be approximated by using an auxiliary
problem that is computationally cheaper to solve, compared to using the mass conserving
``routing`` model (at least at high grid resolutions).

We solve

.. math::
   :label: eq-steady-hydro-aux

   \diff{u}{t} = -\Div (\V u)

on the grounded part of the domain with the initial state `u_0 = \tau M`, where `\tau` is
the scaling of the input rate (:config:`hydrology.steady.input_rate_scaling`) and `M` is
the input rate read from an input file.

The velocity field `\V` approximates the steady state flow direction. In this
implementation `\V` is the normalized gradient of

.. math::

   \psi = \rho_i g H + \rho_w g b + \Delta \psi.

Here `H` is the ice thickness, `b` is the bed elevation, `g` is the acceleration due to
gravity and `\rho_i`, `\rho_w` are densities of ice and water, respectively.

.. note::

   Note that `\psi` ignores subglacial water thickness, producing flow patterns that are
   more concentrated than would be seen in a model that includes it. This effect is
   grid-dependent.

   This implies that this code cannot model the distribution of the flow along the
   grounding line of a particular outlet glacier but it *can* tell us how surface runoff
   is distributed among the outlets.

The term `\Delta \psi` is the adjustment needed to remove internal minima from the "raw"
potential, filling any "lakes" it might have. This modification of `\psi` is performed
using an iterative process; set :config:`hydrology.steady.potential_n_iterations` to
control the maximum number of iterations and :config:`hydrology.steady.potential_delta` to
affect the amount it is adjusted by at each iteration.

The equation :eq:`eq-steady-hydro-aux` is advanced forward in time until `\int_{\Omega}u
< \epsilon\int_{\Omega} u_0`, where `\epsilon` (:config:`hydrology.steady.volume_ratio`)
is a fraction controlling the accuracy of the approximation: lower values of `\epsilon`
would result in a better approximation and a higher computational cost (more iterations).
Set :config:`hydrology.steady.n_iterations` to control the maximum number of these
iterations.

This model restricts the time step length in order to capture the temporal variability of
the forcing: the flux is updated at least once for each time interval in the forcing file.

In simulations with a slowly-varying surface input rate the flux has to be updated every
once in a while to reflect changes in the flow pattern coming from changing geometry. Use
:config:`hydrology.steady.flux_update_interval` (years) to set the update frequency.

See :ref:`sec-steady-hydro` for technical details.

.. _sec-hydrology-routing:

The mass-conserving model: ``-hydrology routing``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this model the water *is* conserved in the map-plane. Water does get put into the till,
with the same maximum value :config:`hydrology.tillwat_max`, but excess water is
horizontally-transported. An additional state variable ``bwat``, the effective thickness
of the layer of transportable water, is used by ``routing``. This transportable water will
flow in the direction of the negative of the gradient of the modeled hydraulic potential.
In the ``routing`` model this potential is calculated by assuming that the transportable
subglacial water is at the overburden pressure :cite:`Siegertetal2007`. Ultimately the
transportable water will reach the ice sheet grounding line or ice-free-land margin, at
which point it will be lost. The amount that is lost this way is reported to the user.

In this model ``tillwat`` also evolves by equation :eq:`eq-tillwatevolve`, but several
additional parameters are used in determining how the transportable water ``bwat`` flows
in the model; see :numref:`tab-hydrologyrouting`. Specifically, the horizontal
subglacial water flux is determined by a generalized Darcy flux relation :cite:`Clarke05`,
:cite:`Schoofetal2012`

.. math::
   :label: eq-flux

   \bq = - k\, W^\alpha\, |\nabla \psi|^{\beta-2} \nabla \psi

where `\bq` is the lateral water flux, `W=` ``bwat`` is the effective thickness of the
layer of transportable water, `\psi` is the hydraulic potential, and `k`, `\alpha`,
`\beta` are controllable parameters (:numref:`tab-hydrologyrouting`).

In the ``routing`` model the hydraulic potential `\psi` is determined by

.. math::
   :label: eq-hydraulicpotential

   \psi = P_o + \rho_w g (b + W)

where `P_o=\rho_i g H` is the ice overburden pressure, `g` is gravity, `\rho_i` is ice
density, `\rho_w` is fresh water density, `H` is ice thickness, and `b` is the bedrock
elevation.

For most choices of the relevant parameters and most grid spacings, the ``routing`` model
is at least two orders of magnitude more expensive computationally than the ``null``
model. This follows directly from the CFL-type time-step restriction on lateral flow of a
fluid with velocity on the order of centimeters to meters per second (i.e. the subglacial
liquid water ``bwat``). (By comparison, much of PISM ice dynamics time-stepping is
controlled by the much slower velocity of the flowing ice.) Therefore the user should
start with short runs of order a few model years. We also recommend ``daily`` or even
``hourly`` reporting for scalar and spatially-distributed time-series to see hydrology
model behavior, especially on fine grids (e.g. `< 1` km).

.. list-table:: Command-line options specific to hydrology model ``routing``
   :name: tab-hydrologyrouting
   :header-rows: 1
   :widths: 5,4

   * - Option
     - Description
   * - :opt:`-hydrology_hydraulic_conductivity` `k`
     - `=k` in formula :eq:`eq-flux`.
   * - :opt:`-hydrology_null_strip` (km)
     - In the boundary strip water is removed and this is reported. This option specifies
       the width of this strip, which should typically be one or two grid cells.
   * - :opt:`-hydrology_gradient_power_in_flux` `\beta`
     - `=\beta` in formula :eq:`eq-flux`.
   * - :opt:`-hydrology_thickness_power_in_flux` `\alpha`
     - `=\alpha` in formula :eq:`eq-flux`.

.. FIXME -hydrology distributed is not documented except by :cite:`BuelervanPelt2015`
