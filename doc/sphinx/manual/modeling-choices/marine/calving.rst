.. include:: ../../../global.txt

.. _sec-calving:

Calving and front retreat
-------------------------

Overview
^^^^^^^^

All mechanisms described below fall into two categories:

- mechanisms computing a *retreat rate* due to calving and using it to update ice geometry
  (:ref:`sec-calving-eigen-calving`, :ref:`sec-calving-vonmises`,
  :ref:`sec-calving-hayhurst`, :ref:`sec-calving-cliff-shear`, :ref:`sec-calving-cliff-tensile`), and
- mechanisms removing ice at a grid point according to a certain criterion
  (:ref:`sec-calving-thickness-threshold`, :ref:`sec-calving-floating-ice`,
  :ref:`sec-prescribed-retreat`).

To select several calving mechanisms, use a comma-separated list of corresponding
keywords. For example,

.. code-block:: none

   -calving eigen_calving,thickness_calving

selects :ref:`sec-calving-eigen-calving` and :ref:`sec-calving-thickness-threshold`.

If more than one retreat-rate-based mechanism is selected, the corresponded rates are
*added* and then used to update ice extent. (In other words: selected calving mechanisms
are applied *together* instead of applying their effects *in sequence*.)

The partially-filled grid cell formulation (section :ref:`sec-part-grid`) provides a
framework suitable to relate calving rates to the mass transport scheme at the ice shelf
terminus. Ice shelf front advance and retreat due to calving are limited to a maximum of
one grid cell length per (adaptive) time step. The combined calving rate (velocity) can be
used to limit the overall timestep of PISM (thus slowing down all of PISM) by using
:config:`geometry.front_retreat.use_cfl`. This "CFL-type" time-step limitation is
definitely recommended in high-resolution runs which attempt to model calving position
accurately. Without this option, under certain conditions where PISM's adaptive time step
happens to be long enough, dendritic structures can appear at the calving front because
the calving mechanism cannot "keep up" with the computed calving rate.

Setting the flag :config:`geometry.front_retreat.wrap_around` to ``true`` allows the front
retreat to "wrap around" the computational domain. (This is appropriate in some regional
synthetic geometry setups.)

.. _sec-calving-rate-modulation:

Scaling calving rates
=====================

Set :config:`calving.rate_scaling.file` to scale the *total* (combined) calving rate from
all selected rate-based mechanisms, e.g. to introduce calving variability corresponding to
seasonal changes in ice melange. The file used with this option should contain the scalar
time-dependent variable :var:`frac_calving_rate` (units: `1`).

.. rubric:: Parameters

Prefix: ``calving.rate_scaling.``

.. pism-parameters::
   :prefix: calving.rate_scaling.

.. _sec-calving-eigen-calving:

Eigen calving
^^^^^^^^^^^^^

PISM-PIK introduced a physically-based 2D-calving parameterization
:cite:`Levermannetal2012`. This calving parameterization is turned on in PISM by option
:opt:`-calving eigen_calving`. Average calving rates, `c`, are proportional to the product
of principal components of the horizontal strain rates, `\dot{\epsilon}_{_\pm}`, derived
from SSA-velocities

.. math::
   :label: eq-calv2

   c &= K\; \dot{\epsilon}_{_+}\; \dot{\epsilon}_{_-},

   \dot{\epsilon}_{_\pm} &> 0.

The rate `c` is in `\text{m}\,\text{s}^{-1}`, and the principal strain rates
`\dot\epsilon_\pm` have units `\text{s}^{-1}`, so `K` has units `\text{m}\,\text{s}`. The
constant `K` incorporates material properties of the ice at the front. It can be set using
:config:`calving.eigen_calving.K`.

The actual strain rate pattern strongly depends on the geometry and boundary conditions
along the confinements of an ice shelf (coast, ice rises, front position). The strain rate
pattern provides information in which regions preexisting fractures are likely to
propagate, forming rifts (in two directions). These rifts may ultimately intersect,
leading to the release of icebergs. This (and other) ice shelf calving models are not
intended to resolve individual rifts or calving events, but it produces
structurally-stable calving front positions which agree well with observations. Calving
rates balance calving-front ice flow velocities on average.

.. rubric:: Parameters

Prefix: ``calving.eigen_calving.``

.. pism-parameters::
   :prefix: calving.eigen_calving.

.. _sec-calving-vonmises:

von Mises stress calving
^^^^^^^^^^^^^^^^^^^^^^^^

While ``eigen_calving`` (section :ref:`sec-calving-eigen-calving`) is appropriate for
Antartic ice shelves, it does not work for outlet glaciers that flow in narrow fjords.
Along valleys with nearly parallel walls the transverse component of the velocity is
close to zero and the transversal strain rate is therefore also close to zero and noisy.

Instead of the product of the eigen strain rates, :cite:`Morlighem2016` propose a calving
law where the calving rate `c` is functionally related to tensile stresses:

.. math::
   :label: eq-calv3

   c = |\mathbf{u}| \frac{\tilde{\sigma}}{\sigma_{max}},

where `\tilde{\sigma}` is the tensile von Mises stress and `\sigma_{max}` is a threshold
that has units `Pa` (see :config:`calving.vonmises_calving.sigma_max`). As the tensile
fracture strength is much smaller than the compressive fracture strength, the effective
tensile strain rate is defined as

.. math::
   :label: eq-calv4

   \tilde{\dot{\epsilon}}_e = \left(\frac{1}{2}\left(\max(0,\dot{\epsilon}_{_+})^2 +
   \max(0,\dot{\epsilon}_{_-})^2\right)\right)^{1/2}.

Following :cite:`Morlighem2016`, `\tilde{\sigma}` is given by

.. math::
   :label: eq-calv5

   \tilde{\sigma} = \sqrt{3} B \tilde{\dot{{\epsilon}}}_e^{1/n},

where `B` is the ice hardness.

.. rubric:: Parameters

Prefix: ``calving.vonmises_calving.``

.. pism-parameters::
   :prefix: calving.vonmises_calving.

.. _sec-calving-hayhurst:

Hayhurst calving
^^^^^^^^^^^^^^^^

The option :opt:`-calving hayhurst_calving` implements the parameterization described in
:cite:`Mercenier2018` (equation 22). This mechanism models calving based on the Hayhurst
creep damage parameter, which accounts for both tensile and shear stresses.

The calving rate is calculated based on the maximum tensile stress at the ice front and
takes into account:

- The water depth relative to ice thickness
- A damage threshold stress
- A damage rate parameter
- A damage law exponent

The calving rate is computed as:

.. math::
   :label: eq-hayhurst

   c = |\mathbf{u}| \frac{\tilde{B} (\sigma_0 - \sigma_{th})^r}{\sigma_{th}}

where:
- `c` is the calving rate
- `|\mathbf{u}|` is the ice velocity magnitude
- `\tilde{B}` is the effective damage rate parameter (default: 65.0 MPa^r/year)
- `\sigma_0` is the maximum tensile stress
- `\sigma_{th}` is the damage threshold stress (default: 0.17 MPa)
- `r` is the damage law exponent (default: 0.43)

The maximum tensile stress `\sigma_0` is approximated as:

.. math::
   :label: eq-hayhurst-sigma

   \sigma_0 = (0.4 - 0.45 (\omega - 0.065)^2) \rho_i g H

where:
- `\omega` is the ratio of water depth to ice thickness
- `\rho_i` is the ice density
- `g` is the gravitational acceleration
- `H` is the ice thickness

.. note::

   The calving parameterization was derived for grounded ice but is also applied to floating ice cells
   (where `\omega > \rho_i/\rho_w`). Without this modification, the formation of an ice shelf would effectively
   stop calving. To prevent this, the ice thickness `H` is adjusted to include the freeboard height and
   `\omega` is recalculated, ensuring continued calving at the grounding line.

.. rubric:: Parameters

Prefix: ``calving.hayhurst_calving.``

.. pism-parameters::
   :prefix: calving.hayhurst_calving.

.. _sec-calving-cliff-shear:

Cliff calving due to shear failure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The option :opt:`-calving cliff_calving_shear` implements a parameterization for calving of ice cliffs
based on shear stress failure :cite:`Schlemm2019`. This mechanism is designed to model calving of
grounded ice cliffs that form when ice shelves collapse, exposing grounded ice to the ocean.

The calving rate is calculated based on the cliff height and water depth, with a maximum rate
limited by mélange buttressing :cite:`Schlemm2021`. The parameterization takes into account:

- The height of the ice cliff above water level
- The ratio of water depth to ice thickness
- A maximum calving rate due to mélange buttressing

The calving rate is computed as:

.. math::
   :label: eq-cliff-shear

   c = \frac{C_0 (F-F_c)^s}{F_s (1 + C_0 (F-F_c)^s / c_{max})}

where:
- `c` is the calving rate
- `C_0` is a scaling factor (default: 90 m/year)
- `F` is the cliff height above water level
- `F_c` is a critical cliff height threshold
- `F_s` is a scaling factor
- `s` is an exponent
- `c_{max}` is the maximum calving rate due to mélange buttressing (default: 3000 m/year)

.. note::

   The maximum calving rate parameter (`c_{max}`) controls the effect of mélange buttressing.
   Setting this parameter to a very large value (e.g., 100 km/year) effectively disables
   mélange buttressing, resulting in nearly unbuttressed cliff calving rates. However, this should
   be used with caution as it can lead to unrealistically high calving rates :cite:`Schlemm2021`.

.. rubric:: Parameters

Prefix: ``calving.cliff_calving_shear.``

.. pism-parameters::
   :prefix: calving.cliff_calving_shear.

.. _sec-calving-cliff-tensile:

Cliff calving due to tensile failure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The option :opt:`-calving cliff_calving_tensile` implements a parameterization for calving of ice cliffs
based on tensile stress failure :cite:`Crawford2021`. This mechanism models calving of grounded ice
cliffs that form when ice shelves collapse, with a focus on tensile failure at the cliff face.

The calving rate is calculated based on the cliff height, with a minimum threshold height required
for calving to occur. The parameterization takes into account:

- The height of the ice cliff above water level
- A minimum cliff height threshold (135 m)
- A power law relationship between cliff height and calving rate

The calving rate is computed as:

.. math::
   :label: eq-cliff-tensile

   c = I H_c^\alpha

where:
- `c` is the calving rate
- `I` is a scaling factor (default: 3.7e-16 m/s)
- `H_c` is the cliff height above water level
- `\alpha` is an exponent (default: 6.9)

Calving only occurs when the cliff height exceeds 135 meters.

.. note::

   The parameters `I` and `\alpha` depend on ice temperature and bed conditions (frozen, slippery, or normal).
   The default values provided are calibrated for cold ice (-20°C) with a frozen bed.
   Different values may be needed for warmer ice or different bed conditions :cite:`Crawford2021`.

.. rubric:: Parameters

Prefix: ``calving.cliff_calving_tensile.``

.. pism-parameters::
   :prefix: calving.cliff_calving_tensile.



.. _sec-calving-thickness-threshold:

Calving of thin floating ice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The option :opt:`-calving thickness_calving` is based on the observation that ice shelf
calving fronts are commonly thicker than about 150--250 m (even though the physical
reasons are not clear yet). Accordingly, any floating ice thinner than `H_{\textrm{cr}}`
is removed along the front, at a rate at most one grid cell per time step. The value of
`H_{\mathrm{cr}}` can be set using the configuration parameter
:config:`calving.thickness_calving.threshold`.

To set a time-and-space dependent ice thickness threshold, set the parameter
:config:`calving.thickness_calving.file`. This file should contain the variable
:var:`thickness_calving_threshold` in meters.

.. rubric:: Parameters

Prefix: ``calving.thickness_calving.``

.. pism-parameters::
   :prefix: calving.thickness_calving.

.. _sec-calving-floating-ice:

Calving of all floating ice
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The option :opt:`-calving float_kill` removes (calves), at each time step of the run, any
ice that satisfies the flotation criterion. Use of this option implies that there are no
ice shelves in the model at all.

Set :config:`calving.float_kill.margin_only` to restrict this to cells at the ice margin.

Sometimes it is useful to preserve a one-cell-wide shelf near the grounding line. To do
this, set :config:`calving.float_kill.calve_near_grounding_line` to false.

.. rubric:: Parameters

Prefix: ``calving.float_kill.``

.. pism-parameters::
   :prefix: calving.float_kill.

.. _sec-prescribed-retreat:

Prescribed front retreat
^^^^^^^^^^^^^^^^^^^^^^^^

Option :opt:`-front_retreat_file` allows prescribing retreat of the ice front. The forcing
file specified using this option should contain :var:`land_ice_area_fraction_retreat` ---
a 2D field, possibly time-dependent, that contains ones in areas that may be covered by
ice and zeros in areas that have to be ice-free. Values between `0` and `1` allow for a
"partial" retreat on coarser grids.

More precisely, :var:`land_ice_area_fraction_retreat` is a mask prescribing the *maximum
ice extent* at a given time throughout a simulation; a certain rate of retreat can be
prescribed by creating a field with an appropriately decreasing maximum extent.

Changes in ice mass resulting from using this mechanism are reported as a part of the
*discharge* (:var:`tendency_of_ice_mass_due_to_discharge`).

.. note::

   This replaces the :literal:`ocean_kill` mechanism available in previous PISM versions.

.. rubric:: Parameters

Prefix: ``geometry.front_retreat.prescribed.``

.. pism-parameters::
   :prefix: geometry.front_retreat.prescribed.
