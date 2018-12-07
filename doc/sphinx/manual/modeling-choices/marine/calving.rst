.. include:: ../../../global.txt

.. _sec-calving:

Calving
-------

.. contents::

The :numref:`tab-calving` summarizes options controlling calving parameterizations
implemented in PISM.

.. list-table:: Options for the calving models in PISM.
   :name: tab-calving
   :header-rows: 1
   :widths: 1,1

   * - Option
     - Description

   * - :opt:`-calving_cfl`
     - Apply CFL-type criterion to reduce (limit) PISM's time step using the horizontal
       calving rate computed by ``eigen_calving`` or ``vonmises_calving``.
    
   * - :opt:`-calving eigen_calving`
     - Physically-based calving parameterization :cite:`Levermannetal2012`,
       :cite:`Winkelmannetal2011`. Whereever the product of principal strain rates is positive,
       the calving rate is proportional to this product.

   * - :opt:`-eigen_calving_K` (`m s`)
     - Sets the proportionality parameter `K` in `\text{m}\,\text{s}`.

   * - :opt:`-calving vonmises_calving`
     - Physically-based calving parameterization :cite:`Morlighem2016` that uses the tensile
       von Mises stresses.

   * - :opt:`-vonmises_calving_sigma_max` (`Pa`)
     - Sets the maximum tensile stress `\tilde{\sigma}` in `\text{Pa}`.

   * - :opt:`-calving thickness_calving`
     - Calve all near-terminus ice which is thinner than ice threshold thickness
       `H_{\textrm{cr}}`.

   * - :opt:`-thickness_calving_threshold` (m)
     - Sets the thickness threshold `H_{\textrm{cr}}` in meters.

   * - :opt:`-thickness_calving_threshold_file`
     - Specifies the file containing the variable :opt:`calving_threshold` to be used as
       the spetially-variable thickness threshold

   * - :opt:`-calving float_kill`
     - All floating ice is calved off immediately.

   * - :opt:`-float_kill_margin_only`
     - At each time step, calve cells at the ice margin only instead of removing all
       floating ice.

   * - :opt:`-float_kill_calve_near_grounding_line`
     - Calve floating ice near the grounding line (this is the default). Disable using
       :opt:`-float_kill_calve_near_grounding_line off`.

   * - :opt:`-calving ocean_kill`
     - All ice flowing into grid cells marked as "ice free ocean", according to the ice
       thickness in the provided file, is calved.

   * - :opt:`-ocean_kill_file`
     - Sets the file with the ``thk`` field used to compute maximum ice extent.

To select several calving mechanisms, use a comma-separated list of keywords mentioned in
:numref:`tab-calving`:

.. code-block:: none

   -calving eigen_calving,thickness_calving,ocean_kill,vonmises_calving

.. _sec-eigen-calving:

Eigen calving
^^^^^^^^^^^^^

PISM-PIK introduced a physically-based 2D-calving parameterization
:cite:`Levermannetal2012`. This calving parameterization is turned on in PISM by option
:opt:`-calving eigen_calving`. Average calving rates, `c`, are proportional to the product
of principal components of the horizontal strain rates, `\dot{\epsilon}_{_\pm}`, derived
from SSA-velocities

.. math::
   :label: eq-calv2

   c = K\; \dot{\epsilon}_{_+}\; \dot{\epsilon}_{_-}\quad\text{and}\quad\dot{\epsilon}_{_\pm}>0\:.

The rate `c` is in `\text{m}\,\text{s}^{-1}`, and the principal strain rates
`\dot\epsilon_\pm` have units `\text{s}^{-1}`, so `K` has units `\text{m}\,\text{s}`. The
constant `K` incorporates material properties of the ice at the front. It can be set using
the :opt:`-eigen_calving_K` option or a configuration parameter
(:config:`calving.eigen_calving.K` in |config-cdl|).

The actual strain rate pattern strongly depends on the geometry and boundary conditions
along the confinements of an ice shelf (coast, ice rises, front position). The strain rate
pattern provides information in which regions preexisting fractures are likely to
propagate, forming rifts (in two directions). These rifts may ultimately intersect,
leading to the release of icebergs. This (and other) ice shelf calving models are not
intended to resolve individual rifts or calving events, but it produces
structurally-stable calving front positions which agree well with observations. Calving
rates balance calving-front ice flow velocities on average.

The partially-filled grid cell formulation (section :ref:`sec-part-grid`) provides a
framework suitable to relate the calving rate produced by ``eigen_calving`` to the mass
transport scheme at the ice shelf terminus. Ice shelf front advance and retreat due to
calving are limited to a maximum of one grid cell length per (adaptive) time step. The
calving rate (velocity) from ``eigen_calving`` can be used to limit the overall timestep
of PISM--thus slowing down all of PISM--by using :opt:`-calving_cfl`. This "CFL"-type
time-step limitation is definitely recommended in high-resolution runs which attempt to
model calving position accurately. Without this option, under certain conditions where
PISM's adaptive time step happens to be long enough, dendritic structures can appear at
the calving front because the calving mechanism cannot "keep up" with the computed calving
rate.

.. _sec-stress-calving:

Von Mises stress calving
^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   This code is experimental and has not yet been thoroughly tested, use at your own risk.
   
While eigen-calving (section :ref:`sec-eigen-calving`) is appropriate for Antartic ice
shelves, it does not work for outlet glaciers that flow in narrow fjords. Along valleys
with nearly parallel walls, the transverse component of the velocity is close to zero, and
the transversal strain rate is therefore also close to zero and noisy.

Instead of the product of the eigen strain rates, :cite:`Morlighem2016` proposes a calving
law where the calving rate `c` is a functionally related to tensile stresses:

.. math::
   :label: eq-calv3

   c = |\mathbf{u}| \frac{\tilde{\sigma}}{\sigma_{max}},

where `\tilde{\sigma}` is the tensile von Mises stress and `\sigma_{max}` is a threshold
that has units `Pa`. It can be set as a configuration parameter
(:config:`calving.vonmises.sigma_max` in |config-cdl|). As the tensile fracture strength
is much smaller than the compressive fracture strenth, the effective tensile strain rate
is defined as

.. math::
   :label: eq-calv4

   \tilde{\dot{\epsilon}}_e = \left(\frac{1}{2}\left(\max(0,\dot{\epsilon}_{_+})^2 +
   \max(0,\dot{\epsilon}_{_-})^2\right)\right)^{1/2}.

Following :cite:`Morlighem2016` `\tilde{\sigma}` is given by

.. math::
   :label: eq-calv5

   \tilde{\sigma} = \sqrt{3} B \tilde{\dot{{\epsilon}}}_e^{1/n},

where `B` is the ice hardness.

Similar to ``eigen_calving``, the calving rate from ``vonmises_calving`` can be used to
limit the overall timestep of PISM --- thus slowing down all of PISM --- by using
:opt:`-calving_cfl`.

.. _sec-additional-calving:

Additional calving methods
^^^^^^^^^^^^^^^^^^^^^^^^^^

PISM also includes three more basic calving mechanisms (:numref:`tab-calving`). The
option :opt:`-calving thickness_calving` is based on the observation that ice shelf
calving fronts are commonly thicker than about 150--250 m (even though the physical
reasons are not clear yet). Accordingly, any floating ice thinner than `H_{\textrm{cr}}`
is removed along the front, at a rate at most one grid cell per time step. The value of
`H_{\mathrm{cr}}` can be set using the :opt:`-thickness_calving_threshold` option or the
:config:`calving.thickness_calving.threshold` configuration parameter.

To set a spatially-variable ice thickness threshold, use the option
:opt:`-thickness_calving_threshold_file` or the parameter
:config:`calving.thickness_calving.threshold_file`. This file should contain the variable
:var:`calving_threshold` in meters (or other compatible units).

Option :opt:`-calving float_kill` removes (calves), at each time step of the run, any ice
that satisfies the flotation criterion. Use of this option implies that there are no ice
shelves in the model at all.

Use the option :opt:`-float_kill_margin_only` to restrict this to cells at the ice margin.

Sometimes it is useful to preserve a one-cell-wide shelf near the grounding line. To do
this, set :config:`calving.float_kill.calve_near_grounding_line` to false.

Option :opt:`-calving ocean_kill` chooses the calving mechanism removing ice in the "open
ocean". It requires the option :opt:`-ocean_kill_file`, which specifies the file
containing the ice thickness field ``thk``. (This can be the input file specified using
``-i``.) Any locations which were ice-free (``thk == 0``) and which had bedrock elevation
below sea level (``topg < 0``), in the provided data set, are marked as ice-free ocean.
The resulting mask is not altered during the run, and is available as diagnostic field
``ocean_kill_mask``. At these places any floating ice is removed at each step of the run.
Ice shelves can exist in locations where a positive thickness was supplied in the provided
data set.
