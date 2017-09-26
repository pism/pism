.. include:: ../../global.rst

.. _sec-marine:

Marine ice sheet modeling
=========================

PISM is often used to model whole ice sheets surrounded by ocean, with attached floating ice shelves, or smaller regions like outlet glaciers flowing into embayments and possibly generating floating tongues.  This section explains the geometry and stress balance mechanisms in PISM that apply to floating ice, at the vertical calving faces of floating ice, or at marine grounding lines.  The physics at calving fronts is very different from elsewhere on an ice sheet, because the flow is nothing like the lubrication flow addressed by the SIA, and nor is the physics like the sliding flow in the interior of an ice domain.  The needed physics at the calving front can be thought of as boundary condition modifications to the mass continuity equation and to the SSA stress balance equation.  The physics of grounding lines are substantially handled by recovering sub-grid information through interpolation.

.. _sec-pism-pik:

PIK options for marine ice sheets
---------------------------------

References :cite:`Albrechtetal2011`, :cite:`Levermannetal2012`, :cite:`Winkelmannetal2011` by the
research group of Prof. Anders Levermann at the Potsdam Institute for Climate Impact
Research ("PIK"), Germany, describe most of the mechanisms covered in this section. These
are all improvements to the grounded, SSA-as-a-sliding law model of :cite:`BBssasliding`. These
improvements make PISM an effective Antarctic model, as demonstrated by
:cite:`Golledgeetal2013`, :cite:`Martinetal2011`, :cite:`Winkelmannetal2012`, among other publications.
These improvements had a separate existence as the "PISM-PIK" model from 2009--2010, but
since PISM stable0.4 are part of PISM itself.

.. list-table:: Options which turn on PIK ice shelf front and grounding line mechanisms. A
                calving law choice is needed in addition to these options.
   :name: tab-pism-pik
   :header-rows: 1

   * - Option
     - Description

   * - :opt:`-cfbc`
     - apply the stress boundary condition along the ice shelf calving front
       :cite:`Winkelmannetal2011`

   * - :opt:`-kill_icebergs`
     - identify and eliminate free-floating icebergs, which cause well-posedness problems
       for the SSA stress balance solver :cite:`Winkelmannetal2011`

   * - :opt:`-part_grid`
     - allow the ice shelf front to advance by a part of a grid cell, avoiding
       the development of unphysically-thinned ice shelves :cite:`Albrechtetal2011` 

   * - :opt:`-subgl`
     - apply interpolation to compute basal shear stress and basal melt near the grounding line :cite:`Feldmannetal2014` 

   * - :opt:`-no_subgl_basal_melt`
     - **don't** apply interpolation to compute basal melt near the grounding line if
       :opt:`-subgl` is set :cite:`Feldmannetal2014`
    
   * - :opt:`-pik`
     - equivalent to option combination ``-cfbc -kill_icebergs -part_grid -subgl``

A summary of options to turn on most of these "PIK" mechanisms is in Table
:numref:`tab-pism-pik`. More information on the particular mechanisms is given in
sub-sections :ref:`sec-cfbc` through :ref:`sec-subgrid-grounding-line` that follow the
Table.

.. note::

   When in doubt, PISM users should set option :opt:`-pik` to turn on all of mechanisms in
   :numref:`tab-pism-pik`. The user should also choose a calving model from
   :numref:`tab-calving`. However, the :opt:`-pik` mechanisms will not be effective if the
   non-default FEM stress balance :opt:`-ssa_method fem` is chosen.

.. _sec-cfbc:

Stress condition at calving fronts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The vertically integrated force balance at floating calving fronts has been formulated by :cite:`Morland` as

.. math::
   :name: eq-cfbc

   \int_{z_s-\frac{\rho}{\rho_w}H}^{z_s+(1-\frac{\rho}{\rho_w})H}\mathbf{\sigma}\cdot\mathbf{n}\;dz
   = \int_{z_s-\frac{\rho}{\rho_w}H}^{z_s}\rho_w g (z-z_s) \;\mathbf{n}\;dz.

with `\mathbf{n}` being the horizontal normal vector pointing from the ice boundary oceanward, `\mathbf{\sigma}` the *Cauchy* stress tensor, `H` the ice thickness and `\rho` and `\rho_{w}` the densities of ice and seawater, respectively, for a sea level of `z_s`. The integration limits on the right hand side of equation :eq:`eq-cfbc` account for the pressure exerted by the ocean on that part of the shelf, which is below sea level (bending and torque neglected). The limits on the left hand side change for water-terminating outlet glacier or glacier fronts above sea level according to the bed topography.  By applying the ice flow law (section :ref:`sec-rheology`), equation :eq:`eq-cfbc` can be rewritten in terms of strain rates (velocity derivatives), as one does with the SSA stress balance itself.

Note that the discretized SSA stress balance, in the default finite difference discretization chosen by :opt:`-ssa_method` ``fd``, is solved with an iterative matrix scheme.  If option :opt:`-cfbc` is set then, during matrix assembly, those equations which are for fully-filled grid cells along the ice domain boundary have terms replaced according to equation :eq:`eq-cfbc`, so as to apply the correct stresses :cite:`Albrechtetal2011`, :cite:`Winkelmannetal2011`.

.. _sec-part-grid:

Partially-filled cells at the boundaries of ice shelves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Albrecht et al :cite:`Albrechtetal2011` argue that the correct movement of the ice shelf calving front on a finite-difference grid, assuming for the moment that ice velocities are correctly determined (see below), requires tracking some cells as being partially-filled (option :opt:`-part_grid`).  If the calving front is moving forward, for example, then the neighboring cell gets a little ice at the next time step.  It is not correct to add that little mass as a thin layer of ice which fills the cell's horizontal extent, as that would smooth the steep ice front after a few time steps.  Instead the cell must be regarded as having ice which is comparably thick to the upstream cells, but where the ice only partially fills the cell.

Specifically, the PIK mechanism turned on by ``-part_grid`` adds mass to the partially-filled cell which the advancing front enters, and it determines the coverage ratio according to the ice thickness of neighboring fully-filled ice shelf cells.  If option ``-part_grid`` is used then the PISM output file will have field ``Href`` which shows the amount of ice in the partially-filled cells as a thickness.  When a cell becomes fully-filled, in the sense that the ``Href`` thickness equals the average of neighbors, then the residual mass is redistributed to neighboring partially-filled or empty grid cells.

The stress balance equations determining the velocities are only sensitive to "fully-filled" cells.  Similarly, advection is controlled only by values of velocity in fully-filled cells.  Adaptive time stepping (specifically: the CFL criterion) limits the speed of ice front propagation so that at most one empty cell is filled, or one full cell emptied, per time step by the advance or retreat, respectively, of the calving front.

.. _sec-kill-icebergs:

Iceberg removal
^^^^^^^^^^^^^^^

Any calving mechanism (see section :ref:`sec-calving`) removes ice along the seaward front of the ice shelf domain.  This can lead to isolated cells either filled or partially-filled with floating ice, or to patches of floating ice (icebergs) fully surrounded by ice free ocean neighbors.  This ice is detached from the flowing and partly-grounded ice sheet.  That is, calving can lead to icebergs.

In terms of our basic model of ice as a viscous fluid, however, the stress balance for an iceberg is not well-posed because the ocean applies no resistance to balance the driving stress.  (See :cite:`SchoofStream`.)  In this situation the numerical SSA stress balance solver will fail.

Option :opt:`-kill_icebergs` turns on the mechanism which cleans this up.  This option is therefore generally needed if there is nontrivial calving.  The mechanism identifies free-floating icebergs by using a 2-scan connected-component labeling algorithm.  It then eliminates such icebergs, with the corresponding mass loss reported as a part of the 2D discharge flux diagnostic (see section :ref:`sec-saving-diagnostics`).

.. _sec-subgrid-grounding-line:

Sub-grid treatment of the grounding line position
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The command-line option :opt:`-subgl` turns on a parameterization of the grounding line
position based on the "LI" parameterization described in :cite:`Gladstoneetal2010` and
:cite:`Feldmannetal2014`. With this option PISM computes an extra flotation mask, available as
the ``gl_mask`` output variable, which corresponds to the fraction of the cell that is
grounded. Cells that are ice-free or fully floating are assigned the value of `0` while
fully-grounded icy cells get the value of `1`. Partially grounded cells, the ones which
contain the grounding line, get a value between `0` and `1`. The resulting field has two
uses:

- It is used to scale the basal friction in cells containing the grounding line in order
  to avoid an abrupt change in the basal friction from the "last" grounded cell to the
  "first" floating cell. See the source code browser for the detailed description and
  section :ref:`sec-MISMIP3d` for an application.
- It is used to adjust the basal melt rate in cells containing the grounding line: in such
  cells the basal melt rate is set to `M_{b,\text{adjusted}} = \lambda
  M_{b,\text{grounded}} + (1 - \lambda)M_{b,\text{shelf-base}}`, where `\lambda` is the
  value of the flotation mask. Use :opt:`-no_subgl_basal_melt` to disable this.


.. _sec-floatmask:

Flotation criterion, mask, and sea level
----------------------------------------

The most basic decision about marine ice sheet dynamics made internally by PISM is whether
a ice-filled grid cell is floating. That is, PISM applies the "flotation criterion"
:cite:`Winkelmannetal2011` at every time step and at every grid location to determine whether
the ice is floating on the ocean or not. The result is stored in the ``mask`` variable.
The ``mask`` variable has ``pism_intent`` = ``diagnostic``, and thus it does *not* need to
be included in the input file set using the ``-i`` option.

The possible values of the ``mask`` are given in :numref:`tab-maskvals`. The mask
does not *by itself* determine ice dynamics. For instance, even when ice is floating (mask
value ``MASK_FLOATING``), the user must turn on the usual choice for ice shelf dynamics,
namely the SSA stress balance, by using options :opt:`-stress_balance ssa` or
:opt:`-stress_balance ssa+sia`.

.. FIXME: this is certainly out of date
   
.. list-table:: The PISM mask, in combination with user options, determines the dynamical
                model.
   :name: tab-maskvals
   :header-rows: 1

   * - Mask value
     - Meaning

   * - 0 = ``MASK_ICE_FREE_BEDROCK``
     - ice free bedrock 

   * - 2 = ``MASK_GROUNDED``
     - ice is grounded 

   * - 3 = ``MASK_FLOATING``
     - ice is floating (the SIA is never applied; the SSA is applied if the ``ssa`` or
       ``ssa+sia`` stress balance model is selected

   * - 4 = ``MASK_ICE_FREE_OCEAN``
     - ice-free ocean 

Assuming that the geometry of the ice is allowed to evolve (which can be turned off by
option ``-no_mass``), and assuming an ocean exists so that a sea level is used in the
flotation criterion (which can be turned off by option :opt:`-dry`), then at each time
step the mask will be updated.

.. _sec-calving:

Calving
-------

.. _sec-eigen-calving:

Eigen calving
^^^^^^^^^^^^^

PISM-PIK introduced a physically-based 2D-calving parameterization :cite:`Levermannetal2012`. This calving parameterization is turned on in PISM by option :opt:`-calving eigen_calving`.  Average calving rates, `c`, are proportional to the product of principal components of the horizontal strain rates, `\dot{\epsilon}_{_\pm}`, derived from SSA-velocities 

.. math::
   :name: eq-calv2

   c = K\; \dot{\epsilon}_{_+}\; \dot{\epsilon}_{_-}\quad\text{and}\quad\dot{\epsilon}_{_\pm}>0\:.

The rate `c` is in `\text{m}\,\text{s}^{-1}`, and the principal strain rates
`\dot\epsilon_\pm` have units `\text{s}^{-1}`, so `K` has units `\text{m}\,\text{s}`. The
constant `K` incorporates material properties of the ice at the front. It can be set using
the :opt:`-eigen_calving_K` option or a configuration parameter (``eigen_calving_K`` in
``src/pism_config.cdl``).

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
   While eigen-calving (section :ref:`sec-eigen-calving`) is appropriate for Antartic
   ice shelves, it does not work for outlet glaciers that flow in narrow fjords. Along
   valleys with nearly parallel walls, the transverse component of the velocity is close
   to zero, and the transversal strain rate is therefore also close to zero and noisy.
   Instead of the product of the eigen strain rates, :cite:`Morlighem2016` proposes a calving
   law where the calving rate `c` is a functionally related to tensile stresses:

.. math::
   :name: eq-calv3

   c = |\mathbf{u}| \frac{\tilde{\sigma}}{\sigma_{max}},

where `\tilde{\sigma}` is the tensile von Mises stress and `\sigma_{max}` is a threshold
that has units `Pa`. It can be set as a configuration parameter
(:config:`calving.vonmises.sigma_max` in ``src/pism_config.cdl``). As the tensile fracture
strength is much smaller than the compressive fracture strenth, the effective tensile
strain rate is defined as

.. math::
   :name: eq-calv4

   \tilde{\dot{\epsilon}}_e = \left(\frac{1}{2}\left(\max(0,\dot{\epsilon}_{_+})^2 +
   \max(0,\dot{\epsilon}_{_-})^2\right)\right)^{1/2}.

Following :cite:`Morlighem2016` `\tilde{\sigma}` is given by

.. math::
   :name: eq-calv5

   \tilde{\sigma} = \sqrt{3} B \tilde{\dot{{\epsilon}}}_e^{1/n},

where `B` is the ice hardness.

Similar to ``eigen_calving``, the calving rate from ``vonmises_calving`` can be used to
limit the overall timestep of PISM--thus slowing down all of PISM--by using
:opt:`-calving_cfl`.

.. _sec-additional-calving:

Additional calving methods
^^^^^^^^^^^^^^^^^^^^^^^^^^

PISM also includes three more basic calving mechanisms (:numref:`tab-calving`). The
option :opt:`-calving thickness_calving` is based on the observation that ice shelf
calving fronts are commonly thicker than about 150--250\,m (even though the physical
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

Option :opt:`-calving ocean_kill` chooses the calving mechanism removing ice in the "open
ocean". It requires the option :opt:`-ocean_kill_file`, which specifies the file
containing the ice thickness field ``thk``. (This can be the input file specified using
``-i``.) Any locations which were ice-free (``thk == 0``) and which had bedrock elevation
below sea level (``topg < 0``), in the provided data set, are marked as ice-free ocean.
The resulting mask is not altered during the run, and is available as diagnostic field
``ocean_kill_mask``. At these places any floating ice is removed at each step of the run.
Ice shelves can exist in locations where a positive thickness was supplied in the provided
data set.

To select several calving mechanisms, use a comma-separated list of keywords mentioned in
:numref:`tab-calving`:

.. code-block:: none

   -calving eigen_calving,thickness_calving,ocean_kill,vonmises_calving

.. list-table:: Options for the four calving models in PISM.
   :name: tab-calving
   :header-rows: 1

   * - Option
     - Description
    
   * - :opt:`-calving eigen_calving`
     - Physically-based calving parameterization :cite:`Levermannetal2012`,
       :cite:`Winkelmannetal2011`. Whereever the product of principal strain rates is positive,
       the calving rate is proportional to this product.

   * - :opt:`-eigen_calving_K` (`m s`)
     - Sets the proportionality parameter `K` in `\text{m}\,\text{s}`.

   * - :opt:`-calving vonmises_calving`
     - Physically-based calving parameterization :cite:`Morlighem2016` that uses the tensile
       von Mises stresses.

   * - :opt:`-calving_cfl`
     - Apply CFL-type criterion to reduce (limit) PISM's time step, according to for
       stress-calving rate.

   * - :opt:`-vonmises_calving_sigma_max` (`Pa`)
     - Sets the maximum tensile stress `\tilde{\sigma}` in `\text{Pa}`.

   * - :opt:`-calving thickness_calving`
     - Calve all near-terminus ice which is thinner than ice threshold thickness
       `H_{\textrm{cr}}`.

   * - :opt:`-thickness_calving_threshold` (m)
     - Sets the thickness threshold `H_{\textrm{cr}}` in meters.

   * - :opt:`-calving float_kill`
     - All floating ice is calved off immediately. 

   * - :opt:`-calving ocean_kill`
     - All ice flowing into grid cells marked as "ice free ocean", according to the ice
       thickness in the provided file, is calved.

   * - :opt:`-ocean_kill_file`
     - Sets the file with the ``thk`` field used to compute maximum ice extent.

.. _sec-model-melange-pressure:

Modeling melange back-pressure
------------------------------


Equation :eq:`eq-cfbc` above, describing the stress boundary condition for ice shelves,
can be written in terms of velocity components:

.. math::
   :name: eq-cfbc-uv

   \newcommand{\psw}{p_{\text{ocean}}}
   \newcommand{\pice}{p_{\text{ice}}}
   \newcommand{\pmelange}{p_{\text{melange}}}
   \newcommand{\n}{\mathbf{n}}
   \newcommand{\nx}{\n_{x}}
   \newcommand{\ny}{\n_{y}}
   
   \begin{array}{lclcl}
     2 \nu H (2u_x + u_y) \nx &+& 2 \nu H (u_y + v_x)  \ny &=& \displaystyle \int_{b}^{h}(\pice - \psw) dz\, \nx,\\
     2 \nu H (u_y + v_x)  \nx &+& 2 \nu H (2v_y + u_x) \ny &=& \displaystyle \int_{b}^{h}(\pice - \psw) dz\, \ny.
   \end{array}

Here `\nu` is the vertically-averaged ice viscosity, `b` is the ice base elevation, `h` is
the ice top surface elevation, and `\psw` and `\pice` are pressures of the column of sea
water and ice, respectively.

We call the integral on the right hand side of :eq:`eq-cfbc-uv` the "pressure imbalance
term". To model the effect of melange :cite:`Amundsonetal2010` on the stress boundary
condition, we assume that the melange back-pressure `\pmelange` does not exceed `\pice -
\psw`. Therefore we introduce `\lambda \in [0,1]` (the melange back pressure fraction)
such that

.. math::

   \pmelange = \lambda (\pice - \psw).

Then melange pressure is added to the ordinary ocean pressure so that the pressure imbalance term scales with `\lambda`:

.. math::
   :name: eq-cfbc-3

   \int_{b}^{h}(\pice - (\psw + \pmelange))\, dz &= \int_{b}^{h}(\pice - (\psw + \lambda(\pice - \psw)))\, dz

   &= (1 - \lambda) \int_{b}^{h} (\pice - \psw)\, dz.

This formula replaces the right hand side of :eq:`eq-cfbc-uv`.

By default, `\lambda` is set to zero, but PISM implements a scalar time-dependent "melange
back pressure fraction offset" forcing in which `\lambda` can be read from a file. Please
see the *PISM's Climate Forcing Manual* for details.
