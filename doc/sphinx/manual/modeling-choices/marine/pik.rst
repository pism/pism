.. include:: ../../../global.txt

.. _sec-pism-pik:

PIK options for marine ice sheets
---------------------------------

.. contents::

References :cite:`Albrechtetal2011`, :cite:`Levermannetal2012`, :cite:`Winkelmannetal2011`
by the research group of Prof. Anders Levermann at the Potsdam Institute for Climate
Impact Research ("PIK"), Germany, describe most of the mechanisms covered in this section.
These are all improvements to the grounded, SSA-as-a-sliding law model of
:cite:`BBssasliding`. These improvements make PISM an effective Antarctic model, as
demonstrated by :cite:`Golledgeetal2013`, :cite:`Martinetal2011`,
:cite:`Winkelmannetal2012`, among other publications. These improvements had a separate
existence as the "PISM-PIK" model from 2009--2010, but since PISM stable0.4 are part of
PISM itself.

A summary of options to turn on most of these "PIK" mechanisms is in
:numref:`tab-pism-pik`. More information on the particular mechanisms is given in
sub-sections :ref:`sec-cfbc` through :ref:`sec-subgrid-grounding-line` that follow the
Table.

.. list-table:: Options which turn on PIK ice shelf front and grounding line mechanisms. A
                calving law choice is needed in addition to these options.
   :name: tab-pism-pik
   :header-rows: 1
   :widths: 1,3

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
     - apply interpolation to compute basal shear stress and basal melt near the grounding
       line :cite:`Feldmannetal2014`

   * - :opt:`-no_subgl_basal_melt`
     - **don't** apply interpolation to compute basal melt near the grounding line if
       :opt:`-subgl` is set :cite:`Feldmannetal2014`
    
   * - :opt:`-pik`
     - equivalent to option combination ``-cfbc -kill_icebergs -part_grid -subgl``

.. note::

   When in doubt, PISM users should set option :opt:`-pik` to turn on all of mechanisms in
   :numref:`tab-pism-pik`. The user should also choose a calving model from
   :numref:`tab-calving`. However, the :opt:`-pik` mechanisms will not be effective if the
   non-default FEM stress balance :opt:`-ssa_method fem` is chosen.

.. _sec-cfbc:

Stress condition at calving fronts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The vertically integrated force balance at floating calving fronts has been formulated by
:cite:`Morland` as

.. math::
   :label: eq-cfbc

   \int_{z_s-\frac{\rho}{\rho_w}H}^{z_s+(1-\frac{\rho}{\rho_w})H}\mathbf{\sigma}\cdot\mathbf{n}\;dz
   = \int_{z_s-\frac{\rho}{\rho_w}H}^{z_s}\rho_w g (z-z_s) \;\mathbf{n}\;dz.

with `\mathbf{n}` being the horizontal normal vector pointing from the ice boundary
oceanward, `\mathbf{\sigma}` the *Cauchy* stress tensor, `H` the ice thickness and `\rho`
and `\rho_{w}` the densities of ice and seawater, respectively, for a sea level of `z_s`.
The integration limits on the right hand side of equation :eq:`eq-cfbc` account for the
pressure exerted by the ocean on that part of the shelf, which is below sea level (bending
and torque neglected). The limits on the left hand side change for water-terminating
outlet glacier or glacier fronts above sea level according to the bed topography. By
applying the ice flow law (section :ref:`sec-rheology`), equation :eq:`eq-cfbc` can be
rewritten in terms of strain rates (velocity derivatives), as one does with the SSA stress
balance itself.

Note that the discretized SSA stress balance, in the default finite difference
discretization chosen by :opt:`-ssa_method` ``fd``, is solved with an iterative matrix
scheme. If option :opt:`-cfbc` is set then, during matrix assembly, those equations which
are for fully-filled grid cells along the ice domain boundary have terms replaced
according to equation :eq:`eq-cfbc`, so as to apply the correct stresses
:cite:`Albrechtetal2011`, :cite:`Winkelmannetal2011`.

.. _sec-part-grid:

Partially-filled cells at the boundaries of ice shelves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Albrecht et al :cite:`Albrechtetal2011` argue that the correct movement of the ice shelf
calving front on a finite-difference grid, assuming for the moment that ice velocities are
correctly determined (see below), requires tracking some cells as being partially-filled
(option :opt:`-part_grid`). If the calving front is moving forward, for example, then the
neighboring cell gets a little ice at the next time step. It is not correct to add that
little mass as a thin layer of ice which fills the cell's horizontal extent, as that would
smooth the steep ice front after a few time steps. Instead the cell must be regarded as
having ice which is comparably thick to the upstream cells, but where the ice only
partially fills the cell.

Specifically, the PIK mechanism turned on by :opt:`-part_grid` adds mass to the
partially-filled cell which the advancing front enters, and it determines the coverage
ratio according to the ice thickness of neighboring fully-filled ice shelf cells. If
option :opt:`-part_grid` is used then the PISM output file will have field
``ice_area_specific_volume`` which tracks the amount of ice in the partially-filled cells
as a "thickness", or, more appropriately, "volume per unit area". When a cell becomes
fully-filled, in the sense that the ``ice_area_specific_volume`` reaches the average of
the ice thickness in neighboring ice-filled cells, then the residual mass is redistributed
to neighboring partially-filled or empty grid cells.

The stress balance equations determining the velocities are only sensitive to
"fully-filled" cells. Similarly, advection is controlled only by values of velocity in
fully-filled cells. Adaptive time stepping (specifically: the CFL criterion) limits the
speed of ice front propagation so that at most one empty cell is filled, or one full cell
emptied, per time step by the advance or retreat, respectively, of the calving front.

.. _sec-kill-icebergs:

Iceberg removal
^^^^^^^^^^^^^^^

Any calving mechanism (see section :ref:`sec-calving`) removes ice along the seaward front
of the ice shelf domain. This can lead to isolated cells either filled or partially-filled
with floating ice, or to patches of floating ice (icebergs) fully surrounded by ice free
ocean neighbors. This ice is detached from the flowing and partly-grounded ice sheet. That
is, calving can lead to icebergs.

In terms of our basic model of ice as a viscous fluid, however, the stress balance for an
iceberg is not well-posed because the ocean applies no resistance to balance the driving
stress. (See :cite:`SchoofStream`.) In this situation the numerical SSA stress balance
solver will fail.

Option :opt:`-kill_icebergs` turns on the mechanism which cleans this up. This option is
therefore generally needed if there is nontrivial calving or significant variations in sea
level during a simulation. The mechanism identifies free-floating icebergs by using a
2-scan connected-component labeling algorithm. It then eliminates such icebergs, with the
corresponding mass loss reported as a part of the 2D discharge flux diagnostic (see
section :ref:`sec-saving-diagnostics`).

.. _sec-subgrid-grounding-line:

Sub-grid treatment of the grounding line position
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The command-line option :opt:`-subgl` turns on a parameterization of the grounding line
position based on the "LI" parameterization described in :cite:`Gladstoneetal2010` and
:cite:`Feldmannetal2014`. With this option PISM computes an extra flotation mask,
available as the :var:`cell_grounded_fraction` output variable, which corresponds to the
fraction of the cell that is grounded. Cells that are ice-free or fully floating are
assigned the value of `0` while fully-grounded icy cells get the value of `1`. Partially
grounded cells, the ones which contain the grounding line, get a value between `0` and
`1`. The resulting field has two uses:

- It is used to scale the basal friction in cells containing the grounding line in order
  to avoid an abrupt change in the basal friction from the "last" grounded cell to the
  "first" floating cell. See the source code browser for the detailed description and
  section :ref:`sec-MISMIP3d` for an application.
- It is used to adjust the basal melt rate in cells containing the grounding line: in such
  cells the basal melt rate is set to `M_{b,\text{adjusted}} = \lambda
  M_{b,\text{grounded}} + (1 - \lambda)M_{b,\text{shelf-base}}`, where `\lambda` is the
  value of the flotation mask. Use :opt:`-no_subgl_basal_melt` to disable this.

