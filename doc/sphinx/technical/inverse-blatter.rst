.. include:: ../global.txt

.. math::

.. _sec-inverse-blatter:

Blatter inversion: technical details
=====================================

.. contents::

.. _sec-inv-blatter-intro:

Introduction
------------

This section documents the inverse modeling framework for the Blatter
(higher-order) stress balance solver (:ref:`sec-blatter-details`). Like the
SSA inversion (:ref:`sec-inverse-ssa`), it estimates basal yield stress
`\tau_c` from observed surface velocities. The key difference is that the
forward model solves the full 3D Blatter-Pattyn equations instead of the
depth-integrated SSA.

Because the Blatter solver produces a 3D velocity field while observations are
2D surface velocities, the forward map includes an explicit **surface
extraction operator** `P` that projects from 3D to 2D. This introduces the
``3D-to-2D`` adjoint structure: the adjoint solve operates on the full 3D
state, but the design variable (tauc) is 2D.

An important simplification is the **incomplete (Picard) adjoint**: dropping the
viscosity-derivative terms from the Jacobian produces a symmetric matrix. This
means ``KSPSolve`` and ``KSPSolveTranspose`` are equivalent, any preconditioner
works (including SOR), and the adjoint gradient is nearly identical to the exact
Newton adjoint. See :ref:`sec-inv-blatter-picard`.

The implementation lives in ``src/inverse/IP_BlatterTaucForwardProblem.{hh,cc}``.
The user-facing driver is ``examples/inverse/pismi.py``.

.. _sec-inv-blatter-notation:

Notation
--------

.. list-table::

   * - `\zeta`
     - parameterized design variable (same as SSA, :ref:`sec-inv-ssa-notation`)
   * - `\tau_c`
     - basal yield stress: `\tau_c = g(\zeta)`
   * - `\uu_{3D}`
     - 3D Blatter velocity: `\uu_{3D} = (u(x,y,z),\, v(x,y,z))`
   * - `\uu_s`
     - 2D surface velocity: `\uu_s = P\, \uu_{3D}`
   * - `P`
     - surface extraction operator (evaluates at `z = s(x,y)`)
   * - `\mathcal{R}`
     - 3D Blatter residual
   * - `J_{\text{State}}`
     - 3D state Jacobian `\partial \mathcal{R} / \partial \uu_{3D}`
   * - `J_{\text{Design}}`
     - design Jacobian `\partial \mathcal{R} / \partial \zeta`
   * - `\Gamma_b`
     - basal boundary (bottom face, `k = 0` in the sigma grid)
   * - `M_z`
     - number of sigma levels (``stress_balance.blatter.Mz``)
   * - `\beta`
     - basal resistance coefficient: `\beta = \beta(\tau_c, |\uu|)`
   * - `\psi`
     - 3D finite-element test function (Q1 hexahedron basis)

.. _sec-inv-blatter-forward:

The forward map
---------------

The Blatter forward map composes the 3D stress-balance solve with surface
extraction:

.. math::
   :label: eq-inv-blatter-F

   F(\zeta) = P\, \uu_{3D}(\zeta),

where `\uu_{3D}` solves the Blatter-Pattyn residual equation
`\mathcal{R}(\uu_{3D}, \zeta) = 0` (see :ref:`sec-bp-intro` for the governing
equations), and `P` extracts the surface velocity from the sigma grid. In the
discretization, `P` simply reads the velocity at the top sigma level
(`k = M_z - 1`):

.. math::

   \uu_s(x_i, y_j) = \uu_{3D}(x_i, y_j, z_{M_z - 1}).

The design parameterization `\tau_c = g(\zeta)` is identical to the SSA
case (see :ref:`sec-inv-ssa-forward`).

.. _sec-inv-blatter-basal:

Where tauc enters the Blatter residual
---------------------------------------

The basal yield stress enters only through the **basal boundary condition**,
implemented in ``Blatter::residual_basal``. For grounded ice, the basal face
integral adds to the residual:

.. math::
   :label: eq-inv-blatter-Rb

   \mathcal{R}_{\text{basal}}^{(t)}
   = \int_{\Gamma_b}
     \beta(\tau_c, |\uu|)\, \uu \cdot \psi_t \, dS,

where `\psi_t` is the 3D finite-element test function evaluated on the basal
face and `\beta` is the basal resistance coefficient from the sliding law. For
the pseudo-plastic law,

.. math::

   \beta(\tau_c, |\uu|) = \tau_c\, f(|\uu|),

where `f(|\uu|)` is a nonlinear function of the velocity magnitude. This
integral is nonzero only at the bottom of the ice column (`k = 0`).

.. _sec-inv-blatter-jacobians:

State and design Jacobians
--------------------------

The **3D state Jacobian** `J_{\text{State}}` is the Newton Jacobian of the
full Blatter system, assembled via ``Blatter::compute_jacobian`` and used by
PETSc's SNES solver during the forward solve. It couples all vertical levels.

The **design Jacobian** `J_{\text{Design}}` maps perturbations of `\zeta`
(2D) into perturbations of the 3D residual. Since `\tau_c` enters only at the
basal boundary :eq:`eq-inv-blatter-Rb`, the design Jacobian is sparse in the
vertical: only bottom-face elements contribute. The nonzero entries are

.. math::
   :label: eq-inv-blatter-Jdesign

   (J_{\text{Design}}\, d\zeta)^{(t)}
   = \int_{\Gamma_b}
     \frac{\partial\beta}{\partial\tau_c}\,
     \uu \cdot \psi_t\,
     g'(\zeta)\, d\zeta \, dS,

where `\partial\beta/\partial\tau_c = f(|\uu|)` (the sliding law evaluated
at unit tauc), and the test functions `\psi_t` on the bottom face are nonzero
for the bottom 4 nodes of each hexahedral element.

The transpose `J_{\text{Design}}^T` maps a 3D adjoint variable `\lambda` to a
2D design perturbation:

.. math::
   :label: eq-inv-blatter-Jdesign-T

   (J_{\text{Design}}^T\, \lambda)_k
   = g'(\zeta_k) \sum_q W_q\,
     \frac{\partial\beta}{\partial\tau_c}\bigg|_q\,
     (\lambda_q \cdot \uu_q)\, \psi_k(q),

where the sum is over basal face quadrature points `q`, and `\lambda_q`,
`\uu_q` are the adjoint and velocity fields evaluated at those points.

.. _sec-inv-blatter-reduced:

Reduced gradient and the 3D-to-2D adjoint
------------------------------------------

The reduced gradient for the composite forward map
:eq:`eq-inv-blatter-F` is

.. math::
   :label: eq-inv-blatter-DF

   DF = P \cdot \bigl(-J_{\text{State}}^{-1}\, J_{\text{Design}}\bigr).

Its transpose, needed for the Tikhonov gradient, is

.. math::
   :label: eq-inv-blatter-DFt

   DF^T = -\,J_{\text{Design}}^T\, J_{\text{State}}^{-T}\, P^T.

To apply `DF^T` to a state-space perturbation `d\uu_s` (2D surface velocity),
the implementation (``IP_BlatterTaucForwardProblem::apply_linearization_transpose``)
proceeds in three steps:

1. **Inject** `d\uu_s` into 3D: compute `\mathbf{r}_{3D} = P^T d\uu_s`,
   which is zero everywhere except at the surface level (`k = M_z - 1`).

2. **Adjoint solve**: solve the 3D linear system

   .. math::

      J_{\text{State}}^T\, \lambda = P^T\, d\uu_s

   Two modes are available, controlled by ``inverse.use_incomplete_adjoint``:

   - **Exact adjoint** (``no``): uses ``KSPSolveTranspose`` on the Newton
     Jacobian. Requires a transpose-compatible preconditioner (SOR does not
     support transpose; use ``-inv_adj_pc_type jacobi``).

   - **Incomplete adjoint** (``yes``, default): uses ``KSPSolve`` on the
     (symmetric) Picard Jacobian. Any preconditioner works. See
     :ref:`sec-inv-blatter-picard`.

   Both modes use a standalone KSP (prefix ``inv_adj_``) rather than the
   SNES's multigrid KSP, avoiding MG hierarchy issues.

3. **Design Jacobian transpose**: compute

   .. math::

      d\zeta = -\,J_{\text{Design}}^T\, \lambda

   by iterating over basal face elements (``apply_jacobian_design_transpose_3d``).

The forward linearization `DF\, d\zeta` (``apply_linearization``) follows the
analogous three steps: apply `J_{\text{Design}}\, d\zeta` (nonzero only at
basal nodes), solve `J_{\text{State}}\, d\uu_{3D} = -J_{\text{Design}} d\zeta`,
then extract the surface: `d\uu_s = P\, d\uu_{3D}`.

.. _sec-inv-blatter-picard:

Incomplete (Picard) adjoint
---------------------------

The Blatter Jacobian in ``jacobian.cc`` consists of two contributions:

.. math::

   J_{\text{State}}^{(ts)} = \underbrace{\eta\, F_{uu}}_{\text{Picard}}
   + \underbrace{\eta_u\, F_u}_{\text{Newton correction}},

where `\eta` is the effective viscosity, `\eta_u = d\eta/d\gamma \cdot
d\gamma/du` involves the viscosity derivative, and `F_{uu}`, `F_u` are
strain-rate terms (see ``jacobian_f`` in ``jacobian.cc``). The Picard terms
are **symmetric** in the test/trial function indices (`s \leftrightarrow t`),
while the Newton correction terms are not.

The **incomplete adjoint** approximation drops the Newton correction terms,
yielding a symmetric Jacobian:

.. math::

   J_{\text{Picard}}^{(ts)} = \eta\, F_{uu}.

This approximation was introduced by :cite:`Morlighem2013` for the Antarctic
ice sheet and shown to produce adjoint gradients within 4% of the exact
Newton adjoint. The key advantages:

1. **Symmetry**: `J_{\text{Picard}} = J_{\text{Picard}}^T`, so ``KSPSolve``
   and ``KSPSolveTranspose`` are equivalent.
2. **Any preconditioner**: no transpose-compatibility requirement, so the
   standard MG+SOR smoother works.
3. **Lower cost**: no need to compute `d\eta/d\gamma`.

In PISM, the forward SNES always assembles the upper triangle and mirrors it
to the lower triangle (line "fill the lower-triangular part" in
``compute_jacobian``). This means the assembled Jacobian is effectively the
symmetrized Newton Jacobian — close to the Picard Jacobian. For the adjoint
solve, using ``KSPSolve`` on this symmetric matrix is equivalent to the
incomplete adjoint.

The configuration flag ``inverse.use_incomplete_adjoint`` (default ``yes``)
selects between:

- ``yes``: ``KSPSolve`` on the (symmetric) SNES Jacobian
- ``no``: ``KSPSolveTranspose`` on the same matrix (transpose-compatible
  preconditioner required via ``-inv_adj_pc_type jacobi``)

The ISMIP-HOM twin experiment (``examples/inverse/ismiphom_twin.py``)
confirms that both methods produce identical convergence histories and
recovered `\tau_c` fields.

.. _sec-inv-blatter-tikhonov:

The Tikhonov objective
----------------------

The Tikhonov formulation is identical to the SSA case
(:ref:`sec-inv-ssa-tikhonov`):

.. math::
   :label: eq-inv-blatter-J

   \mathcal{J}(\zeta) =
   \mathcal{J}_{\text{state}}\bigl(F(\zeta) - \uu_{\text{obs}}\bigr)
   + \frac{1}{\eta}\,\mathcal{J}_{\text{design}}(\zeta - \zeta_0),

with the only difference being that `F(\zeta)` now involves a 3D Blatter
solve instead of a 2D SSA solve. The functionals `\mathcal{J}_{\text{state}}` and
`\mathcal{J}_{\text{design}}` are the same 2D functionals described in
:ref:`sec-inv-ssa-functionals`.

.. _sec-inv-blatter-impl:

Implementation notes
--------------------

Solver configuration
^^^^^^^^^^^^^^^^^^^^

The Blatter forward solve is configured via PETSc command-line options with the
``bp_`` prefix. Recommended settings for the inversion:

.. list-table::
   :header-rows: 1

   * - Option
     - Value
     - Purpose
   * - ``-bp_pc_type mg``
     - multigrid
     - preconditioner for the 3D Blatter Jacobian
   * - ``-bp_pc_mg_levels 3``
     -
     - number of multigrid levels
   * - ``-bp_mg_coarse_ksp_type preonly``
     -
     - direct solve on the coarsest level
   * - ``-bp_mg_coarse_pc_type lu``
     -
     - LU factorization on coarsest level
   * - ``-bp_mg_levels_ksp_type chebyshev``
     -
     - smoother type
   * - ``-bp_snes_rtol 0.001``
     -
     - SNES relative tolerance
   * - ``-bp_ksp_rtol 0.001``
     -
     - KSP relative tolerance
   * - ``-inv_adj_ksp_type cg``
     -
     - adjoint KSP type (CG works well for symmetric Picard)
   * - ``-inv_adj_pc_type gamg``
     -
     - adjoint preconditioner (algebraic multigrid)

.. note::

   With the default incomplete adjoint (``inverse.use_incomplete_adjoint yes``),
   the forward SNES can use any MG smoother (including SOR). The adjoint
   solve uses a separate KSP with its own preconditioner (``inv_adj_`` prefix).
   Only when using the exact adjoint (``inverse.use_incomplete_adjoint no``)
   does the adjoint KSP require a transpose-compatible preconditioner
   (e.g., ``-inv_adj_pc_type jacobi``).

Element assembly
^^^^^^^^^^^^^^^^

The design Jacobian assembly loops only over **basal elements** (`k = 0` in
the column loop), since `\tau_c` enters only through the basal boundary
integral. This is implemented in ``apply_jacobian_design_3d`` and
``apply_jacobian_design_transpose_3d``. These methods use local (ghosted)
vectors for the 3D DMDA array assembly and scatter back to global using
``DMLocalToGlobal`` with ``ADD_VALUES``.

The adjoint solve uses a standalone KSP (prefix ``inv_adj_``) operating on the
SNES Jacobian that was already assembled during the forward solve. This avoids
reusing the SNES's multigrid KSP (swapping operators on the MG KSP triggers
``PCSetUp_MG`` issues). Configure via ``-inv_adj_ksp_type``,
``-inv_adj_pc_type``, etc.

Key files
^^^^^^^^^

- ``include/pism/inverse/IP_BlatterTaucForwardProblem.hh`` — class definition
- ``src/inverse/IP_BlatterTaucForwardProblem.cc`` — forward, adjoint, and
  linearization implementation
- ``include/pism/inverse/IP_BlatterTaucTaoTikhonovProblem.hh`` — Tikhonov
  specialization with tauc bounds
- ``site-packages/PISM/invert/blatter.py`` — Python forward-run setup
- ``site-packages/PISM/invert/blatter_tao.py`` — Python TAO solver wrapper
- ``examples/inverse/pismi.py`` — unified inversion driver (SSA and Blatter)
- ``examples/inverse/ismiphom_twin.py`` — ISMIP-HOM twin experiment comparing
  incomplete vs exact adjoint

.. _sec-inv-blatter-limitations:

Limitations
-----------

- **Design variable**: only `\tau_c` is currently supported (not hardness).
- **State space**: observations are matched against 2D surface velocity only,
  not depth-resolved velocity profiles.
- **H1 regularization and periodic BCs**: the ``IPGroundedIceH1NormFunctional2S``
  does not wrap around periodic boundaries, causing edge artifacts in inversions
  on periodic domains (e.g., ISMIP-HOM). Use L2-only regularization
  (``-inverse.design.cH1 0 -inverse.design.cL2 1``) for periodic problems, at
  the cost of slower convergence.
- **Computational cost**: each TAO iteration requires a full 3D Blatter forward
  solve (`\sim 5\text{--}15` SNES iterations) plus an adjoint solve.
  Inversions are significantly more expensive than the SSA equivalent.

.. _sec-inv-blatter-references:

References
----------

The Blatter inversion extends the SSA inverse framework of
:cite:`Maxwelletal2008` and :cite:`Habermannetal2013` to the higher-order
Blatter-Pattyn equations (:cite:`BrownSmithAhmadia2013`, :cite:`Tezaur2015`).
The 3D adjoint approach for ice-sheet inversions is also discussed in
:cite:`Goldberg2011`. The incomplete (Picard) adjoint approximation is
described and validated in :cite:`Morlighem2013`.
