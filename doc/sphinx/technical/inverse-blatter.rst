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
`\tau_c` or vertically-averaged ice hardness `B` from observed surface
velocities, and the two can be alternated (:ref:`sec-inv-blatter-alternating`).
The key difference is that the forward model solves the full 3D
Blatter-Pattyn equations instead of the depth-integrated SSA.

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

The implementation lives in ``src/inverse/IP_BlatterForwardProblem.{hh,cc}``
(design-variable-agnostic part: forward solve, surface extraction, adjoint
solve) and its two derived classes ``IP_BlatterTaucForwardProblem`` and
``IP_BlatterHardavForwardProblem`` (design Jacobians). The user-facing driver
is ``pismi`` (``site-packages/PISM/pismi.py``).

.. _sec-inv-blatter-notation:

Notation
--------

.. list-table::

   * - `\zeta`
     - parameterized design variable (same as SSA, :ref:`sec-inv-ssa-notation`)
   * - `\tau_c`
     - basal yield stress: `\tau_c = g(\zeta)`
   * - `B`
     - vertically-averaged ice hardness (``hardav``, `B = A^{-1/n}`):
       `B = g(\zeta)`, column-constant
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
   * - `\eta`
     - effective viscosity `\eta = E\, \frac12 B\, (\epsilon + \gamma)^{(1-n)/(2n)}`
       (`\gamma` the second invariant of the strain rate, `E` the enhancement factor)

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

.. _sec-inv-blatter-hardav:

Ice hardness as the design variable
-----------------------------------

With ``-inv_design hardav`` the design variable is the vertically-averaged
ice hardness `B(x,y) = g(\zeta)`. As in the SSA case the hardness is a
single value per column; the Blatter solver replicates it over the sigma
grid (``Blatter::init_averaged_ice_hardness``, selected by
``stressbalance::Inputs::averaged_hardness``) instead of deriving the
hardness from enthalpy.

The hardness enters the residual only through the effective viscosity in the
**volume** term (``Blatter::residual_f``),

.. math::
   :label: eq-inv-blatter-Rf

   \mathcal{R}_{f}^{(t)} = \int_\Omega \eta(B, \gamma)\,
   \begin{pmatrix}
   \psi_{t,x}(4u_x + 2v_y) + \psi_{t,y}(u_y + v_x) + \psi_{t,z} u_z \\
   \psi_{t,x}(u_y + v_x) + \psi_{t,y}(2u_x + 4v_y) + \psi_{t,z} v_z
   \end{pmatrix} dV
   \equiv \int_\Omega \eta\, F(\uu, \psi_t)\, dV,

so, unlike `\tau_c`, it couples to *every* element of the ice column.
Because `\eta` is linear in `B`, the design Jacobian needs no new flow-law
derivative:

.. math::
   :label: eq-inv-blatter-Jdesign-B

   J_{\text{Design}}\, d\zeta = \int_\Omega \eta(dB, \gamma)\, F(\uu, \psi_t)\, dV,
   \qquad dB = g'(\zeta)\, d\zeta,

i.e. the viscous residual evaluated with `B` replaced by `dB` (this is the
same trick ``IP_SSAHardavForwardProblem`` uses). Its transpose maps the 3D
adjoint `\lambda = (\lambda_u, \lambda_v)` to the 2D design space by
**integrating over the whole ice column**:

.. math::
   :label: eq-inv-blatter-Jdesign-B-T

   (J_{\text{Design}}^T \lambda)_k
   = g'(\zeta_k) \sum_{\text{elements in column } k}\ \sum_q W_q\,
     \frac{\eta_q}{B_q}\,
     \Bigl[\lambda_{u,x}(4u_x + 2v_y) + \lambda_{u,y}(u_y + v_x) + \lambda_{u,z} u_z
     + \lambda_{v,x}(u_y + v_x) + \lambda_{v,y}(2u_x + 4v_y) + \lambda_{v,z} v_z\Bigr]_q
     \chi_k(q),

where `\lambda_{u,x} = \sum_t \lambda_{u,t}\, \psi_{t,x}` etc. are the
gradients of the interpolated adjoint field, `\eta_q / B_q = E\,\frac12
(\epsilon + \gamma_q)^{(1-n)/(2n)}`, and all eight nodes of every element
contribute to the 2D node of their column (``apply_jacobian_design_transpose_3d``
in ``IP_BlatterHardavForwardProblem.cc``). The surface extraction, adjoint
solve and reduced gradient are exactly as for `\tau_c` (next section).

The parameterization `g` uses ``inverse.design.param_hardav_scale`` (default
`10^8` Pa s\ :sup:`1/3`) and ``inverse.design.param_hardav_eps``; bound
constraints (``tikhonov_blmvm``) use ``inverse.stress_balance.hardav_min`` and
``hardav_max``. The ``exp`` parameterization is recommended: it keeps `B > 0`,
which the Blatter SNES requires. The ``zeta_fixed_mask`` computed by ``pismi``
frees the hardness on all icy cells (including floating ice, where it is the
only control on the velocity).

``examples/inverse/blatter_inverse_checks.py`` verifies the implementation with
an adjoint (dot-product) test `\langle DF\,\delta, w\rangle = \langle
\delta, DF^T w\rangle` and a finite-difference gradient test, for both design
variables and in parallel.

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

   Three methods are available, controlled by ``inverse.adjoint.method``
   (or ``-inv_adjoint_method``):

   - **approximate** (default): ``KSPSolve`` on the SNES Jacobian. The
     forward Newton Jacobian is symmetrized by the upper-triangle mirror
     in ``compute_jacobian``, so ``KSPSolve`` is a good approximation to
     ``KSPSolveTranspose``. Fast — reuses the existing matrix, no
     reassembly. Use GMRES (``-inv_adj_ksp_type gmres``).

   - **incomplete**: ``KSPSolve`` on a separately assembled Picard
     Jacobian (drops viscosity derivative terms via ``compute_picard_jacobian``).
     The matrix is truly symmetric, so CG works
     (``-inv_adj_ksp_type cg -inv_adj_pc_type gamg``). See
     :ref:`sec-inv-blatter-picard`.

   - **exact**: ``KSPSolveTranspose`` on the Newton Jacobian. Requires a
     transpose-compatible preconditioner (``-inv_adj_pc_type jacobi``);
     SOR and GAMG do not support transpose.

   All three methods use a standalone KSP (prefix ``inv_adj_``) rather than
   the SNES's multigrid KSP, avoiding MG hierarchy issues.

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

PISM provides three adjoint methods via the configuration parameter
``inverse.adjoint.method`` (command-line: ``-inv_adjoint_method``):

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 55

   * - Method
     - KSP call
     - Matrix
     - Notes
   * - ``approximate``
     - ``KSPSolve``
     - SNES Jacobian (Newton, symmetrized by upper-triangle mirror)
     - Default. Fastest — reuses existing matrix. Use GMRES + GAMG.
   * - ``incomplete``
     - ``KSPSolve``
     - Separate Picard Jacobian (assembled via ``compute_picard_jacobian``)
     - Truly symmetric. CG + GAMG works. Requires reassembly each iteration.
   * - ``exact``
     - ``KSPSolveTranspose``
     - SNES Jacobian (Newton)
     - Mathematically exact. Requires transpose-compatible PC (Jacobi).

The **approximate** method exploits the fact that the forward SNES assembles
only the upper triangle of the element Jacobian and mirrors it to the lower
triangle (see ``compute_jacobian``). The resulting matrix is symmetric
regardless of whether the Newton correction terms are present, so
``KSPSolve`` gives the same result as ``KSPSolveTranspose``.

The **incomplete** method assembles a separate Picard Jacobian by temporarily
setting ``m_use_picard = true`` and calling ``compute_picard_jacobian``. This
drops the `\eta_u F_u` terms in ``jacobian_f``, producing an element matrix
that is symmetric *before* the mirror step. The matrix is truly SPD, so CG
converges and GAMG works as a preconditioner.

The ISMIP-HOM twin experiment (``examples/inverse/ismiphom_twin.py``)
confirms that all three methods produce identical convergence histories and
recovered `\tau_c` fields, consistent with :cite:`Morlighem2013`.

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

.. _sec-inv-blatter-alternating:

Alternating tauc / hardav inversion
-----------------------------------

Surface velocities constrain basal drag and ice stiffness jointly, so
``pismi`` can alternate between the two design variables in one invocation:

.. code-block:: none

   pismi -i STATE.nc -inv_data OBS.nc -o OUT.nc -stress_balance.model blatter \
         -inverse.alternating_cycles 3 -inverse.alternating_misfit_tol 0.01 ...

Each cycle runs a `\tau_c` inversion with the current hardness held fixed,
followed by a hardness inversion with the new `\tau_c` held fixed. The
phases hand their results to each other through the output file:

- The first (`\tau_c`) phase uses the column-constant hardness ``hardav``
  from the input file if present, and otherwise computes it from enthalpy
  (``rheology::averaged_hardness_vec``) so that both phases see the same
  hardness model. That hardness is also saved as ``hardav_prior``.
- Every phase writes the physical fields ``tauc`` and ``hardav`` and its
  parameterized solution ``zeta_inv_tauc`` / ``zeta_inv_hardav``. The next
  phase reads the *other* variable from the output file and starts its own
  variable from the previous cycle's result; the Tikhonov priors
  (``tauc_prior``, ``hardav_prior``) stay anchored to the original inputs.
- Iteration histories are stored as ``inv_misfit_c<cycle>_<var>`` etc., the
  final misfit of each phase as the global attribute ``pismi_misfit_c<cycle>_<var>``,
  and the last completed phase as ``pismi_alternation_completed``.
  ``-inv_restart`` resumes from that phase.

The loop stops early when the relative misfit improvement over a full cycle
drops below ``inverse.alternating_misfit_tol``. Alternation requires
``stress_balance.model = blatter``: only the Blatter forward problems accept
both design variables as inputs (``IP_BlatterTaucForwardProblem`` uses a
``hardav`` field from ``Grid::variables()`` when available;
``IP_BlatterHardavForwardProblem`` uses ``tauc``).

Using the inverted hardness in forward runs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A forward ``pism`` run with the Blatter stress balance uses a prescribed
column-constant hardness instead of the enthalpy-derived one when
``stress_balance.averaged_hardness.enabled`` (``-use_averaged_hardness``) is
set. The field ``hardav`` is then a model state variable: it is read from the
input file on restart, can be regridded from an inversion output file, and is
written to output files (the ``hardav`` diagnostic reports it as well). A
typical forward leg following an alternating inversion is

.. code-block:: none

   pism -i STATE.nc -stress_balance.model blatter -use_averaged_hardness \
        -input.regrid.file OUT.nc -input.regrid.vars tauc,hardav \
        -basal_yield_stress.model constant ...

If ``hardav`` is neither in the input file nor listed in ``-input.regrid.vars``
the run stops with an error, and so does a run that sets the flag with a
stress balance other than Blatter (the SSA solvers always derive hardness
from enthalpy).

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

   The forward SNES can always use any MG smoother (including SOR). The adjoint
   solve uses a separate KSP (``inv_adj_`` prefix). Recommended adjoint settings:

   - ``approximate`` (default): ``-inv_adj_ksp_type gmres -inv_adj_pc_type gamg``
   - ``incomplete``: ``-inv_adj_ksp_type cg -inv_adj_pc_type gamg``
   - ``exact``: ``-inv_adj_ksp_type gmres -inv_adj_pc_type jacobi``

Element assembly
^^^^^^^^^^^^^^^^

The `\tau_c` design Jacobian assembly loops only over **basal elements**
(`k = 0` in the column loop), since `\tau_c` enters only through the basal
boundary integral; the hardness design Jacobian loops over **all** elements
of every column. Both are implemented in ``apply_jacobian_design_3d`` and
``apply_jacobian_design_transpose_3d`` of the respective forward problem.
These methods use local (ghosted) vectors for the 3D DMDA array assembly and
scatter back to global using ``DMLocalToGlobal`` with ``ADD_VALUES``; the
transposes accumulate into owned 2D nodes only. The state Jacobian is
re-assembled at the converged solution before it is used for linearizations
and adjoint solves.

The adjoint solve uses a standalone KSP (prefix ``inv_adj_``) operating on the
SNES Jacobian that was already assembled during the forward solve. This avoids
reusing the SNES's multigrid KSP (swapping operators on the MG KSP triggers
``PCSetUp_MG`` issues). Configure via ``-inv_adj_ksp_type``,
``-inv_adj_pc_type``, etc.

Key files
^^^^^^^^^

- ``src/inverse/IP_BlatterForwardProblem.{hh,cc}`` — shared base: forward
  solve, surface extraction, adjoint solve, reduced linearization
- ``src/inverse/IP_BlatterTaucForwardProblem.{hh,cc}`` — `\tau_c` design
  Jacobian (basal face)
- ``src/inverse/IP_BlatterHardavForwardProblem.{hh,cc}`` — hardness design
  Jacobian (volume / column integral)
- ``src/inverse/IP_Blatter{Tauc,Hardav}TaoTikhonovProblem.hh`` — Tikhonov
  specializations with bounds
- ``site-packages/PISM/invert/blatter.py`` — Python forward-run setup
- ``site-packages/PISM/invert/blatter_tao.py`` — Python TAO solver wrapper
- ``site-packages/PISM/pismi.py`` — unified inversion driver (SSA and
  Blatter; single design variable or alternation)
- ``examples/inverse/ismiphom_twin.py`` — ISMIP-HOM twin experiment comparing
  incomplete vs exact adjoint
- ``examples/inverse/blatter_inverse_checks.py`` — adjoint and gradient
  consistency checks

.. _sec-inv-blatter-limitations:

Limitations
-----------

- **Design variables**: `\tau_c` and column-constant hardness; a
  depth-varying hardness is not supported.
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
