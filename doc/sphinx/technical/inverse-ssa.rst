.. include:: ../global.txt

.. math::

.. _sec-inverse-ssa:

SSA inversion: technical details
================================

.. contents::

.. _sec-inv-ssa-intro:

Introduction
------------

PISM provides an inverse modeling framework for estimating the basal yield
stress `\tau_c` (or, alternatively, the vertically-averaged ice hardness)
from observed surface velocities, using the Shallow Shelf Approximation (SSA)
as the forward model. The implementation is described in
:cite:`Maxwelletal2008` and :cite:`Habermannetal2013`.

This section describes the mathematical formulation. The implementation lives
in ``src/inverse/`` and is exposed to Python via ``PISM.invert.ssa``. The user-
facing driver is ``examples/inverse/pismi_ssa.py``.

For the analogous Blatter inversion, see :ref:`sec-inverse-blatter`.

.. _sec-inv-ssa-notation:

Notation
--------

.. list-table::

   * - `\zeta`
     - parameterized design variable (the unknown)
   * - `\tau_c`
     - basal yield stress (the physical design variable)
   * - `g`
     - parameterization function: `\tau_c = g(\zeta)`
   * - `\pi`
     - projection that fixes `\zeta` at locations where `\tau_c` is known
   * - `\uu`
     - SSA velocity (the state variable)
   * - `\uu_{\text{obs}}`
     - observed surface velocity
   * - `F`
     - forward map: `F(\zeta) = \uu_{\text{SSA}}(g(\pi(\zeta)))`
   * - `\mathcal{R}`
     - SSA residual: `\mathcal{R}(\uu, \zeta) = 0` defines the forward solution
   * - `J_{\text{State}}`
     - state Jacobian `\partial \mathcal{R} / \partial \uu`
   * - `J_{\text{Design}}`
     - design Jacobian `\partial \mathcal{R} / \partial \zeta`
   * - `DF`
     - reduced gradient `\partial F / \partial \zeta`
   * - `\mathcal{J}_{\text{state}}`
     - state misfit functional (e.g., mean-square)
   * - `\mathcal{J}_{\text{design}}`
     - design regularization functional (e.g., `H^1`)
   * - `\eta`
     - Tikhonov penalty weight
   * - `\beta`
     - basal resistance coefficient: `\beta = \beta(\tau_c, |\uu|)`

.. _sec-inv-ssa-forward:

The forward map and design parameterization
-------------------------------------------

The basal yield stress must satisfy `\tau_c \ge 0`. Rather than enforcing this
as a constraint, we parameterize `\tau_c` by an unconstrained variable `\zeta`
through a non-negative function `g`:

.. math::
   :label: eq-inv-ssa-zeta

   \tau_c = g(\zeta).

PISM provides four parameterizations (in ``IPDesignVariableParameterization.hh``):

.. list-table::
   :header-rows: 1

   * - Name
     - `g(\zeta)`
     - `g'(\zeta)`
   * - ``ident``
     - `s\,\zeta`
     - `s`
   * - ``square``
     - `s\,\zeta^2`
     - `2 s\,\zeta`
   * - ``exp``
     - `s\,\exp(\zeta)`
     - `s\,\exp(\zeta)`
   * - ``trunc``
     - smooth `s\,\zeta` (linear for large `|\zeta|`, square near zero)
     - (smooth)

where `s` is a scale factor read from the configuration parameter
``inverse.design.param_tauc_scale``. The ``exp`` parameterization is the
default and guarantees `\tau_c > 0` strictly.

A second projection `\pi` clamps `\zeta` to known values at locations where
`\tau_c` is fixed (e.g., zero on floating ice; a high value on bare bedrock).
This is controlled by an integer mask ``zeta_fixed_mask`` and the configuration
flag ``inverse.use_zeta_fixed_mask``.

The forward map combining these is

.. math::
   :label: eq-inv-ssa-F

   F(\zeta) = \uu_{\text{SSA}}\bigl(g(\pi(\zeta))\bigr).

.. _sec-inv-ssa-jacobians:

State and design Jacobians
--------------------------

The SSA solver computes `\uu_{\text{SSA}}` by solving the residual equation

.. math::
   :label: eq-inv-ssa-residual

   \mathcal{R}(\uu;\, \tau_c,\, \text{other parameters}) = 0,

typically with Newton's method (the SNES inside ``SSAFEM``). Define the
implicit residual function

.. math::

   \mathcal{R}(\uu, \zeta) = \mathcal{R}(\uu;\, g(\pi(\zeta)),\, \dots).

Two Jacobians appear in the inverse problem:

The **state Jacobian** is `\partial\mathcal{R}/\partial\uu`,

.. math::
   :label: eq-inv-ssa-Jstate

   (J_{\text{State}})_{ij} = \frac{\partial \mathcal{R}_i}{\partial U_j},

where `U_j` are the velocity degrees of freedom. This is exactly the Newton
Jacobian assembled by SSAFEM during a forward solve, available through
``IP_SSATaucForwardProblem::assemble_jacobian_state``.

The **design Jacobian** is `\partial\mathcal{R}/\partial\zeta`,

.. math::
   :label: eq-inv-ssa-Jdesign

   (J_{\text{Design}})_{ik} = \frac{\partial \mathcal{R}_i}{\partial Z_k}.

Since `\tau_c` enters the SSA residual only through the basal drag term
`\beta(\tau_c, |\uu|)\,\uu`, the chain rule gives

.. math::
   :label: eq-inv-ssa-Jdesign-chain

   \frac{\partial \mathcal{R}}{\partial \zeta_k}
   = \frac{\partial \mathcal{R}}{\partial \tau_c}\, g'(\zeta_k) \,\phi_k
   = \frac{\partial \beta}{\partial \tau_c}\, \uu \, g'(\zeta_k)\, \phi_k

at each finite-element node `k` (with basis function `\phi_k`). For the common
pseudo-plastic sliding law `\beta(\tau_c, |\uu|) = \tau_c\, f(|\uu|)`, this
simplifies to `\partial\beta/\partial\tau_c = f(|\uu|) = \beta/\tau_c`. The
implementation uses the equivalent expression
``m_basal_sliding_law->drag(1.0, u, v)`` which evaluates `f(|\uu|)` directly.

.. _sec-inv-ssa-reduced:

Reduced gradient
----------------

The forward map satisfies the implicit equation `\mathcal{R}(F(\zeta), \zeta) = 0`.
Differentiating with respect to `\zeta` gives

.. math::

   J_{\text{State}}\, DF\, d\zeta + J_{\text{Design}}\, d\zeta = 0,

so the **reduced gradient** is

.. math::
   :label: eq-inv-ssa-DF

   DF = -\,J_{\text{State}}^{-1}\, J_{\text{Design}}.

For the inverse problem we need both `DF\, d\zeta` (forward linearization,
``apply_linearization``) and `DF^T\, d\uu` (adjoint linearization,
``apply_linearization_transpose``):

.. math::
   :label: eq-inv-ssa-DFt

   DF^T = -\,J_{\text{Design}}^T\, J_{\text{State}}^{-T}.

The transpose `J_{\text{State}}^{-T}` is realized by ``KSPSolveTranspose``
(PETSc) on the assembled state Jacobian.

.. _sec-inv-ssa-tikhonov:

The Tikhonov objective
----------------------

The inverse problem is formulated as Tikhonov-regularized least squares:

.. math::
   :label: eq-inv-ssa-tikhonov

   \min_{\zeta}\quad
   \mathcal{J}(\zeta)
   \;=\;
   \mathcal{J}_{\text{state}}\bigl(F(\zeta) - \uu_{\text{obs}}\bigr)
   \;+\;
   \frac{1}{\eta}\,
   \mathcal{J}_{\text{design}}(\zeta - \zeta_0),

where `\zeta_0` is the prior estimate (typically derived from the model's
initial `\tau_c`) and `\eta > 0` is the penalty weight. Larger `\eta` puts
more weight on data fit; smaller `\eta` favors the prior.

The gradient of `\mathcal{J}` with respect to `\zeta` is

.. math::
   :label: eq-inv-ssa-grad

   \nabla \mathcal{J}
   = DF^T \,\nabla \mathcal{J}_{\text{state}}\bigl(F(\zeta) - \uu_{\text{obs}}\bigr)
   + \frac{1}{\eta}\,\nabla \mathcal{J}_{\text{design}}(\zeta - \zeta_0).

Computing `DF^T \cdot v` requires one adjoint solve per gradient evaluation
(the dominant cost of each TAO iteration).

.. _sec-inv-ssa-functionals:

Functionals
-----------

State (misfit) functionals
^^^^^^^^^^^^^^^^^^^^^^^^^^

Three options for `\mathcal{J}_{\text{state}}` are available, controlled by
``inverse.state_func``:

- **``meansquare``** (default): weighted mean-square misfit normalized by
  ``inverse.stress_balance.velocity_scale``,

  .. math::

     \mathcal{J}_{\text{state}}(\uu) = \frac{1}{V^2 |\Omega|} \int_\Omega w(x)\, |\uu|^2\, dA,

  where `V` is the velocity scale and `w(x)` is an optional misfit weight
  (e.g., to mask out unobserved areas).

- **``log_ratio``**: relative-error functional based on logarithm of speed
  ratios; useful when the velocity range spans many orders of magnitude.

- **``log_relative``**: relative-error functional that smoothly transitions
  between absolute and relative error.

Design (regularization) functionals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two options for `\mathcal{J}_{\text{design}}`, controlled by ``inverse.design.func``:

- **``sobolevH1``** (default): a weighted `H^1` Sobolev norm,

  .. math::
     :label: eq-inv-ssa-H1

     \mathcal{J}_{\text{design}}(\zeta) = \frac{1}{|\Omega|} \int_\Omega
     \Bigl( c_{L^2}\, \zeta^2 + c_{H^1}\, L^2\, |\nabla \zeta|^2 \Bigr) dA,

  where `c_{L^2}` and `c_{H^1}` are weights (configuration parameters
  ``inverse.design.cL2``, ``inverse.design.cH1``) and `L` is a length scale
  (``inverse.stress_balance.length_scale``). The `L^2` factor on the gradient
  term makes the two terms dimensionally consistent and allows interpretation
  as a smoothing scale.

- **``tv``**: total-variation regularization, useful when the recovered field
  is expected to have sharp transitions.

For tauc inversion, an alternative ``IPGroundedIceH1NormFunctional`` restricts
integration to grounded-ice points only (used when
``-inv_ssa_grounded_ice_tauc`` is set).

.. _sec-inv-ssa-algorithms:

Algorithms
----------

Several minimization algorithms are available for solving
:eq:`eq-inv-ssa-tikhonov`, selected via ``inverse.stress_balance.method``:

.. list-table::
   :header-rows: 1

   * - Method
     - Class
     - Description
   * - ``tikhonov_lmvm``
     - PETSc TAO
     - L-BFGS quasi-Newton (recommended default)
   * - ``tikhonov_cg``
     - PETSc TAO
     - Conjugate gradient
   * - ``tikhonov_blmvm``
     - PETSc TAO
     - Bound-constrained L-BFGS (uses
       ``inverse.stress_balance.tauc_min/max``)
   * - ``tikhonov_lcl``
     - PETSc TAO
     - Linearly constrained Lagrangian (full-space method)
   * - ``tikhonov_gn``
     - PISM (custom)
     - Gauss-Newton with adaptive Tikhonov parameter selection
   * - ``sd``, ``nlcg``, ``ign``
     - siple
     - Iterative gradient methods (steepest descent, NLCG, inexact GN);
       require the optional ``siple`` Python package

The TAO-based methods are implemented through the ``IPTaoTikhonovProblem``
class template (``include/pism/inverse/IPTaoTikhonovProblem.hh``), which
provides PETSc with `\mathcal{J}(\zeta)` and `\nabla \mathcal{J}` callbacks.

.. _sec-inv-ssa-impl:

Implementation
--------------

Key files:

- ``include/pism/inverse/IP_SSATaucForwardProblem.hh`` and
  ``src/inverse/IP_SSATaucForwardProblem.cc`` implement the forward problem
  and the adjoint (`apply_jacobian_design`,
  `apply_jacobian_design_transpose`, `apply_linearization`,
  `apply_linearization_transpose`).
- ``IPDesignVariableParameterization.hh/cc`` implements the parameterizations
  `g(\zeta)`.
- ``IPTaoTikhonovProblem.hh`` is a class template that wraps any forward
  problem into a TAO-compatible Tikhonov minimization.
- ``IP_SSATaucTaoTikhonovProblem.hh/cc`` is the concrete instantiation for
  tauc inversion (with bound constraints when `tikhonov_blmvm` is used).
- The Python user interface is in ``site-packages/PISM/invert/ssa.py`` and
  the driver script is ``examples/inverse/pismi_ssa.py``.

A second forward problem ``IP_SSAHardavForwardProblem`` provides analogous
functionality for inverting the vertically-averaged ice hardness
`\bar B` rather than `\tau_c`.

.. _sec-inv-ssa-references:

References
----------

The PISM SSA inversion framework is based on the work of David Maxwell
(:cite:`Maxwelletal2008`) and was applied to Greenland outlet glaciers
in :cite:`Habermannetal2013` and :cite:`Habermannetal2017`. Earlier ice-sheet adjoint inversion work
includes :cite:`MacAyealtutorial` and :cite:`Goldberg2011`.
