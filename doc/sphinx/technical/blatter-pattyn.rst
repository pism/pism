.. include:: ../global.txt

.. math::

.. _sec-blatter-details:

Blatter stress balance solver: technical details
================================================

.. contents::

.. _sec-bp-notation:

Notation
--------

.. list-table::

   * - `u`
     - `x`-component of ice velocity
   * - `v`
     - `y`-component of ice velocity
   * - `\uu`
     - 2D vector `(u, v)`
   * - `\n`
     - outward-pointing normal vector at domain boundaries
   * - `\eta`
     - ice viscosity (see :eq:`eq-bp-viscosity`)
   * - `\rho`
     - ice density
   * - `\rho_w`
     - water (usually sea water) density
   * - `g`
     - gravitational acceleration
   * - `z_{\text{sea level}}`
     - sea level elevation
   * - `H`
     - ice thickness
   * - `s`
     - ice surface elevation
   * - `B`
     - ice hardness
   * - `A`
     - ice softness
   * - `\varepsilon_0`
     - regularization parameter
   * - `\beta`
     - basal resistance coefficient
   * - `n`
     - Glen exponent

.. _sec-bp-intro:

Introduction
------------

This implementation is based on the PETSc example ``snes/tutorials/ex48.c`` (see
:cite:`BrownSmithAhmadia2013`) and is inspired by :cite:`Tezaur2015`,
:cite:`Lipscomb2019`, :cite:`Hoffman2018`, :cite:`Tuminaro2016`, and :cite:`Perego2012`.

Define

.. math::
   \E_1 &= \left( 2 u_x + v_y, \quad \frac12 (u_y + v_x), \quad \frac12 u_z \right),

   \E_2 &= \left( \frac12 (u_y + v_x), \quad u_x + 2 v_y, \quad \frac12 v_z \right),

   \E &= (\E_1, \E_2).

Using this notation, the Blatter-Pattyn (BP) stress balance equations are

.. math::
   :label: eq-bp

   - \nabla \cdot (2\, \eta\, \E_1) + \rho\, g\, s_x &= 0,

   - \nabla \cdot (2\, \eta\, \E_2) + \rho\, g\, s_y &= 0.

Here "`\nabla\cdot`" is the three-dimensional divergence operator and the regularized ice
viscosity `\eta` is defined by

.. math::
   :label: eq-bp-viscosity

   \eta &= \frac{B}{2} \p{\gamma_{\text{BP}} + \frac{\varepsilon_{0}}{2}}^{\exponent},

   \gamma_{\text{BP}} &= u_x^{2} + v_y^{2} + u_x v_y + \frac14 (u_y + v_x)^{2} + \frac14 u_z^{2} + \frac14 v_z^{2},

where `\gamma_{\text{BP}}` is the Blatter-Pattyn approximation of the second invariant of
the strain rate tensor:

.. math::
   :label: eq-bp-strain-rate-tensor

   \dot\epsilon_{ij} &= \frac12\left( \diff{u_i}{x_j} + \diff{u_j}{x_i} \right)

   \gamma &= \frac12 \left( \mathop{trace}(\dot\epsilon^2) - \mathop{trace}(\dot \epsilon)^2 \right).

It is also assumed that

- ice is incompressible and `\mathop{trace}(\dot \epsilon) = 0`, and
- `\diff{w}{x} \ll \diff{u}{x}` and `\diff{w}{z} \ll \diff{u}{y}`.

The BP approximation of the second invariant `\gamma_{\text{BP}}` is obtained by omitting
`w_x` and `w_y` and expressing `w_z` using incompressibility.

.. note::

   There are at least three equivalent expressions for ice viscosity in the literature:
   the form in :eq:`eq-bp-viscosity` appears in :cite:`BrownSmithAhmadia2013` while
   :cite:`Tezaur2015` and :cite:`Lipscomb2019` define the effective strain rate
   `\dot\epsilon_e`:

   .. math::
      \newcommand{\edot}{\dot\epsilon}
      \edot_e^2 = \edot_{xx}^2 + \edot_{yy}^2 + \edot_{xx}\edot_{yy} + \edot_{xy}^2 +
      \edot_{xz}^2 + \edot_{yz}^2.

   Meanwhile :cite:`Dukowiczetal2010` have

   .. math::

      \dot\epsilon_{\text{BP}}^2 &= \left( \diff{u}{x} \right)^2 +
      \left( \diff{v}{y} \right)^2 +
      \left( \diff{u}{x} + \diff{v}{y} \right)^2

      &+
      \frac12 \left( \diff{u}{y} + \diff{v}{x} \right)^2 +
      \frac12 \left( \diff{u}{z} \right)^2 +
      \frac12 \left( \diff{v}{z} \right)^{2}.

.. _sec-bp-weak-form:

Weak form
---------

.. note::

   Recall the product rule

   .. math::
      \nabla \cdot (f X) &= \nabla f \cdot X + f \nabla \cdot X, \text{ or}

      f\nabla \cdot X &= \nabla\cdot(f X) - \nabla f \cdot X

   and the divergence theorem:

   .. math::
      \int_{\Omega} \nabla\cdot (f X) = \int_{\partial \Omega} (f X)\cdot \n\, ds.

.. Below
   \Id: integral over the domain
   \Ib: integral over the boundary
   \f: body force

We omit discussions of function spaces; see :cite:`Tezaur2015`, :cite:`Perego2012` and
other references mentioned above for details.

To obtain the weak form, we multiply both equations in :eq:`eq-bp` by a *scalar* test
function `\psi` and integrate over the domain `\Omega`. For the first equation, we get

.. math::
   :label: eq-bp-weak-form-x

   \Id \psi \left( - \nabla\cdot (2\, \eta\, \E_1) + \rho\, g\, s_x \right) &= 0,

   -\Id \psi \nabla\cdot\left( 2\, \eta\, \E_1 \right) +  \Id \psi\, \rho\, g\, s_x &= 0,

   -\Id \left( \nabla \cdot \left( \psi\, 2\, \eta\, \E_1 \right) - \nabla \psi \cdot 2\, \eta\, \E_1 \right) + \Id \psi\, \rho\, g\, s_x &= 0,

   -\Id \nabla \cdot \left(\psi\, 2\, \eta\, \E_1\right) + \Id \nabla \psi \cdot 2\, \eta\, \E_1 + \Id \psi\, \rho\, g\, s_x &= 0,

   {-\Ib \left(\psi\, 2\, \eta\, \E_1\right) \cdot \n\, ds} + {\Id \nabla \psi \cdot 2\, \eta\, \E_1} + {\Id \psi\, \rho\, g\, s_x} &= 0.

Similarly, multiplying the second equation by `\psi` and integrating by parts yields

.. math::
   :label: eq-bp-weak-form-y

   {-\Ib \left(\psi\, 2\, \eta\, \E_2\right) \cdot \n\, ds} + {\Id \nabla \psi \cdot 2\, \eta\, \E_2} + {\Id \psi\, \rho\, g\, s_y} = 0.

We combine these and say that the weak form of :eq:`eq-bp` is

.. math::
   :label: eq-bp-weak-form

   {-\Ib \left(\psi\, 2\, \eta\, \E\right) \cdot \n\, ds} + {\Id \nabla \psi \cdot 2\, \eta\, \E} + {\Id \psi\, \rho\, g\, \nabla s} = 0,

where `\nabla s = ( s_x, s_y )`.

The first term corresponds to :ref:`Neumann and Robin boundary conditions <sec-bp-bc>`; it
vanishes for "natural" BC `\left(2\, \eta\, \E \right) \cdot \n = 0`. In the basal and
lateral cases this stress is nonzero. The third one corresponds to the gravitational
driving stress and is replaced by a compensatory term in :ref:`verification tests
<sec-bp-testing-verification>` that use manufactured solutions.

.. _sec-bp-bc:

Boundary conditions
-------------------

The domain boundary consists of three disjoint parts:

1. The interface between ice and the underlying bed or (if the ice is floating) water.

   At this interface PISM uses Robin BC implementing basal sliding.
2. The interface between ice and the air\ [#f1]_ above it.

   The integral over this part of the boundary is *zero* because we assume that natural
   ("no stress") boundary conditions apply. (We ignore atmospheric pressure.)

   Vertical "cliffs" at grounded margins are interpreted as approximations of the very
   steep, but not vertical, upper surface. Following this interpretation we use natural
   boundary conditions at grounded lateral margins.

3. The *vertical* interface between ice and air (above sea level) or water (below sea
   level) at marine ice margins.

   At this interface PISM uses Neumann BC corresponding to the difference between the
   cryostatic pressure of the ice on one side and the hydrostatic pressure of water on the
   other side of the interface.

In addition to this we support Dirichlet boundary conditions for :ref:`verification
<sec-bp-testing-verification>` and to :ref:`de-couple unknowns at ice-free locations
<sec-bp-ice-extent>`.

.. note::

   Our implementation supports Dirichlet boundary conditions, but this feature is not
   exposed to the rest of PISM.

   Unlike :cite:`BrownSmithAhmadia2013` we *do not* support "no-slip" BC at the base. This
   allows us to avoid Jacobian scaling tricks they needed to achieve good multigrid
   performance.

For each node that belongs to the Dirichlet boundary we assemble "trivial" equations

.. math::

   u &= u_0,

   v &= v_0,

where `(u_0, v_0)` are given.

The implementation avoids adding contributions from adjacent elements to residual and
Jacobian entries corresponding to Dirichlet locations. Prescribed velocity values are
substituted into equations that depend on them (see section 3.2 of
:cite:`BrownSmithAhmadia2013` for details).

.. _sec-bp-bc-basal:

Basal boundary
##############

The boundary condition corresponding to basal sliding is

.. math::
   :label: eq-bp-sliding

   2\, \eta\, \E \cdot \n = - \beta \, \uu.

In the weak form :eq:`eq-bp-weak-form` this corresponds to replacing the first term:

.. math::
   :label: eq-bp-sliding-weak-form

   -\Ibase (\psi\, 2\, \eta\, \E) \cdot \n\, ds = \Ibase \beta\, \uu\, ds.

Here `\beta` has the same meaning as in :eq:`eq-sliding-linear`.

Where ice is grounded `\beta` is determined as described in :ref:`sec-basestrength`. It is
assumed to be zero (corresponding to no drag) elsewhere.

The grounding line (if present) divides bottom faces of some elements into grounded parts
that experience drag and floating parts that do not. This implementation uses a low-order
quadrature with *many* equally-spaced points (a Newton-Cotes quadrature) to integrate over
the part of the basal boundary containing the grounding line. Here `\beta` is computed at
each quadrature point, depending on whether the ice is grounded or floating at its
location. This is similar to the SEP3 parameterization described in :cite:`Seroussi2014`.

.. note::

   * It may be a good idea to implement scaling of the `\beta` coefficient according to
     the fraction of an element face that is grounded, similar to SEP1 (equation 7 in
     :cite:`Seroussi2014`). An approximation of the grounded fraction for an element
     column can be computed using existing code in PISM.

   * In general the bottom face of an element is *not* planar and *not* parallel to the
     plane `z = 0`. This means that integrals over the basal boundary should be evaluated
     using parameterizations of element faces (i.e. as surface integrals) and *not* using
     2D FEM machinery.

.. _sec-bp-bc-marine:

Marine margins
##############

We assume that marine ice margins consist of *vertical* cliffs, i.e. the outward-pointing
normal vector has the form `\n = (\cdot,\cdot,0)`.

The Neumann boundary condition at marine margins is

.. math::
   :label: eq-bp-lateral-bc

   2\, \eta\, \E \cdot \n &= p_{\text{ice}} - p_{\text{water}},\;\text{where}

   p_{\text{ice}} &= \rho\, g\, (s - z),

   p_{\text{water}} &= \rho_w\, g\, \max(z_{\text{sea level}} - z, 0).

In other words, this boundary condition corresponds to the difference between the
cryostatic pressure of the ice on one side of the interface and the hydrostatic pressure
of water on the other. The atmospheric pressure is ignored. Equation
:eq:`eq-bp-lateral-bc` is a generalization of equation 18 in :cite:`Tezaur2015`.

Just like in the implementation of the basal boundary condition near grounding lines, we
use a low order quadrature with many equally-spaced points to approximate integrals over
element faces intersected by the sea level. This *should* improve the quality of the
approximation of this boundary condition (note that the right hand side of
:eq:`eq-bp-lateral-bc` is continuous but not continuously differentiable).

.. note::

   We need to evaluate the importance of the quadrature choice described above. Does using
   depth-dependent BC matter? (We could simplify the code and use depth-averaged BC it if
   it does not.) Should we use lots of quadrature points?

.. note::

   CISM 2.1 :cite:`Lipscomb2019` uses a depth-averaged lateral boundary condition *without
   justification*.

   The lateral BC described in :cite:`Tezaur2015` is equivalent to the one described here,
   but the implementation in Albany/FELIX (and therefore in MALI) matches the one in CISM.

Implementations in CISM and MALI *do not* use this boundary condition at grounded margins
because doing so appears to produce over-estimates of the ice speed near grounded ice
margins. Our experiments show the same behavior.

.. _sec-bp-solver:

Solver implementation
---------------------

.. _sec-bp-discretization:

Discretization
##############

To create a non-linear algebraic system of equations approximating :eq:`eq-bp-weak-form`,
we create a hexahedral mesh on the domain `\Omega` and use `Q_1` Galerkin finite elements.

Let `\phi_j` be the scalar trial function associated with the node `j`, then the FE
approximation of the solution `\uu` has the form

.. math::
   :label: eq-bp-basis-expansion

   u &= \sum_j \phi_j u_j,

   v &= \sum_j \phi_j v_j.

Then the problem is

  Find `u_j`, `v_j` (`j = 1,\dots,N`) so that

  .. math::

     {-\Ib \left(\psi_i 2\, \eta\, \E\right) \cdot \n\, ds} +
     {\Id \nabla \psi_i \cdot 2\, \eta\, \E} +
     {\Id \psi_i\, \rho\, g\, \nabla s} = 0

holds for all `i = 1,\dots,N`, where `N` is the number of nodes in the mesh, subject to
the :ref:`boundary conditions <sec-bp-bc>`.

As in section 3 of :cite:`BrownSmithAhmadia2013`, we write the discretization of
:eq:`eq-bp-weak-form` as an algebraic system `F(U) = 0` with Jacobian `J(U)` and solve
this nonlinear system using Newton iterations requiring approximations of `\delta U` in

.. math::
   :label: eq-bp-newton-step

   J(U)\, \delta U = -F(U),

where

.. math::
   :label: eq-bp-residual

   F(U) &= F^1(U) + F^2(U) + F^3(U),

   F^1 &= -\Ib (\psi\, 2\, \eta\, \E) \cdot \n\, ds,

   F^2 &= \Id \nabla \psi \cdot 2\, \eta\, \E,

   F^3 &= \Id \psi\, \f.

(compare to :eq:`eq-bp-weak-form`) and `J(U)` is a square sparse matrix containing one row
per node in the mesh and at most 54 non-zero entries per row (there are `2` unknowns per
node and each node belongs to at most 8 elements forming a `3\cdot3\cdot3 = 27` node
"neighborhood").

.. _sec-bp-residual:

Residual evaluation
###################

The residual evaluation is performed in the usual manner, by iterating over all the
elements containing ice (see :ref:`sec-bp-mesh`) and adding element
contributions to the "global" residual vector.

The residual itself can be broken up into the following parts:

#. The *basal boundary* term implementing the basal drag (part of `F^1`)
#. The "main" part (`F^2`)
#. The *source term* corresponding to the driving stress (`F^3`)
#. The *top boundary* part (zero in actual simulations because we use the natural BC at the top
   boundary; can be non-zero in verification tests, part of `F^1`)
#. The *marine boundary* part implementing stress (Neumann) BC at the calving front (part of `F^1`).
#. Residual at *Dirichlet nodes*.

This decomposition makes it possible to use source terms dictated by the choice of a
manufactured solution while keeping (and testing) the rest of the code.

.. note::

   We integrate over the whole domain (`\Omega`) below (see :eq:`eq-bp-residual-ii`) for
   simplicity. In actuality each integral is over the *intersection of supports of test
   and trial functions* appearing in the integrand.

.. _sec-bp-residual-main:

Main residual contribution
%%%%%%%%%%%%%%%%%%%%%%%%%%

.. math::
   :label: eq-bp-residual-ii

   F^2_{i,u} &= \Id \nabla \psi_i \cdot 2\, \eta\,  \E_1

   &= \Id \eta\, \left(\dpsi{x}\, \p{4u_x + 2v_y} + \dpsi{y}\,\p{u_y + v_x} + \dpsi{z}\, u_z \right),

   F^2_{i,v} &= \Id \nabla \psi_i \cdot 2\, \eta\,  \E_2

   &= \Id \eta\, \left(\dpsi{y}\, \p{2u_x + 4v_y} + \dpsi{x}\,\p{u_y + v_x} + \dpsi{z}\, v_z \right).

Here `F_{i, u}` if the contribution to the `u`-component of the residual at the `i`-th node, etc.

.. _sec-bp-residual-driving-stress:

Driving stress contribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. math::
   :label: eq-bp-residual-iii

   F^3_{i, u} &= \Id \psi_i\, \rho\, g\, s_x,

   F^3_{i, v} &= \Id \psi_i\, \rho\, g\, s_y.

.. _sec-bp-residual-basal:

Basal contribution
%%%%%%%%%%%%%%%%%%

.. math::
   :label: eq-bp-residual-i-base

   F^1_{i, u} &= \Ibase \psi_i\, \beta\, u,

   F^1_{i, v} &= \Ibase \psi_i\, \beta\, v.

.. _sec-bp-residual-marine:

Marine boundary contribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. math::
   :label: eq-bp-residual-i-margin

   F^1_{i, u} &= - \Ib \psi_i (p_{\text{ice}} - p_{\text{water}})\, \n_x,

   F^1_{i, v} &= - \Ib \psi_i (p_{\text{ice}} - p_{\text{water}})\, \n_y,

where

.. math::

   p_{\text{ice}} &= \rho\, g\, (s - z),

   p_{\text{water}} &= \rho_w\, g\, \max(z_{\text{sea level}} - z, 0).

.. _sec-bp-residual-dirichlet:

Residual at Dirichlet BC locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. math::
   :label: eq-bp-residual-dirichlet

   F_{i, u} &= u - u_0,

   F_{i, v} &= v - v_0,

where `(u_0, v_0)` is the prescribed velocity.

.. _sec-bp-jacobian:

Jacobian evaluation
###################

We use an analytical (as opposed to approximated using finite differences or automatic
differentiation) Jacobian computed using formulas listed below.

Here we focus on the derivation of the Jacobian contribution corresponding to the "main"
part of the residual (`F^2`, see :eq:`eq-bp-residual` and :eq:`eq-bp-residual-ii`). The
only other non-trivial contribution comes from the basal boundary condition.

This Jacobian contribution has four parts:

.. math::
   \K uu &= \diff{F^2_{i, u}}{u_j},

   \K uv &= \diff{F^2_{i, u}}{v_j},

   \K vu &= \diff{F^2_{i, v}}{u_j},

   \K vv &= \diff{F^2_{i, v}}{v_j},

with `F^2_{\cdot, \cdot}` defined by :eq:`eq-bp-residual-ii`.

Let `G_{i, k} = \nabla\psi_i \cdot 2\, \E_{k}`, then

.. math::

   F^2_{i, u} = \Id \eta\, G_{i, 1}

and (using the product rule) we get

.. math::
   :label: eq-bp-jacobian-useful

   \K uu &= \Id \eta \diff{G_{i, 1}}{u_j} + \diff{\eta}{u_j} G_{i, 1},

   \K uv &= \Id \eta \diff{G_{i, 1}}{v_j} + \diff{\eta}{v_j} G_{i, 1},

   \K vu &= \Id \eta \diff{G_{i, 2}}{u_j} + \diff{\eta}{u_j} G_{i, 2},

   \K vv &= \Id \eta \diff{G_{i, 2}}{v_j} + \diff{\eta}{v_j} G_{i, 2}.

The derivatives of `\eta` are computed using the chain rule:

.. math::

   \diff{\eta}{u_j} &= \diff{\eta}{\gamma}\, \diff{\gamma}{u_j},

   \diff{\eta}{v_j} &= \diff{\eta}{\gamma}\, \diff{\gamma}{v_j}.

Taking derivatives of `\eta` and `\gamma` :eq:`eq-bp-viscosity` gives

.. math::
   :label: eq-bp-deta-dgamma

   \diff{\eta}{\gamma} &= \frac{B}{2} \cdot \exponent \cdot \p{\gamma + \frac{\varepsilon_{0}}{2} }^{\exponent-1}

   &= \eta \cdot \exponent \cdot \p{\gamma + \frac{\varepsilon_{0}}{2} }^{-1}.

and

.. math::
   :label: eq-bp-dgamma-duj

   \diff{\gamma}{u_j} &=
   2u_x \diff{u_x}{u_j} + 2v_y \diff{v_y}{u_j} + \diff{u_x}{u_j}v_y + u_x\diff{v_y}{u_j}
   +\frac 14 \cdot 2 \p{u_y + v_x}\p{\diff{u_y}{u_j} + \diff{v_x}{u_j} }

   &+ \frac 14 \cdot 2 u_z \diff{u_z}{u_j} + \frac 14 \cdot 2 v_z \diff{v_z}{u_j}

   &= 2u_x\dphi[j]{x} + v_y\dphi[j]{x} + \frac 12 \dphi[j]{y}\p{u_y + v_x} + \frac 12 u_z \dphi[j]{z},

   \diff{\gamma}{v_j} &=
   2u_x \diff{u_x}{v_j} + 2v_y \diff{v_y}{v_j} + \diff{u_x}{v_j}v_y + u_x\diff{v_y}{v_j}
   +\frac 14 \cdot 2 \p{u_y + v_x}\p{\diff{u_y}{v_j} + \diff{v_x}{v_j} }

   &+ \frac 14 \cdot 2 u_z \diff{u_z}{v_j} + \frac 14 \cdot 2 v_z \diff{v_z}{v_j}

   &= 2v_y\dphi[j]{y} + u_x\dphi[j]{y} + \frac 12 \dphi[j]{x}\p{u_y + v_x} + \frac 12 v_z \dphi[j]{z}.

The derivatives of `G_{\cdot, \cdot}` are

.. math::
   :label: eq-bp-dg-duj

   \diff{G_{i, 1}}{u_j} &= 4\dpsi{x}\dphi[j]{x} + \dpsi{y}\dphi[j]{y} + \dpsi{z}\dphi[j]{z},

   \diff{G_{i, 1}}{v_j} &= 2\dpsi{x}\dphi[j]{y} + \dpsi{y}\dphi[j]{x},

   \diff{G_{i, 2}}{u_j} &= 2\dpsi{y}\dphi[j]{x} + \dpsi{x}\dphi[j]{y},

   \diff{G_{i, 2}}{v_j} &= 4\dpsi{y}\dphi[j]{y} + \dpsi{x}\dphi[j]{x} + \dpsi{z}\dphi[j]{z}.

To compute `\diff{\gamma}{u_j}` :eq:`eq-bp-dgamma-duj` and `\diff{G_{\cdot,
\cdot}}{\cdot}` :eq:`eq-bp-dg-duj` we use FE basis expansions of `u` and `v`
:eq:`eq-bp-basis-expansion`, which imply:

.. math::

   \diff{u}{u_j} &= \phi_j,

   \diff{u_x}{u_j} &= \diff{\phi_j}{x},

and so on.

.. _sec-bp-jacobian-basal:

Basal contribution to the Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. math::
   :label: eq-bp-jacobian-basal

   \K uu &= \Ibase \psi_i\, \phi_j\, \left( \beta + \dbeta\, u^2 \right),

   \K uv &= \Ibase \psi_i\, \phi_j\, \dbeta\, u\, v,

   \K vu &= \Ibase \psi_i\, \phi_j\, \dbeta\, v\, u,

   \K vv &= \Ibase \psi_i\, \phi_j\, \left( \beta + \dbeta\, v^2 \right).

Here `\dbeta` is the derivative of `\beta` with respect to `\alpha = \frac12 |\uu|^2 =
\frac12 \left( u^2 + v^2 \right)` (one of the outputs of PISM's basal resistance parameterizations).

Note that

.. math::

   \diff{\alpha}{u_j} &= \frac 12 \cdot (2u) \cdot \diff{u}{u_j}

   &= u\, \phi_j,

   \diff{\alpha}{v_j} &= v\, \phi_j.

Recall (see :eq:`eq-bp-residual-i-base`) that the basal sliding contribution to the
residual is equal to

.. math::

   F^1_u &= \Ibase \psi_i\, \beta\, u,

   F^1_v &= \Ibase \psi_i\, \beta\, v,

and so the basal contribution to the Jacobian consists of partial derivatives of these
with respect to `u_j` and `v_j`.

For example,

.. math::

   \K uu &= \diff{\left( \Ibase\psi_i \beta\, u \right)}{u_j}

   &= \Ibase \psi_i\, \diff{\left( \beta\, u \right)}{u_j}

   &= \Ibase \psi_i\, \left( \beta\, \diff{u}{u_j} + \diff{\beta}{u_j}\, u \right)

   &= \Ibase \psi_i\, \left( \beta\, \phi_j + \dbeta \cdot\diff{\alpha}{u_j}\, u \right)

   &= \Ibase \psi_i\, \left( \beta\, \phi_j + \dbeta\, u^2\, \phi_j \right)

   &= \Ibase \psi_i\, \phi_j \left( \beta + \dbeta\, u^2 \right).

.. _sec-bp-jacobian-dirichlet:

Jacobian at Dirichlet locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The Jacobian at Dirichlet locations is set to `1`. Together with the :ref:`residual at
Dirichlet locations <sec-bp-residual-dirichlet>` this completes the assembly of trivial
equations at these locations.

More specifically (due to the interleaved ordering of the unknowns: `u_j`, `v_j`,
`u_{j+1}`, `v_{v+1}`, ...) this requires setting the corresponding block of the Jacobian
to the `2 \times 2` identity matrix.

.. note::

   It may be interesting to see if varying the scaling of Jacobian entries at Dirichlet
   nodes affects the condition number of the Jacobian matrix. See the variable ``scaling``
   in the code and set

   .. code-block:: bash

      -bp_ksp_view_singularvalues

   to see if it has an effect.

.. _sec-bc-iceberg-elimination:

Iceberg elimination
###################

As described in :cite:`Tuminaro2016`, isolated patches of ice with low basal resistance
*and* patches connected only via a single node (a "hinge") are problematic because the
system :eq:`eq-bp` determines ice velocity up to rigid rotations and translations.

This is not a new issue: both FD and FEM solvers of the SSA stress balance require some
form of iceberg elimination. We use a connected component labeling algorithm to identify
patches of *floating* ice. In the FEM context this requires inspecting *elements*: two
elements are considered connected if they share a boundary.

.. note::

   We could improve this mechanism by implementing a version of the method described in
   :cite:`Tuminaro2016`: instead of removing *floating ice* not connected to grounded
   ice (PISM's approach) they

   1. Identify all the connected patches of ice.
   2. Remove patches of ice which have no mesh *nodes* with

      .. math::

         \beta > \beta_{\text{threshold}}.

.. _sec-bp-pc:

Preconditioning
###############

Back in 2013 Brown et al :cite:`BrownSmithAhmadia2013` showed that multigrid can act as an
effective preconditioner for BP stress balance solvers. However, all test cases considered
in that work use periodic boundary conditions and therefore avoid all the issues related
to the moving ice margin present in complete ice sheet models. Our tests suggests that a
naive implementation assembling trivial equations at "ice free" nodes combined with
standard geometric multigrid (coarsening in all 3 directions) is not likely to succeed and
we need a different approach.

So, we use semi-coarsening in the vertical direction *even though Brown et al state that
"semi-coarsening is unattractive"*. One of the arguments against semi-coarsening is the
larger number of multigrid levels needed: semi-coarsening gives a smaller reduction in the
number of unknowns from one level to the next (factor of 2 instead of 8 in the full
multigrid approach). Our tests show that "aggressive" semi-coarsening (i.e. using
coarsening factors larger than 2 and as high as 8) appears to be effective, allowing one
to achieve similar reductions in the number of unknowns from one level to the next in a
multigrid hierarchy.

The second argument against semi-coarsening is deeper: spatial variations in the sliding
parameter `\beta` may lead to the kind of anisotropy that cannot be addressed by
coarsening in the vertical direction (see chapter 7 in :cite:`BriggsHensonMcCormick` for a
discussion). Still, we are encouraged by results published by Tuminaro et al
:cite:`Tuminaro2016`, who used a similar mesh structure and an approach equivalent to
using geometric multigrid with semi-coarsening in the vertical direction for the finer
part of the hierarchy and algebraic multigrid for coarser levels. [#f2]_

Inspired by :cite:`Tuminaro2016`, we use *geometric multigrid* to build a mesh hierarchy
with the coarsest level containing a small number (2 or 3) of vertical levels combined
with *algebraic* multigrid as a preconditioner for the solver on the coarsest level.
(Semi-coarsening in the vertical direction cannot reduce the number of unknowns in the `x`
and `y` directions and the coarsest problem is likely to be too large for the redundant
direct solver, which is PETSc's default.)

In addition to this, we follow :cite:`BrownSmithAhmadia2013` in ordering unknowns so that
columns are contiguous (and `u` and `v` are interleaved), allowing ILU factorization to
compute a good approximation of ice velocities in areas where SIA is applicable.

.. _sec-bp-pc-grid-coarsening:

Vertical grid sizes compatible with coarsening
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ideally, the coarsest mesh in the hierarchy should have 2 nodes in the `z` direction, i.e.
be *one element thick*. If `N` is the coarsening factor, the total number of vertical
levels (`M_z`) has to have the form

.. math::
   :label: eq-bp-mz

   M_z = A\cdot N^k + 1

for some positive integer `A` (ideally `A=1`), so that the mesh hierarchy containing `k`
levels will include levels with

.. math::

   A + 1,\, A\cdot N + 1,\, A\cdot N^2 + 1,\, \dots,\, A\cdot N^k + 1

nodes in the `z` direction.

This means that for a given :config:`stress_balance.blatter.coarsening_factor` and number
of multigrid levels ``-bp_pc_mg_levels k`` the value of `M_z` cannot be chosen at random.

.. _sec-bp-pc-options:

Controlling using PETSc options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The PETSc ``SNES`` object solving the BP system uses the ``bp_`` prefix for all
command-line options.

To choose a preconditioner, set

.. code::

   -bp_pc_type XXX

where ``XXX`` is the name of a preconditioner.

Run

.. code-block:: bash

   pismr -stress_balance blatter [other options] -help | grep "-bp"

to get the list of all PETSc options controlling this solver.

To use a geometric multigrid preconditioner with `N` levels, set

.. code-block:: bash

   -bp_pc_type mg -bp_pc_mg_levels N

An "aggressive" (i.e. greater than 2) coarsening factor may work well. Use
:config:`stress_balance.blatter.coarsening_factor` to set it.

See :ref:`sec-bp-pc-grid-coarsening` for the discussion of the relationship between the
number of vertical levels, number of multigrid levels, and the coarsening factor.

Set

.. code-block:: bash

   -bp_mg_coarse_pc_type gamg

to use PETSc's GAMG on the coarsest multigrid level.

.. note::

   It would be interesting to compare different preconditioning options on the coarsest MG
   level (GAMG, Hypre BoomerAMG, ...).

.. _sec-bp-pc-implementation:

Additional code needed to support geometric multigrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

To support the multigrid preconditioner we need to be able re-discretize the system on the
mesh provided *by PETSc* in our to residual and Jacobian evaluators. In general, this
requires

* re-computing grid-related constants (`\dx`, etc) using the grid (avoid using
  hard-wired constants, e.g. computed using the fine grid), and

* additional code to restrict gridded inputs from the fine grid mesh to coarser meshes.

This solver does not support coarsening in horizontal directions, so gridded
two-dimensional inputs can be used on all multigrid levels. The grid spacing (`\dx`,
`\dy`) remains the same as well.

To transfer the one three-dimensional gridded input field (ice hardness), we create
interpolation matrices mapping from a coarse level to the next (finer) level in the
hierarchy. The transpose of this matrix is used as a restriction operator.

.. _sec-bp-parameter-continuation:

Parameter continuation as a recovery mechanism
##############################################

As in :cite:`Tezaur2015`, we can start with a large `\varepsilon_0`, find an approximate
solution, then use it as an initial guess for the next solve with a reduced
`\varepsilon_0`.

.. note::

   Not implemented yet.

.. _sec-bp-mesh:

Model domain and mesh structure
###############################

The domain is

.. math::

   x &\in [x_{\text{min}}, x_{\text{max}}],

   y &\in [y_{\text{min}}, y_{\text{max}}],

   z &\in [z_{\text{min}}, z_{\text{min}} + H],

where `[x_{\text{min}}, x_{\text{max}}] \times [y_{\text{min}}, y_{\text{max}}]` is the
"map plane" domain corresponding to the maximum ice extent, `z_{\text{min}}` is the bottom
ice surface elevation (equal to bed elevation where ice is grounded and determined using
sea level, ice thickness, and the floatation criterion where floating) and `H` is the ice
thickness.

Coordinates of the mesh nodes have the form

.. math::
   :label: eq-bp-mesh-nodes

   x_i &= x_{\text{min}} + i \cdot \dx,

   y_j &= y_{\text{min}} + j \cdot \dy,

   z_k &= z_{\text{min}}(x_i, y_j) + k \cdot \frac{H(x_i, y_j)}{M_z - 1}.

Each element's projection onto the plane `z = 0` is a rectangle with sides `\dx` and
`\dy`, but the spacing between nodes in the vertical direction is *not* constant:
each vertical column of nodes contains `M_z` nodes with the spacing of `H / (M_z - 1)`.
This mesh structure is *exactly the same* as the one used in :cite:`BrownSmithAhmadia2013`
and CISM 2.1 :cite:`Lipscomb2019`.

.. _sec-bp-ice-extent:

Supporting evolving ice extent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The ice extent changes as the model runs and the solver implementation has to allow for
these changes.

To simplify the logic used to identify the interior of the ice volume and its lateral
boundaries we compute the "type" of each node in the mesh. Given a threshold `\Hmin` (see
:config:`stress_balance.ice_free_thickness_standard`), we say that

* an element contains ice if ice thickness at all its nodes equals to or exceeds `\Hmin`,
* a node is *interior* (within the ice) if all the elements it belongs to contain ice,
* a node is *exterior* (outside the ice volume) if no element it belongs to contains ice,
* a node that is neither interior nor exterior is a *boundary* node.

   Only elements containing ice are included in the residual and Jacobian evaluation.

We prescribe zero `(0, 0)` velocity at exterior nodes by marking them as Dirichlet BC
locations.

In addition to this, we need to identify vertical faces of elements at :ref:`lateral
boundaries <sec-bp-bc-marine>`.

   An element face is a part of the lateral boundary if all four of its nodes are
   *boundary* nodes.

.. _sec-bp-inputs-outputs:

Solver inputs and outputs
#########################

The BP solver uses the following inputs:

- basal yield stress (`\tau_c`),
- ice thickness,
- bed elevation,
- sea level elevation,
- ice enthalpy (used to compute ice hardness)

And provides the following outputs:

- `u` and `v` components of ice velocity at the :ref:`nodes of the mesh <sec-bp-mesh>`
  (saved to output files to be used as an initial guess later)
- vertically-averaged `u` and `v` (used in the mass continuity step, i.e. to update ice
  geometry)
- basal frictional heating (used to couple stress balance to energy conservation)

.. _sec-bp-computation-steps:

Steps performed by the solver
#############################

#. Compute ice bottom elevation.
#. Compute floatation function `f` (`f \le 0` if ice is grounded, `f > 0` if floating).
#. Compute :ref:`node type <sec-bp-ice-extent>`.
#. Compute ice hardness at the :ref:`nodes of the mesh <sec-bp-mesh>`.
#. Call PETSc's ``SNESSolve()``.
#. Extract basal velocity and compute basal frictional heating.
#. Compute vertically-averaged ice velocity.

.. _sec-bp-pism-integration:

Integration with the rest of PISM
#################################


.. _sec-bp-pism-energy-conservation:

Conservation of energy
%%%%%%%%%%%%%%%%%%%%%%

Coupling to PISM's energy balance models requires

- all 3 components of ice velocity on PISM's grid, and
- strain heating.

These are computed by using piecewise-linear interpolation in the vertical direction to
put `u`, `v` on PISM's grid, after which vertical velocity `w` and strain heating are
computed using existing code.

.. _sec-bp-testing-verification:

Testing and verification
------------------------

.. _sec-bp-verificaion-xy:

Verification test XY
####################

:Exact solution:
   .. math::

      u &= \exp(x) \sin(2 \pi y)

      v &= \exp(x) \cos(2 \pi y)

:Domain: `x \in [0, 1]`, `y \in [0, 1]`
:Boundary conditions: Dirichlet BC corresponding to the exact solution on all lateral
   boundaries (`x = 0`, `x = 1`, `y = 0`, `y = 1`). Natural BC at the top and bottom
   boundaries.

:Comments: This test uses a constant ice hardness, `n = 3`, and has no variation in the
   `z` direction. It is similar to the `x-y` MMS verification test in :cite:`Tezaur2015`,
   section 4.1 (we use Dirichlet boundary conditions along the whole lateral boundary
   instead of Robin conditions derived from the chosen exact solution).

   The compensatory term is computed using SymPy_; please see the code in
   ``src/stressbalance/blatter/verification``.

.. _sec-bp-verificaion-xz:

Verification test XZ
####################

:Exact solution:
   .. math::

      u &= \frac{2 A (g \rho)^{n} \left((s - z)^{n + 1} - H^{n + 1}\right) \left|s_x\right|^{n - 1} s_x}{n + 1} - \frac{H g \rho s_x}{\beta}

      v &= 0
:Domain:
   .. math::

      x &\in [-L, L],

      z &\in [s(x) - H, s(x)],

      s(x) &= s_0 - \alpha\, x^2.

:Boundary conditions:
   - Dirichlet BC corresponding to the exact solution at `x = -L` and `x = L`.
   - Basal (`z = s(x) - H`) BC is a combination (i.e. sum) of the BC in
     :ref:`sec-bp-bc-basal` and a compensatory term derived using the exact solution.
   - Top surface (`z = s(x)`) BC is derived using the exact solution.
:Comments:
   This test

   - uses a constant ice hardness,
   - Glen exponent `n = 3`,
   - has no variation (and is periodic) in the `y` direction,
   - uses a constant basal resistance coefficient `\beta`.

   It is similar to the `x-z` MMS verification test :cite:`Tezaur2015`, section 4.2
   (again, we use Dirichlet BC at lateral boundaries instead of Robin conditions stated in
   the paper).

   See :cite:`Tezaur2015` and the code in ``src/stressbalance/blatter/verification`` for
   details.

.. _sec-bp-verificaion-xz-cfbc:

Verification test XZ-CFBC
#########################

This setup tests the "calving front boundary condition" (see :ref:`sec-bp-bc-marine`).

:Exact solution:
   .. math::

      u &= \frac{(\rho - \rho_w) g L}{2 B \pi} \sin\p{\frac{\pi x}L} z,

      v &= 0
:Domain: `x \in [0, L]`, `z \in [-H, 0]`
:Boundary conditions:
   - Dirichlet BC at `x = 0`.
   - Uses the BC described in :ref:`sec-bp-bc-marine` at `x = L`.
:Comments: This test uses the Glen exponent of `1` (constant viscosity) and has no
   variation in the `y` direction.

   The sea level is set to `0`, overriding the floatation criterion to ensure that *all*
   the ice is submerged.

.. _sec-bp-verification-xz-vanderveen:

Verification test XZ-VV (van der Veen profile)
##############################################

This setup tests the implementation of the basal sliding boundary condition (see
:ref:`sec-bp-bc-basal`) using the van der Veen shelf profile :cite:`vanderVeen`:

.. math::
   :label: eq-bp-van-der-Veen

   H(x) &= \left[ \frac{4 C x}{Q_0} + H_0^{-4} \right]^{-\frac14},

   u(x) &= \frac{Q_0}{H(x)},

   v(x) &= 0,

   C &= \left( \frac{\alpha g \rho_i}{2 B} \right)^3,

where `Q_0` is the flux at the left boundary and `H_0` is the corresponding ice thickness.

The surface elevation `s` and bed elevation `b` are defined by

.. math::

   s(x) &= \alpha H(x),

   b(x) &= (\alpha - 1) H(x)

for some positive constant `\alpha`. (A free-floating shelf corresponds to `\alpha = 1 -
\rho_i / \rho_w`).

.. note::

   Functions `(u, v)` in :eq:`eq-bp-van-der-Veen` solve :eq:`eq-bp` *exactly* in the
   interior of the domain

   .. math::

      0 \le x_{\text{min}} &\le x \le x_{\text{max}},

      b(x) &\le z \le s(x)

   if the Glen exponent `n = 3`. No compensatory source term is needed.

We use a Dirichlet BC at the left boundary:

.. math::

   u(x_{\text{min}}) &= \frac{Q_0}{H_0},

   v(x_{\text{min}}) &= 0

and a stress BC at the right boundary:

.. math::
   :label: eq-bp-van-der-Veen-right-BC

   2\, \eta\, \E \cdot \n_{\text{right}} = \left(\alpha g \rho_i H(x) ,\, 0\right)

The stress BC at the top surface is

.. math::
   :label: eq-bp-van-der-Veen-top-BC

   2\, \eta\, \E \cdot \n_{\text{top}} = \left( \frac{2 B C^{\frac{4}{3}} \alpha H^{6}(x)}
   {\sqrt{C^{2} \alpha^{2} H^{10}(x) + Q_{0}^{2}}},\, 0\right)

The boundary condition at the bottom surface has the form

.. math::

   2\, \eta\, \E \cdot \n_{\text{bottom}} &= - \beta \, \uu,

   \beta &= \frac{2 B C^{\frac{4}{3}} \left(\alpha - 1\right) H^{7}(x)}
   {Q_{0} \sqrt{C^{2} \left(\alpha - 1\right)^{2} H^{10}(x) + Q_{0}^{2}}}

.. _sec-bp-known-issues:

Known issues and future work
----------------------------

- Eliminate "wiggles" near areas with steep surface slopes.
- Implement drag at lateral boundaries in fjords and alpine valleys.
- Implement parameter continuation.
- Couple to melange back pressure parameterizations by replacing `p_{\text{water}}` in
  :ref:`sec-bp-bc-marine`.

.. rubric:: Footnotes

.. [#f1] Strictly speaking the top surface of the ice may be in contact with firn or snow
   as well as air, but these details are not relevant here.

.. [#f2] The code developed by :cite:`Tuminaro2016` uses the algebraic multigrid framework
   *throughout*, i.e. even for the part of the hierarchy where the mesh structure allows
   one to use "geometric" coarsening.
