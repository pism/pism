.. include:: ../global.txt

.. _sec-steady-hydro:

Computing steady-state subglacial water flux
============================================

The "routing" subglacial hydrology model is described by equations

.. math::
   :label: eq-hydrology-routing

   \diff{W}{t} + \diff{\Wtill}{t} + \Div \bq = \frac{m}{\rho_w}

   \diff{\Wtill}{t} = \frac{m}{\rho_w} - C_d

   \bq = -k W^{\alpha} |\nabla \psi|^{\beta - 2} \nabla \psi

   \psi = P_o + \rho_w g (b + W)

on a an ice covered area `\Omega`. We assume zero flux boundary conditions on the inflow
part of the boundary and no boundary condition on the outflow boundary. See
:cite:`BuelervanPelt2015` (equations 1, 2, 6, 16, 26) for details and the notation. Here
we also assume that `m \ge 0`.

Our goal is to estimate `Q = \bq \cdot \n`, the flux through the outflow part of the
boundary of `\Omega` corresponding to the steady state of :eq:`eq-hydrology-routing` using
a method that is computationally cheaper than using the explicit in time approximation
of :eq:`eq-hydrology-routing` described by :cite:`BuelervanPelt2015`.

Pick a contiguous section `\omega` of `\partial \Omega` (the terminus of an outlet
glacier, for example). Let `B` be the union of all the trajectories of the vector field
`\bq` in `\Omega` that pass through `\omega`. The area `B` is the "drainage basin"
corresponding to `\omega`.

Let `\gamma = \partial B \setminus \omega`. Note that if a point `P` is in `\gamma` then
one of the following conditions is satisfied.

#. `|\bq| = 0` (it is the origin of a trajectory that starts within `\Omega`) or
#. `P \in \partial \Omega` (specifically, `P` is a part of the inflow part of the boundary
   of `\Omega`)
#. `\bq \cdot \n = 0` (`P` is not at the end of a trajectory, and so the normal to the
   boundary is orthogonal to `\bq`).

Therefore `\bq \cdot \n = 0` on `\gamma` and

.. math::

   \oint_{\partial B} \bq \cdot \n\; ds &= \int_{\omega} \bq \cdot \n\; ds + \int_{\gamma} \bq \cdot \n\; ds

   &= \int_{\omega} \bq \cdot \n ds.

Assuming the steady state (and setting time derivatives in :eq:`eq-hydrology-routing` to
zero), integrating over `B`, and applying the divergence theorem gives

.. math::
   :label: eq-steady-hydro-1

   \int_{\omega} \bq \cdot \n\; ds = \int_{B} \frac{m}{\rho_w},

i.e. *in a steady state the flux through a terminus is equal to the total rate at which
water is added to the corresponding drainage basin due to the source term*.

Next, consider a related initial boundary value problem

.. math::
   :label: eq-emptying-problem

   \diff{u}{t} = -\Div (\V u)

on `B` with `u(x, y, 0) = u_0(x, y)` (`u_0 \ge 0`), `\V = -k(x, y) \nabla \psi`, zero flux
on the inflow boundary, and no boundary condition on the outflow boundary.

Here `\psi` is the hydraulic potential corresponding to the steady state of
:eq:`eq-hydrology-routing` and `k(x, y)` is a strictly positive but otherwise arbitrary
conductivity function.

Note that since `\psi` is a steady state hydraulic potential all trajectories of the
vector field `\V` leave `B` and for `\epsilon > 0` there is a time `T > 0` such that

.. math::

   \int_B u(T) = \epsilon \int_B u_0.

Integrating over time from `0` to `T`, we get

.. math::

   \int_0^T \diff{u}{t}\, dt = - \int_0^T \Div (\V u),\, \text{or}

.. math::

   u_0 = u(T) + \int_0^T \Div (\V u).

Integrating over `B` and using the divergence theorem gives

.. math::

   \int_B u_0 &= \int_B u(T) + \int_B \int_0^T \Div (\V u)

   &= \epsilon \int_B u_0 + \int_0^T \int_B \Div (\V u)

   &= \epsilon \int_B u_0 + \int_0^T \oint_{\partial B} (\V u) \cdot \n

   &= \epsilon \int_B u_0 + \int_0^T \int_{\omega} (\V u) \cdot \n.

Finally,

.. math::
   :label: eq-steady-hydro-2

   \int_B u_0 = \frac{1}{1 - \epsilon} \int_0^T \int_{\omega} (\V u) \cdot \n

Combining :eq:`eq-steady-hydro-1` and :eq:`eq-steady-hydro-2` and choosing `u_0 = \tau m\,
/\, \rho_w` for some `\tau > 0`\ [#]_ gives us a way to estimate the flux through the
outflow boundary *if we know the direction of the steady state flux*:

.. math::
   :label: eq-steady-hydro-3

   \int_{\omega} \bq \cdot \n\; ds = \frac{1}{\tau(1 - \epsilon)} \int_0^T \int_{\omega} (\V u) \cdot \n\; ds.

Here the right hand side of :eq:`eq-steady-hydro-3` can be estimated by advancing an
explicit-in-time approximation of :eq:`eq-emptying-problem` until `\int_B u` drops below
a chosen threshold.

However, the direction of the steady state flux `\bq` depends on steady state
distributions of `W` and `\Wtill` and these quantities are expensive to compute.

To avoid this issue we note that `W \ll H` and so `\psi` is well approximated by `\psi_0 =
P_o + \rho_w g b` everywhere except the vicinity of subglacial lakes. Moreover, if
`|\nabla W|` is small then `\nabla \psi_0` is a reasonable approximation of `\nabla \psi`.

We approximate `\psi` by `\tilde \psi = P_o + \rho_w g b + \delta` where `\delta > 0` is
an adjustment needed to ensure that `\tilde \psi` has no local minima in the interior of
`\Omega` and `|\nabla \tilde \psi| > 0` everywhere on `\Omega` except possibly on a set of
measure zero (no "plateaus").

The approximation of `\tilde \psi` on a computational grid is computed as follows.

1. Set `k = 0`, `\tilde \psi_{0} = \psi`.
2. Iterate over all grid points. If a grid point `(i, j)` is at a local minimum, set
   `\tilde \psi_{k + 1}(i, j)` to the average of neighboring values of `\tilde \psi_{k}`
   plus a small increment `\Delta \psi`, otherwise set `\tilde \psi_{k + 1}(i, j)` to
   `\tilde \psi_{k}(i, j)`.
3. If step 2 found no local minima, stop. Otherwise increment `k` and proceed to step 2.

Next, note that it is not necessary to identify the drainage basin `B` for a terminus
`\omega`: it is defined by `\psi` and therefore an approximation of
:eq:`eq-emptying-problem` will automatically distribute water inputs from the ice surface
(or melting) along the ice margin.

.. _sec-steady-hydro-algorithm:

The algorithm
^^^^^^^^^^^^^

Using an explicit time stepping approximation of :eq:`eq-emptying-problem` we can estimate
`\int_{\omega} \bq \cdot \n \; ds` as follows.

#. Given ice thickness `H` and bed elevation `b` compute `\tilde \psi` by filling "dips"
   as described above.

#. Choose the stopping criterion `\epsilon > 0` and the scaling for the source term `\tau > 0`.

#. Set

   .. math::

      u &\leftarrow \frac{\tau m}{\rho_w},

      t &\leftarrow 0,

      Q &\leftarrow (0, 0).

#. Compute the CFL time step `\dt` using `u` and `\V`.

#. Perform an explicit step from `t` to `t + \dt`, updating `u`.

#. Accumulate this step's contribution to `Q`:

   .. math::

      Q \leftarrow Q + \dt \cdot \V u.

#. Set `t \leftarrow t + \dt`

#. If `\int_{\Omega} u\; dx\, dy > \epsilon`, go to 4.

#. Set

   .. math::

      Q \leftarrow \frac{1}{t (1 - \epsilon^{*})}\; Q,

   where

   .. math::

      \epsilon^{*} = \frac{\int_{\Omega} u}{\int_{\Omega} u_{0}}.

.. rubric:: Footnotes

.. [#] The constant `\tau` is needed to get appropriate units, but its value is irrelevant.
