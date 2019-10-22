.. include:: ../global.txt
.. include:: ../math-definitions.txt

.. _sec-steady-hydro:

Computing steady-state subglacial water flux
============================================

Consider the simplified mass conservation equation omitting water diffusion

.. math::
   :label: eq-water-conservation

   \diff{W}{t} + \div (\V W) = M,

defined on an ice covered area `A`. Here `W` is water thickness in meters, `\V` is the
water velocity (independent of `W`), and `M(x,y) \ge 0` is the input rate in meters per
second. We assume zero flux boundary condition on the inflow part of the boundary and no
boundary condition on the outflow part of the boundary.

Our goal is to estimate the flux `Q` through the outflow part of the boundary
corresponding to the steady state distribution of `W`:

.. math::
   :label: eq-steady-flux

    Q = \V W \cdot \n.

We assume that `\V` is a gradient of a continuously differentiable hydraulic potential
`\psi`

.. math::
   :label: eq-water-velocity

   \V = - k\; \nabla \psi

and `k` is the hydraulic conductivity. We assume that `\psi` has no local minima in `A`.

If :eq:`eq-water-conservation` has a steady state, we have

.. math::
   :label: eq-water-steady-state

   \div (\V W) = M.

Pick a contiguous section `\omega` of `\partial A` (say, at the terminus of an outlet
glacier). Let `B` be the union of all trajectories of the vector field `\V` in `A` that
pass through `\omega`. `B` is the "drainage basin" of an outlet glacier for which `\omega`
is a part of the terminus.

Let `\delta = \partial B \setminus \omega`. Note that if a point `P \in \delta`, then one
of the following conditions is satisfied:

#. `|\V| = 0` (it is the origin of a trajectory that starts within `A`) or
#. `P \in \partial A` (specifically, `P` is a part of the inflow part of the boundary of
   `A`)
#. `\V \cdot \n = 0` (`P` is not at the end of a trajectory, and so the normal to the
   boundary is orthogonal to `\V`).

Now

.. math::
   :label: eq-int-over-terminus

   \oint_{\partial B} \V W \cdot \n &= \int_{\omega} \V W \cdot \n \; ds + \int_{\delta} \V W \cdot \n \; ds

   &= \int_{\omega} \V W \cdot \n \; ds

because the second term on the right hand side vanishes since `W \V \cdot \n = 0` on `\delta`.

So, the flux through the terminus (using the divergence theorem and
:eq:`eq-int-over-terminus`) is

.. math::
   :label: eq-terminus-flux

   \int_{\omega} \V W \cdot \n \; ds = \int_{B} M.

.. _sec-aux-ibvp:

The auxiliary IBVP
^^^^^^^^^^^^^^^^^^

Next, consider a related IBVP

.. math::
   :label: eq-dump-and-wait

   \diff{u}{t} = -\div (\V u)

on `B` with the initial state `u(x, y, 0) = u_0(x, y)`. Pick a time `T` such that `\int_B
u(T) = \epsilon \int_B u_0` for some `0 < \epsilon < 1`. Then

.. math::
   :label: eq-time-integral

   \int_0^T \diff{u}{t} dt &= -\int_0^T \div (\V u)

   u(T) - u(0) &= -\int_0^T \div (\V u)

   u_0 &= u(T) + \int_0^T \div (\V u).

Using the divergence theorem again, we get

.. math::

   \int_B u_0 &= \int_B \left(u(T) + \int_0^T \div (\V u) \right)

   &= \epsilon \int_B u_0 + \int_0^T \int_B \div (\V u)

   &= \epsilon \int_B u_0 + \int_0^T \oint_{\partial B} (\V u) \cdot \n

   &= \epsilon \int_B u_0 + \int_0^T \int_{\omega} (\V u) \cdot \n \; ds.

This gives

.. math::
   :label: eq-time-and-space-integral

   (1 - \epsilon) \int_B u_0 &= \int_0^T \int_{\omega} (\V u) \cdot \n \; ds,

   \int_B u_0 &= \frac{1}{1 - \epsilon} \int_0^T \int_{\omega} (\V u) \cdot \n \; ds.

Combining :eq:`eq-time-and-space-integral` and :eq:`eq-terminus-flux` and setting
`u_0 = M \tau` for some constant time `\tau > 0` we get [#]_

.. math::
   :label: eq-final

   \int_{\omega} \V W \cdot \n \; ds = \frac{1}{\tau (1 - \epsilon)}\int_0^T \int_{\omega} (\V u) \cdot \n \; ds

.. _sec-steady-flux:

The algorithm
^^^^^^^^^^^^^

Using an explicit time stepping approximation of :eq:`eq-dump-and-wait` we can estimate
`\int_{\omega} \V W \cdot \n \; ds` as follows.

#. Choose the stopping criterion `\epsilon > 0`.

#. Set

   .. math::

      u &\leftarrow u_0,

      t &\leftarrow 0,

      Q &\leftarrow 0.

#. Compute the CFL time step `\Delta t` using `u` and `\V`.

#. Perform an explicit step from `t` to `t + \Delta t`, updating `u`.

#. Accumulate this step's contribution to `Q`:

   .. math::

      Q \leftarrow Q + \Delta t \cdot \Delta x \cdot \V u.

#. Set `t \leftarrow t + \Delta t`

#. If `\int_{A} u\, dx dy > \epsilon`, go to 3.

#. Set

   .. math::

      Q \leftarrow \frac{1}{t (1 - \epsilon)}\; Q

.. [#] The constant `\tau` is needed to get appropriate units, but its value is irrelevant.
