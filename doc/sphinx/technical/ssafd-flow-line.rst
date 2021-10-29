.. include:: ../global.txt

.. _sec-flow-line-ssa:

Notes about the flow-line SSA
=============================

Using the same notation as in the rest of the manual, the flow-line case the shallow shelf
approximation reads

.. math::
   :label: eq-ssa-flow-line

   \left(4 \nu H u' \right)' = \rho_i\, g\, H\, h'.

Here `u` is the ice velocity, `\nu = \nu(u')` is the ice viscosity, `H` is the ice
thickness, and `h` is the ice surface elevation.

Let `\alpha = (1 - \rhoi / \rhow)` and note that `h = \alpha\, H` whenever ice is
floating, so

.. math::
   :label: eq-ssa-flow-line-2

   \left(4 \bar \nu H u' \right)' = \frac12 \alpha\, \rho_i\, g\, (H^2)'.

We can easily integrate this, getting

.. math::
   :label: eq-ssa-flow-line-3

   4 \bar \nu H u' = \frac12 \alpha\, \rho_i\, g\, H^2.

.. note::

   This is a non-linear *first-order* ODE for the ice velocity `u`.

.. _sec-ssa-flow-line-invariance:

An observation
--------------

Consider a boundary-value problem for the velocity `u` of a floating ice shelf on an
interval `[0, L]`:

.. math::
   :label: eq-ssa-flowline-bvp

   \left(4 \nu H u' \right)' &= \rho_i\, g\, H\, h',

   u(0) &= u_{0},

   \left. \left( 4\nu H u' \right) \right|_{x=L} &= \tau_{\text{stat}}(L),

   \tau_{\text{stat}}(x) &= \frac12 \rhoi\, g\, \alpha\, H(x)^2.

Now consider a similar BVP on `[0, a]` for some positive `a < L`:

.. math::
   :label: eq-ssa-flowline-bvp-2

   \left(4 \nu H v' \right)' &= \rho_i\, g\, H\, h',

   v(0) &= u_{0},

   \left. \left( 4\nu(v') H v' \right) \right|_{x=a} &= \tau_{\text{stat}}(a).

Because :eq:`eq-ssa-flowline-bvp` is a first-order ODE, the solutions `u` of
:eq:`eq-ssa-flowline-bvp` and `v` of :eq:`eq-ssa-flowline-bvp-2` *coincide* on `[0, a]`,
i.e. `u(x) = v(x)` for all `x \in [0, a]`.

.. note::

   This implies that the velocity `u(x)` of a floating flow-line ice shelf modeled by
   :eq:`eq-ssa-flowline-bvp` is not sensitive to ice geometry perturbations at locations
   downstream from `x`.

.. _sec-ssa-flow-line-discretized:

Discrete analog of this property
--------------------------------

Let `N > 2` be an integer and define the `N`\-point grid `x_{j} = 0 + j\dx`, `\dx =
L/(N-1)`.

Discretizing :eq:`eq-ssa-flowline-bvp` on this grid results in a non-linear system
with `N-1` unknowns (call this system `S_{L}`).

Let `a = L - \dx` and write down a system of equations by discretizing
:eq:`eq-ssa-flowline-bvp-2` on the grid consisting of the first `N-1` of `x_{j}`. This
system (call it `S_{a}`) has `N-2` unknowns.

The property we would like our discretization to have is this:

   If `u_{1},\dots,u_{N-1}` is the solution of `S_{L}` and `v_{1},\dots,v_{N-2}` solves
   `S_{a}`, then `u_{k} = v_{k}` for all `k = 1,\dots,N-2`.

Note that the first `N-3` equations in `S_{L}` and `S_{a}` are the same.

.. _sec-ssa-flow-line-discrete:

Discretization
--------------

Let `[X]_{j}` be an approximation of `X` at a grid location `j`.

The standard approach *within* the domain is to use centered finite differences and linear
interpolation to approximate staggered-grid values, i.e. averaging values at immediate
regular grid point neighbors:

.. math::
   :label: eq-ssa-flow-line-fd-basics

   {}[f]_{i + \frac12} &= \frac12 (f_i + f_{i+1}),

   {}[f']_{i} &= \frac1{\dx}([f]_{i+\frac12} + [f]_{i-\frac12}).

.. _sec-ssafd-flow-line-generic:

Interior
~~~~~~~~

.. math::
   :label: eq-ssa-flow-line-interior

   \frac1{\dx} \left( [4\nu(u')H u']_{i+\frac12} - [4\nu(u')H u']_{i-\frac12} \right)
   = [\rhoi\, g\, \alpha\, H\, H']_{i}.

.. _sec-ssa-flow-line-ice-front:

Ice front
~~~~~~~~~

Let `n` be the last grid point with non-zero ice thickness, i.e. assume that `H_{k} = 0`
for `k > n`.

The implementation of the stress boundary condition at the ice front amounts to adding one
more equation (see :eq:`ssafd-cfbc-vertintbdry`):

.. math::
   :label: eq-ssa-flow-line-cfbc

   {}[4\nu(u') H u']_{n+\frac12} = [\tau_{\text{stat}}]_{n+\frac12}.

We can then combine :eq:`eq-ssa-flow-line-cfbc` with :eq:`eq-ssa-flow-line-interior` (with
`i` replaced by `n`) to get the discretization at the ice front:

.. math::
   :label: eq-ssa-flow-line-cf

   \frac1{\dx}\left( [\tau_{\text{stat}}]_{n+\frac12} - [4\nu(u')H u']_{n-\frac12} \right)
   &= \rhoi\, g\, \alpha\, [H\, H']_{n},

   - [4\nu(u')H u']_{n-\frac12}
   &= \dx\rhoi\, g\, \alpha\, [H\, H']_{n} - [\tau_{\text{stat}}]_{n+\frac12}.

.. _sec-ssa-flow-line-fd:

Choosing FD approximations
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rubric:: Assuming that the ice front is at `n`

If we assume that the ice front is at `n`, the last equation in the system looks like
:eq:`eq-ssa-flow-line-cf`:

.. math::
   :label: eq-ssa-flow-line-cf-n

   - [4\nu(u')H u']_{n-\frac12}
   = \dx\rhoi\, g\, \alpha\, [H\, H']_{n} - [\tau_{\text{stat}}]_{n+\frac12}.

.. rubric:: Assuming that the ice front is at `n+1`

If we assume that the ice front is at `n+1`, the `n`\-th equation looks like a "generic"
interior equation :eq:`eq-ssa-flow-line-interior` and we have one more equation
(:eq:`eq-ssa-flow-line-cf-n`) shifted by `1`:

.. math::
   :label: eq-ssa-flow-line-cf-n-plus-1

   \frac1{\dx} \left( [4\nu(u')H u']_{n+\frac12} - [4\nu(u')H u']_{n-\frac12} \right)
   &= [\rhoi\, g\, \alpha\, H\, H']_{n}.

   - [4\nu(u')H u']_{(n+1)-\frac12}
   &= \dx\rhoi\, g\, \alpha\, [H\, H']_{n+1} - [\tau_{\text{stat}}]_{n+1+\frac12}.

Note that the second equation in :eq:`eq-ssa-flow-line-cf-n-plus-1` is the same as
:eq:`eq-ssa-flow-line-cf-n`, but with the index shifted by `1`. Both correspond to
locations *at the ice front*.

.. rubric:: The goal

We want to choose FD approximations `[\tau_{\text{stat}}]_{*}` and `[H\, H']_{*}`
in a way that would make it possible to obtain :eq:`eq-ssa-flow-line-cf-n` by
transforming equations :eq:`eq-ssa-flow-line-cf-n-plus-1`.

We propose using *constant extrapolation* to approximate `H_{n+\frac12}`

.. math::
   :label: eq-ssa-flow-line-fd-approx

   {}[H]_{i+\frac12} &=
   \begin{cases}
   H_{i}, &i \text{ is at the ice front}\\
   \frac12(H_{i} + H_{i+1}), &\text{otherwise}.
   \end{cases}

This gives us the following approximation of derivatives:

.. math::

   {}[H']_{i} &=
   \begin{cases}
   \frac1{\dx}(H_{i} - [H]_{i-\frac12}), &i \text{ is at the ice front}\\
   \frac1{\dx}([H]_{i+\frac12} - [H]_{i-\frac12}), &\text{otherwise}.
   \end{cases}

After substituting :eq:`eq-ssa-flow-line-fd-basics` this becomes

.. math::
   :label: eq-ssa-flow-line-fd-approx-derivative

   {}[H']_{i} &=
   \begin{cases}
   \frac1{2\dx}(H_{i} - H_{i-1}), &i \text{ is at the ice front}\\
   \frac1{2\dx}(H_{i+1} - H_{i-1}), &\text{otherwise}.
   \end{cases}

.. note::

   The ice front case in :eq:`eq-ssa-flow-line-fd-approx-derivative` is the **one half**
   of the standard one-sided finite-difference approximation of `H'`.

.. rubric:: Checking if :eq:`eq-ssa-flow-line-fd-approx` is the right choice

Consider the first equation in :eq:`eq-ssa-flow-line-cf-n-plus-1` and note that it
corresponds to the case in which `n` is **not** at the ice front.

.. math::

   \frac1{\dx} \left( [4\nu(u')H u']_{n+\frac12} - [4\nu(u')H u']_{n-\frac12} \right)
   = [\rhoi\, g\, \alpha\, H\, H']_{n},

Multiplying by `\dx` and moving one of the terms to the right hand side, we get

.. math::
   :label: eq-ssa-flow-line-cf-n-plus-1-first

   - [4\nu(u')H u']_{n-\frac12}
   &= \dx[\rhoi\, g\, \alpha\, H\, H']_{n} - [4\nu(u')H u']_{n+\frac12}.

   - [4\nu(u')H u']_{n-\frac12}
   &= \dx\, \rhoi\, g\, \alpha\, H_{n}\frac{H_{n+1} - H_{n-1}}{2\dx} - [4\nu(u')H u']_{n+\frac12}.

   &= \frac12\, \rhoi\, g\, \alpha\, H_{n}\, (H_{n+1} - H_{n-1}) - [4\nu(u')H u']_{n+\frac12}.

Now consider the second equation in :eq:`eq-ssa-flow-line-cf-n-plus-1`. Note that here
`n+1` **is** at the ice front.

.. math::
   :label: eq-ssa-flow-line-cf-n-plus-1-second

   - [4\nu(u')H u']_{(n+1)-\frac12}
   &= \dx\rhoi\, g\, \alpha\, [H\, H']_{n+1} - [\tau_{\text{stat}}]_{n+1+\frac12}.

   &= \dx\rhoi\, g\, \alpha\, H_{n+1}[H']_{n+1} - \frac12 \rhoi\, g\, \alpha\, [H]_{n+1+\frac12}^{2}.

   &= \dx\rhoi\, g\, \alpha\, H_{n+1}\frac{H_{n+1} - H_{n}}{2\dx} - \frac12 \rhoi\, g\,
   \alpha\, H_{n+1}^{2}

   &= \frac12\, \rhoi\, g\, \alpha( H_{n+1}(H_{n+1} - H_{n}) - H_{n+1}^{2} )

   &= \frac12\, \rhoi\, g\, \alpha( H_{n+1}^{2} - H_{n+1}H_{n} - H_{n+1}^{2} )

   &= -\frac12\, \rhoi\, g\, \alpha\, H_{n+1}H_{n}.

Put together, :eq:`eq-ssa-flow-line-cf-n-plus-1-first` and
:eq:`eq-ssa-flow-line-cf-n-plus-1-second` read as follows:

.. math::
   :label: eq-ssa-flow-line-cf-n-plus-1-rewritten

   - [4\nu(u')H u']_{n-\frac12}
   &= \frac12\, \rhoi\, g\, \alpha\, H_{n}\, (H_{n+1} - H_{n-1}) - [4\nu(u')H u']_{n+\frac12},

   - [4\nu(u')H u']_{n+\frac12}
   &= -\frac12\, \rhoi\, g\, \alpha\, H_{n+1}H_{n}.

Substituting the second equation into the first produces

.. math::
   :label: eq-ssa-flow-line-cf-n-plus-1-final

   - [4\nu(u')H u']_{n-\frac12}
   &= \frac12\, \rhoi\, g\, \alpha\, H_{n}\, (H_{n+1} - H_{n-1}) - \frac12\, \rhoi\, g\, \alpha\, H_{n+1}H_{n}

   &= \frac12\, \rhoi\, g\, \alpha\, (H_{n}\, (H_{n+1} - H_{n-1}) - H_{n+1}H_{n})

   &= -\frac12\, \rhoi\, g\, \alpha\, H_{n}\, H_{n-1}

Compare :eq:`eq-ssa-flow-line-cf-n-plus-1-final` to
:eq:`eq-ssa-flow-line-cf-n-plus-1-second`. Note that *they are the same*, except for the
index shift. In other words, :eq:`eq-ssa-flow-line-cf-n-plus-1-final` is the same as
:eq:`eq-ssa-flow-line-cf-n`, as desired.

This confirms that finite difference approximations :eq:`eq-ssa-flow-line-fd-approx` and
:eq:`eq-ssa-flow-line-fd-approx-derivative` result in a discretization with the property
we seek:

   modeled ice velocity at a given location `x` along a flow-line ice shelf is not
   sensitive to geometry perturbations downstream from `x`.
