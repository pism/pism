.. include:: ../global.txt

.. only:: html

   .. When building PDFs this will be included already.

   .. include:: ../math-definitions.txt

.. _sec-bombproof:

BOMBPROOF, PISM's numerical scheme for conservation of energy
=============================================================


Introduction
------------

One of the essential goals for any thermomechanically-coupled numerical ice sheet model is
a completely bombproof numerical scheme for the advection-conduction-reaction problem for
the conservation of energy within the ice. "Bombproof" means being as stable as possible
in as many realistic modeling contexts as possible. PISM's scheme is observed to be highly
robust in practice, but it is also provably stable in a significant range of
circumstances. The scheme is special to shallow ice sheets because it makes specific
tradeoff choices with respect to vertical velocity. It is generic and low order in how it
treats horizontal velocity. In this page we state the scheme, prove its stability
properties, and address the basal boundary condition.

The scheme is conditionally-stable. The length of the time step is limited only by the
maximum magnitude of the horizontal velocities within the ice, i.e. horizontal CFL. This
condition for stability is included in the PISM adaptive time-stepping technique.

Accuracy is necessarily a second goal. Our shallow scheme has truncation error `O(\Delta
z^2)` in many circumstances, though it reverts to lower order when it "detects trouble" in
the form of large vertical velocities. Overall the scheme has order `O(\Delta t,\Delta
x,\Delta y, \Delta z^2)` in circumstances where the vertical ice flow velocity is small
enough, relative to conductivity, and otherwise reverts to
`O(\Delta t,\Delta x,\Delta y, \Delta z^1)`.

The conservation of energy problem for the ice is in terms of an enthalpy field
:cite:`AschwandenBuelerKhroulevBlatter`. The current scheme supercedes the cold-ice,
temperature-based scheme described in the appendices of :cite:`BBL` and in
:cite:`BBssasliding`. Compared to a cold-ice scheme, the enthalpy formulation does a
better job of conserving energy, has a more-physical model for basal melt rate and
drainage, and can model polythermal ice with any CTS topology
:cite:`AschwandenBuelerKhroulevBlatter`. The finite difference implementation of the
enthalpy method is robust and avoids the CFL condition on vertical advection which was
present in the older, cold-ice scheme.

The bedrock thermal problem is solved by splitting the timestep into an update of the
bedrock temperature field, assuming the ice base is as constant temperature, and an update
of the ice enthalpy field, by the BOMBPROOF scheme here, assuming the upward heat flux
from the bedrock layer is constant during the timestep. For more on the implementation,
see the ``BedThermalUnit`` class.

The region in which the conservation of energy equation needs to be solved changes over
time. This is an essential complicating factor in ice sheet modeling. Also relevant is
that the velocity field has a complicated provenance as it comes from different stress
balance equations chosen at runtime. These stress balances, especially with transitions in
flow type, for instance at grounding lines, are incompletely understood when
thermomechanically-coupled. (The ``ShallowStressBalance`` instance owned by
``IceModel`` could be the SIA, the SSA, a hybrid of these, or other stress balances in
future PISM versions.) We will therefore not make, in proving stability, assumptions about
the regularity of the velocity field in space or time other than boundedness.

Nor do we want the numerical scheme for advection to need any information about the
velocity except its value at the beginning of the time step. Thus the conservation of
energy timestep is assumed to be split from the mass continuity time step and its
associated stress balance solve. We have not considered implementing a scheme which
requires the Jacobian of the velocity field with respect to changes in enthalpy, for
example. At very least such a fully-implicit scheme would require blind iteration (e.g.
with no guarantee of convergence of the iteration). The scheme we propose involves no such
iteration.


Conservation of energy in a shallow ice sheet
---------------------------------------------

In an enthalpy formulation :cite:`AschwandenBuelerKhroulevBlatter` (and references therein),
the ice sheet is regarded as a mixture of two phases, solid and liquid, so that both cold
and temperate ice with liquid ice matrix can be modeled. The specific enthalpy field of
the ice mixture is denoted `E(x,y,z,t)` and has units `J / kg`. (Within
the PISM documentation the symbol `H` is used for ice thickness so we use `E` for enthalpy
here and in the PISM source code versus "H" in :cite:`AschwandenBuelerKhroulevBlatter`.)
The conservation of energy equation is

.. math::
   :label: basicEnergy

   \rho \frac{dE}{dt} = -\nabla \cdot \mathbf{q} + Q,

where `\rho` is the mixture density. The mixture density is assumed to be the same as ice
density even if there is a nonzero liquid fraction, and the mixture is assumed to be
incompressible :cite:`AschwandenBuelerKhroulevBlatter`. The left and right sides of
equation :eq:`basicEnergy`, and thus the quantity `Q`, have units
`J\,\text{s}^{-1}\,\text{m}^{-3} = \text{W}\,\text{m}^{-3}`.

Neglecting the dependence of conductivity and heat capacity on temperature
:cite:`AschwandenBuelerKhroulevBlatter`, the heat flux in cold ice and temperate ice is

.. math::
   :label: heatflux

   \mathbf{q} =
   \begin{cases}
      - \frac{k_i}{c_i} \nabla E, & \text{cold ice}, \\
      - K_0 \nabla E, & \text{temperate ice},\\
   \end{cases}

where `k_i,c_i,K_0` are constant :cite:`AschwandenBuelerKhroulevBlatter`. The nonzero flux
in the temperate ice case, may be conceptualized as a regularization of the "real"
equation, or as a flux of latent heat carried by liquid water. Also, `dE/dt` stands for
the material derivative of the enthalpy field, so the expanded form of :eq:`basicEnergy`
is

.. math::

   \rho \left(\frac{\partial E}{\partial t} + \mathbf{U}\cdot \nabla E\right)
   = \nabla \cdot \left(\left\{\begin{matrix} \frac{k_i}{c_i} \\ K_0 \end{matrix}\right\}
   \nabla E \right) + Q,

where `\mathbf{U}` is the three dimensional velocity, thus advection is included.

The additive quantity `Q` is the dissipation (strain-rate) heating,

.. math::

   Q = \sum_{i,j=1}^3 D_{ij} \tau_{ij}

where `D_{ij}` is the strain rate tensor and `\tau_{ij}` is the deviatoric stress tensor.
Reference :cite:`BBssasliding` addresses how this term is computed in PISM, according to
the shallow stress balance approximations; see
``StressBalance::compute_volumetric_strain_heating()``. (`Q` is called `\Sigma` in
:cite:`BBL`, :cite:`BBssasliding` and in many places in the source code.)

Friction from sliding also is a source of heating. It has units of `W / m^2 = J / (m^2 s)`, that
is, the same units as the heat flux `\mathbf{q}` above. In formulas we write

.. math::

   F_b = - \tau_b \cdot \mathbf{u}_b,

where `\tau_b` is the basal shear stress and `\mathbf{u}_b` is the basal sliding velocity;
the basal shear stress is oppositely-directed to the basal velocity. For example, in the
plastic case `\tau_b = - \tau_c \mathbf{u}_b / |\mathbf{u}_b|` where `\tau_c` is a
positive scalar, the yield stress. See method
``StressBalance::basal_frictional_heating()``. The friction heating is concentrated at
`z=0`, and it enters into the basal boundary condition and melt rate calculation,
addressed in section :ref:`sec-melt` below.

We use a shallow approximation of equation :eq:`basicEnergy` which lacks horizontal
conduction terms  :cite:`Fowler`.  For the initial analysis of the core BOMBPROOF
scheme, we specialize to cold ice.  Within cold ice, the coefficient in the heat
flux is constant, so

.. math::

   \nabla \cdot \mathbf{q} = - \frac{k_i}{c_i} \frac{\partial^2 E}{\partial z^2}.

Therefore the equation we initially analyze is

.. math::
   :label: basicShallow

   \rho_i \left(\frac{\partial E}{\partial t} + \mathbf{U}\cdot \nabla E\right) = \frac{k_i}{c_i} \frac{\partial^2 E}{\partial z^2} + Q,

We focus the analysis on the direction in which the enthalpy has largest derivative,
namely with respect to the vertical coordinate `z`.  Rewriting equation
:eq:`basicShallow` to emphasize the vertical terms we have

.. math::
   :label: vertProblem

    \rho_i \left(\frac{\partial E}{\partial t} + w \frac{\partial E}{\partial z}\right) 
         = \frac{k_i}{c_i}  \frac{\partial^2 E}{\partial z^2} + \Phi

where

.. math::

   \Phi = Q - \rho_i \left(u \frac{\partial E}{\partial x}
   + v \frac{\partial E}{\partial y}\right)

We assume that the surface enthalpy `E_s(t,x,y)` (K) and the geothermal flux `G(t,x,y)`
(`W / m^2`) at `z=0` are given. (The latter is the output of the ``energy::BedThermalUnit``
object, and it may come from an evolving temperature field within the upper crust, the
bedrock layer. If a surface temperature is given then it will be converted to enthalpy by
the ``EnthalpyConverter`` class.) The boundary conditions to problem :eq:`vertProblem`
are, therefore,

.. math::
   :label: columnbcs

   E(t,x,y,z=H) &=E_s(t,x,y),\\
   -\frac{k_i}{c_i} \frac{\partial E}{\partial z}\Big|_{z=0} &= G.

For a temperate ice base, including any ice base below which there is liquid water,
the lower boundary condition is more interesting.  It is addressed below in section
:ref:`sec-melt`.

The core BOMBPROOF scheme
-------------------------

For the discussion of the numerical scheme below, let `E_{ijk}^n` be our approximation to
the exact enthalpy `E` at the grid point with coordinates `(x_i,y_j,z_k)` at time `t_n`.
When `i,j` are uninteresting we suppress them and write `E_k^n`, and we will use similar
notation for numerical approximations to the other quantities. We put the horizontal
advection terms in the source term `\Phi` because we treat them explicitly, evaluating at
time `t_n`. (Implicit or semi-implicit treatment of horizontal advection would require a
coupled system distributed across processors, a difficulty which is currently avoided.)

The scheme we use for horizontal advection is explicit first-order upwinding. There is a
CFL condition for the scheme to be stable, in the absence of conduction, based on the
magnitude of the horizontal velocity components. To state the upwind scheme itself, let

.. math::

   \Up{f_{\bullet}}{\alpha} =
   \begin{cases}
     f_i-f_{i-1}, & \alpha \ge 0, \\
     f_{i+1}-f_i, & \alpha < 0.
   \end{cases}

The approximate horizontal advection terms, and thus the approximation to the whole
term `\Phi`, are

.. math::

   \Phi_{ijk}^n = \Sigma_{ijk}^n - \rho_i
   \left( u_{ijk}^n\,\frac{\Up{E_{\bullet jk}^n}{u_{ijk}^n}}{\Delta x}
   + v_{ijk}^n\,\frac{\Up{E_{i\bullet k}^n}{v_{ijk}^n}}{\Delta y} \right).

The CFL stability condition for this part of the scheme is


.. math::
   :label: CFL

   \Delta t \,\left( \left|\frac{u_{ijk}^n}{\Delta x}\right|
   + \left|\frac{v_{ijk}^n}{\Delta y}\right| \right) \le 1.

The routine ``max_timestep_cfl_3d()`` computes the maximum of velocity
magnitudes. This produces a time step restriction based on the above CFL condition. Then
``IceModel::max_timestep()`` implements adaptive time-stepping based on this and
other stability criteria.

In the analysis below we assume an equally-spaced grid `z_0,\dots,z_{M_z}`
with `\Delta z = z_{k+1} - z_k`.  In fact PISM has a remapping scheme in each
column, wherein the enthalpy in a column of ice is stored on an unequally-spaced 
vertical grid, but is mapped to a fine, equally-spaced grid for the conservation
of energy computation described here.  (Similar structure applies to the age
computation.  See classes ``EnthalpyModel`` and ``AgeModel``.)

The `z` derivative terms in :eq:`vertProblem` will be approximated implicitly. Let
`\lambda` be in the interval `0 \le \lambda \le 1`. Suppressing indices `i,j`, the
approximation to :eq:`vertProblem` is


.. math::
   :label: bombone

    \rho_i &\left(
           \frac{E_k^{n+1} - E_k^n}{\Delta t}
           + \lambda w_k^n \frac{E_{k+1}^{n+1} - E_{k-1}^{n+1}}{2 \Delta z}
           + (1-\lambda) w_k^{n} \frac{\Up{E_{\bullet}^{n+1}}{w_k^{n}}}{\Delta z} \right) \\
         &= \frac{k_i}{c_i}\, \frac{E_{k+1}^{n+1} - 2 E_{k}^{n+1} + E_{k-1}^{n+1}}{\Delta z^2} + \Phi_k^n.

Equation :eq:`bombone`, along with a determination of `\lambda` by
:eq:`lambdachoice` below, is the scheme BOMBPROOF.  It includes two approximations
of vertical advection,  implicit centered difference  (`\lambda = 1`) and
implicit first-order upwinding (`\lambda=0`).  They are combined using
nonnegative coefficients which sum to one, a convex combination.  The centered
formula has higher accuracy,

.. math::

   w_k^n \frac{E_{k+1}^{n+1} - E_{k-1}^{n+1}}{2 \Delta z}
   = w \frac{\partial E}{\partial z} + O(\Delta t,\Delta z^2),

while the first order upwind formula has lower accuracy,

.. math::

   w_k^{n} \frac{\Up{E_{\bullet}^{n+1}}{w_k^{n}}}{\Delta z}
   = w \frac{\partial E}{\partial z} + O(\Delta t,\Delta z).

Thus we prefer to use the centered formula when possible, but we apply (implicit)
upwinding when it is needed for its added stability benefits.

We now rewrite :eq:`bombone` for computational purposes as one of a system of equations
for the unknowns `\{E_k^{n+1}\}`.  In this system the coefficients will be
scaled so that the diagonal entries of the matrix have limit one as
`\Delta t\to 0`.  Let

.. math::

   \nu &= \frac{\Delta t}{\Delta z},\\
   R &= \frac{k_i \Delta t}{\rho_i c_i \Delta z^2}.

Now multiply equation :eq:`bombone` by `\Delta t`, divide it by `\rho_i`,
and rearrange:

.. math::
   :label: bombtwo

    \left(-R - \nu w_k^n \uppair{1-\lambda/2}{\lambda/2}\right) E_{k-1}^{n+1}\\ +
    \left(1 + 2 R + \nu w_k^n (1-\lambda) \uppair{+1}{-1}\right) E_k^{n+1}\\ +
    \left(-R + \nu w_k^n \uppair{\lambda/2}{1-\lambda/2} \right) E_{k+1}^{n+1}\\
    = E_k^n + \Delta t \rho_i^{-1}\Phi_k^n

Here `\uppair{a}{b} = a` when `w_k^n\ge 0` and `\uppair{a}{b} = b` when `w_k^n < 0`.

Equation :eq:`bombtwo` has coefficients which are scaled to have no units.  It is
ready to be put in the system managed by ``enthSystemCtx``.

One way of stating the stability of first-order upwinding is to say it satisfies
a  "maximum principle" :cite:`MortonMayers`.  An example of a maximum principle
for this kind of finite difference scheme is that if `U_{k-1}^n,U_k^n,U_{k+1}^n`
are adjacent gridded values of some abstract quantity at time step `t_n`, and
if the next value satisfies the scheme

.. math::
   :label: abstractexplicit

   U_k^{n+1} = C_{-1} U_{k-1}^n + C_0 U_k^n + C_{+1} U_{k+1}^n

for *nonnegative* coefficients `C_i` summing to one, `C_{-1} + C_0 + C_{+1} = 1`, then it
follows by the triangle inequality that

.. math::

   \min\{|U_{k-1}^n|, |U_k^n|, |U_{k+1}^n|\}
   \le |U_k^{n+1}| \le \max\{|U_{k-1}^n|, |U_k^n|, |U_{k+1}^n|\}.

Thus a "wiggle" cannot appear in `\{U_k^{n+1}\}` if previous values `\{U_k^n\}` were
smoother. The proof below shows the corresponding "wiggle-free" property for scheme :eq:`bombtwo`.

However, the pure implicit centered difference scheme (`\lambda=1`), namely

.. math::
   :label: centered

    &\left(-R - \nu w_k^n/2\right) E_{k-1}^{n+1} + \left(1 + 2 R\right) E_k^{n+1} \\
    &+ \left(-R + \nu w_k^n/2\right) E_{k+1}^{n+1}
    = E_k^n + \Delta t \rho_i^{-1}\Phi_k^n

is *less stable* than implicit first-order upwinding. It is less stable in the same sense
that Crank-Nicolson is a less stable scheme than backwards Euler for the simplest heat
equation `u_t = u_{xx}` :cite:`MortonMayers`. In fact, although oscillatory modes cannot
grow exponentially under equation :eq:`centered`, those modes *can* appear when none are
present already, even in the homogeneous case `\Phi_k^n=0`.

Stability properties of the BOMBPROOF scheme
--------------------------------------------

We want to be precise about the phrase "unconditionally stable" for BOMBPROOF.
To do so we consider somewhat simplified cases which are amenable to analysis, and
we prove two stability properties.  These stability properties identify the
precise advantages of BOMBPROOF.

Theorem (stating the stability properties).
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assume, for the precise but limited assertion of this theorem, that the surface
temperature `T_s` and the geothermal flux `G` are constant in time. Assume also that the
entire source function `\Phi` is identically zero (but see comments below). Fix an
equally-spaced vertical grid `z_0=0 < z_1 < \dots < z_N=H`, so that the upper grid point
coincides with the surface of the ice. With these assumptions, if

.. math::
   :label: lambdachoice

   \lambda = \min
   \left\{
   1, \min_{k=0,\dots,N} \left\{ \frac{2 k_i}{|w_k^n| \rho_i c_i \Delta z} \right\}
   \right\},

reset at each time step `n`, then scheme :eq:`bombone`, :eq:`bombtwo` is
unconditionally-stable in the following two senses:

1. A maximum principle applies without further assumptions.
 
2. Suppose we freeze the coefficients of the problem to have constant values in time and
   space. (Concretely, we assume that `\lambda` is chosen independently of the time step
   `n`, and that `\Delta t` is the same for each time step. We assume constant vertical
   velocity `w_k^n=w_0`. We also consider a spatially-periodic or unbounded version of our
   problem, with no boundary conditions.) Then a von Neumann analysis of the constant
   coefficient problem yields a growth factor less than one for all modes on the grid.

Remarks
^^^^^^^

The phrases *maximum principle* and *von Neumann analysis* will be precisely illustrated
in the following proof. Both approaches are in :cite:`MortonMayers`. There is additional
information on the von Neumann analysis of implicit finite difference methods for
advection in :cite:`Strikwerda`.

These statements also apply in case `k_i=0`, in which case :eq:`lambdachoice` implies
`\lambda=0`, and the method reduces to implicit first-order upwinding. (Implicit
first-order upwinding has properties 1 and 2 :cite:`Strikwerda`.) The case `k_i=0` is
relevant because it applies to the least-transport model of temperate ice in which there
is zero enthalpy conduction. (One reasonable model for temperate ice is to assume no
transport of the liquid fraction, whether diffusive transport or otherwise, and to ignore
conduction along the temperature gradient, because the gradient is only from
pressure-melting temperature differences.)

Proof of 1
^^^^^^^^^^

In the case considered for the maximum principle, with `\Phi_k^n=0`,
we can rewrite :eq:`bombtwo` as

.. math::
   :label: formax

    &\left(1 + 2 R + \nu w_k^n (1-\lambda) \uppair{+1}{-1}\right) E_k^{n+1} \\
    &= E_k^n + \left(R + \nu w_k^n \uppair{1-\lambda/2}{\lambda/2}\right) E_{k-1}^{n+1}
                    + \left(R - \nu w_k^n \uppair{\lambda/2}{1-\lambda/2}\right) E_{k+1}^{n+1}.

We claim that with choice :eq:`lambdachoice` for `0 \le \lambda \le 1`, all
coefficients in :eq:`formax` are nonnegative.  At one extreme, in
the upwinding case (`\lambda=0`), all the coefficients are nonnegative.  Otherwise, note that
`\nu w_k^n (1-\lambda) \uppair{+1}{-1}` is nonnegative for any valid value
of `\lambda` and for any value of `w_k^n`, noting the meaning of the `\uppair{+1}{-1}`
symbol.  Thus the coefficient on the left is always nonnegative.  The coefficient of 
`E_{k-1}^{n+1}` is clearly nonnegative for any valid value of `\lambda` if `w_k^n \ge 0`.
The coefficient of `E_{k+1}^{n+1}` is clearly nonnegative for any valid value of `\lambda` if
`w_k^n \le 0`.

Therefore the only concerns are for the coefficient of `E_{k-1}^{n+1}` when `w_k^n\le 0` and the
coefficient of `E_{k+1}^{n+1}` when `w_k^n\ge 0`.  But if `\lambda` is smaller than
`2k_i/(|w_k^n| \rho_i c_i \Delta z)` then

.. math::

   R - \nu |w_k^n| (\lambda/2) = \frac{k_i \Delta t}{\rho_i c_i \Delta z^2} - \frac{\Delta t |w_k^n|}{\Delta z} \frac{\lambda}{2} \ge \frac{k_i \Delta t}{\rho_i c_i \Delta z^2}
    - \frac{\Delta t |w_k^n|}{\Delta z} \frac{k_i}{|w_k^n| \rho_i c_i \Delta z} = 0.

Thus all the coefficients in :eq:`formax` are nonnegative.  On the other hand, in equation :eq:`formax`, all coefficients on the right side sum to

.. math::

   1+2R+\nu w_k^n \uppair{1-\lambda}{-1+\lambda} = 1+2R+\nu w_k^n (1-\lambda) \uppair{+1}{-1},

which is exactly the coefficient on the left side of :eq:`formax`.  It follows that

.. math::

   E_k^{n+1} = a_k E_k^n + b_k E_{k-1}^{n+1} + c_k E_{k+1}^{n+1}

where `a_k,b_k,c_k` are positive and `a_k+b_k+c_k=1`. Thus a maximum principle applies
:cite:`MortonMayers`. **END OF PROOF OF 1.**

Proof of 2
^^^^^^^^^^

As a von Neumann analysis is much more restrictive than the analysis above, we will be brief.
Let's assume the velocity is downward, `w_0<0`; the other case is similar.  Equation
:eq:`bombtwo` becomes

.. math::
   :label: prevon

    &\left(-R - \nu w_0 (\lambda/2)\right) E_{k-1}^{n+1}
    + \left(1 + 2 R - \nu w_0 (1-\lambda)\right) E_k^{n+1} \\
    & + \left(-R + \nu w_0 (1-\lambda/2) \right) E_{k+1}^{n+1}  = E_k^n.

The heart of the von Neumann analysis is the substitution of a growing or decaying
(in time index `n`) oscillatory mode on the grid of spatial wave number `\mu`:

.. math::

   E_k^n = \sigma^n e^{i\mu\,(k\Delta z)}.

Here `k\Delta z = z_k` is a grid point. Such a mode is a solution to :eq:`prevon` if and
only if

.. math::

    \sigma\Big[  &(-R - \nu w_0(\lambda/2)) e^{-i\mu\Delta z}
    + (1 + 2 R - \nu w_0 (1-\lambda)) \\
    &+ (-R  + \nu w_0 (\lambda/2)) e^{+i\mu\Delta z}
    + \nu w_0 (1-\lambda) e^{+i\mu\Delta z} \Big] = 1.

This equation reduces by standard manipulations to

.. math::

   \sigma = \frac{1}{1 + \left(4 R - 2 \nu w_0 (1-\lambda)\right)\cos^2(\mu \Delta z/2)
   + i\,\nu w_0 (1-\lambda/2)\sin(\mu\Delta z)}.

Note `4 R - 2 \nu w_0 (1-\lambda) \ge 0` without restrictions on
numerical parameters `\Delta t`, `\Delta z`, because `w_0<0` in the
case under consideration.  Therefore

.. math::

   |\sigma|^2 = \frac{1}{\left[1 + \left(4 R - 2 \nu w_0 (1-\lambda)\right)
   \cos^2(\mu \Delta z/2)\right]^2
   + \left[\nu w_0 (1-\lambda/2)\sin(\mu\Delta z)\right]^2}.

This positive number is less than one, so `|\sigma| < 1`.  It follows that all
modes decay exponentially.
**END OF PROOF OF 2.**

Remark about our von Neumann stability analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The constant `\lambda` is carefully chosen in :eq:`lambdachoice` so that the maximum
principle 1 applies. On the other hand, both the implicit first-order upwind and the
implicit centered difference formulas have unconditional stability in the von Neumann
sense. The proof of case 2 above is thus a formality, merely showing that a convex
combination of unconditionally stable (von Neumann sense) schemes is still unconditionally
stable in the same sense.

Convergence: a consequence of the maximum principle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we define the pointwise numerical error `e_k^n = E_k^n - E(t_n,x_i,y_j,z_k)`,
where `E(\dots)` is the unknown exact solution (exact enthalpy field) :cite:`MortonMayers`,
then :eq:`formax` implies an equality of the form

.. math::

   A e_k^{n+1} = e_k^n + B_- e_{k-1}^{n+1} + B_+ e_{k+1}^{n+1} + \Delta t\, \tau_k^n

where `\tau_k^n` is the truncation error of the scheme and `A,B_\pm` are nonnegative
coefficients, which need no detail for now other than to note that `1 + B_- + B_+ = A`.
Letting `{\bar e}^n = \max_k |e_k^n|` we have, because of the positivity of coefficients,

.. math::
   :label: prebound

    A |e_k^{n+1}| \le {\bar e}^n + \left(B_- + B_+\right){\bar e}^{n+1} + \Delta t\,\bar\tau^n

for all `k`, where `\bar\tau^n = \max_k |\tau_k^n|`.  Now let `k` be the index for
which `|e_k^{n+1}| = {\bar e}^{n+1}`.  For that `k` we can replace `|e_k^{n+1}|` in
equation :eq:`prebound` with `{\bar e}^{n+1}`.  Subtracting the same quantity from
each side of the resulting inequality gives

.. math::

   {\bar e}^{n+1} \le {\bar e}^n + \Delta t\,\bar\tau^n,

It follows that `\bar e^n \le C \Delta t`, for some finite `C`, if `\bar e^0 = 0`
:cite:`MortonMayers`.  Thus a maximum principle for BOMBPROOF implies convergence
in the standard way :cite:`MortonMayers`.  This convergence proof has the same
assumptions as case 1 in the theorem, and thus it only *suggests* convergence
in any broad range of glaciologically-interesting cases.


Remark on nonzero source term
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now recall we assumed in Theorem 1 that the entire "source" `\Phi_k^n` was identically zero.
Of course this is not realistic.  What we understand is provable, however, is that if a
numerical scheme for a linear advection/conduction equation

.. math::

   u_t + A u_x = B u_{xx}

is stable in the general sense of numerical schemes for partial differential equations
(e.g. as defined in subsection 5.5 of :cite:`MortonMayers`) then the same scheme is stable
in the same general sense when applied to the equation with (linear) lower order terms:

.. math::

   u_t + A u_x = B u_{xx} + C u + D.

A precise statement of this general fact is hard to find in the literature, to put it
mildly, but theorem 2.2.3 of :cite:`Strikwerda` is one interesting case (`B=0` and `D=0`).
But even the form we state with linear term (`C u + D`) is not adequate to the job because
of the strongly-nonlinear dependence of `\Phi` on the temperature `T` :cite:`BBL`.

Nonetheless the maximum principle is a highly-desirable form of stability because we can
exclude "wiggles" from the finite difference approximations of the conductive and
advective terms, even if the complete physics, with strain heating in particular, is not
yet shown to be non-explosive. Because the complete physics includes the appearance of the
famous "spokes" of EISMINT II, for example, a maximum principle cannot apply too
literally. Indeed there is an underlying fluid instability :cite:`BBL`, one that means the
solution of the continuum equations can include growing "wiggles" which are fluid features
(though not at the grid-based spatial frequency of the usual numerical wiggles). Recall
that, because we use first-order upwinding on the horizontal advection terms, we can
expect maximum principle-type stability behavior of the whole three-dimensional scheme.

.. _sec-melt:

Temperate basal boundary condition, and computing the basal melt rate
---------------------------------------------------------------------

At the bottom of grounded ice, a certain amount of heat comes out of the earth and either
enters the ice through conduction or melts the base of the ice. On the one hand, see the
documentation for ``BedThermalUnit`` for the model of how much comes out of the
earth. On the other hand, :cite:`AschwandenBuelerKhroulevBlatter` includes a careful
analysis of the subglacial layer equation and the corresponding boundary conditions and
basal melt rate calculation, and the reader should consult that reference.

Regarding the floating case
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The shelf base temperature `T_{sb}` and the melt rate `M` are supplied by ``OceanModel``.
Note that we make the possibly-peculiar physical choice that the shelf base temperature is
used as the temperature at the *top of the bedrock*, which is actually the bottom of the
ocean. This choice means that there should be no abrupt changes in top-of-bedrock heat
flux as the grounding line moves. This choice also means that the conservation of energy
code does not need to know about the bedrock topography or the elevation of sea level. (In
the future ``OceanModel`` could have a ``subshelf_bed_temperature()`` method.)
