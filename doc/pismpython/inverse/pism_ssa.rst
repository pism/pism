.. _PISM_SSA:

PISM's SSA Model
=================

Continuous Model
------------------------

The PISM model for SSA velocities :math:`\vU=(U_1,U_2)` is

.. math::
  -\sum_{i=1}^2 \partial_i \left( 2[\nu H+\epsilon_{\mathrm{SSA}}]
  \left[D_{ij}\vU + \Div \vU \;\delta_{ij}\right]\right) + \tau_{b,j} = f_j
  :label: SSA

for :math:`j=1,2`, and where

 * :math:`\nu` is the effective viscosity,
 * :math:`\mathbf{\tau}_{b}` is the basal shear stress,
 * :math:`H` is ice thickness :ncvar:`thk`,
 * :math:`\epsilon_{\mathrm{SSA}}` is a regularization parameter :cfg:`epsilon_ssa`
 * :math:`\mathbf{f}` is the driving stress :math:`-\rho g H \nabla h`, where

   * :math:`\rho` is ice density :cfg:`ice_density`
   * :math:`g` is gravitational acceleration, :cfg:`standard_gravity`
   * :math:`h` is ice surface elevation :ncvar:`usurf`

and where 

.. math::
  D_{ij} \vU =  \frac{1}{2}\left[ \frac{\partial U_i}{\partial x_j} + \frac{\partial U_j}{\partial x_i}\right].

The viscosity :math:`\nu` and basal shear stress :math:`\tau_{b}`
are themselves functions of :math:`\vU`.  Viscosity is given by
any of a few flow laws, one choice being thermocoupled Glen ice

.. math::
  \nu(\vU) = \frac{B}{2}\left( \epsilon_{\mathrm{Schoof}}^2+\frac{1}{2} |\vD \vU|^2 + \frac{1}{2} (\Div \vU)^2 \right)^{(p-2)/2}

where 

  * :math:`p=1+(1/n)` where :math:`n` is the Glen flow exponent :cfg:`Glen_exponent`.
  * :math:`B` is vertically averaged ice hardness,
  * :math:`\epsilon_{\mathrm{Schoof}}` is a regularizing parameter
    (:cfg:`Schoof_regularizing_velocity` divided by 
    :cfg:`Schoof_regularizing_length`).

All of the flow laws contain a linear factor :math:`B`, which is one of
the potential inversion design variables.

The basal shear stress is determined by a pseudo-plastic till law

.. math::

  \mathbf{\tau}_b(\vU) = \tau_c\; U_{\mathrm{th}}^q (\epsilon_{\mathrm{reg}}^2+ |\vU|^2)^{\frac{q-1}{2}}\vU

where

  * :math:`\tau_c` is effective yield stress, a potential inversion design variable,
  * :math:`U_{\mathrm{th}}` is a threshold velocity (:cfg:`pseudo_u_threshold`) 
  * :math:`\epsilon_{\mathrm{reg}}` is a regularization parameter (:cfg:`plastic_regularization`)
  * :math:`0\le q \le 1`   
    controls till plasticity; 1 is plastic and 0 is linearly-viscous
    (:cfg:`pseudo_plastic_q`).

The surface elevation :math:`h` (:ncvar:`usurf`) is not part of 
PISM's model state, but is computed from a combination of 
bedrock elevation, :ncvar:`topg`,
ice thickness :math:`H=`\ :ncvar:`thk`, and sea level.
For grounded ice, :math:`h` is the sum of thickness and bedrock elevation.
Otherwise :math:`h` is the surface elevation of floating ice with the given thickness. See :ref:`SSADiscrete` for a discussion of how the distinction
between grounded and floating ice is made.

Vertically averaged ice hardness :math:`B` is also not part of the model
state.  It is computed from the PISM 3-d enthalpy (:ncvar:`enthalpy`) 
variable; the current ice flow law has the responsibility of 
converting enthalpy and pressure at depth into hardness, which is 
then averaged over an ice column to determine :math:`B`.

Boundary Conditions
-------------------

PISM solves the SSA on a rectangular domain. For regional models (:cfg:`-regional`), the domain is non-periodic, whereas the domain
is otherwise periodic on both pairs of edges.  Ice need not be present
over the entire domain, though the SSA is applied in ice-free 
regions as discussed in :ref:`SSADiscrete`.

Dirichlet boundary conditions (i.e. locations where :math:`\vU` is known)
can be turned on with the :cfg:`-ssa_dirichlet_bc` flag, in which case
the known velocities are taken from the NC variable
:ncvar:`vel_ssa_bc`.  For regional models, the Dirichlet locations are specified
indirectly via the NC mask variable :ncvar:`no_model_mask`, otherwise the NC mask variable :ncvar:`bc_mask` determines these locations.

PISM supports a calving front boundary condition :cite:`AlbrechtLevermann2012` 
that modifies the stress balance at the ice/ocean interface (config variable 
:cfg:`calving_front_stress_boundary_condition`). This boundary condition is 
**not supported**, however, by PISM's SSA inversion algorithms.


.. _SSADiscrete:

Discretization Considerations
-----------------------------

PISM supports two discretization schemes for solving the SSA: 
finite-differences 
(:cfg:`-ssa_method fd`) and finite-elements (:cfg:`-ssa_method fem`).  The
finite difference version contains support for the calving front
boundary condition which is not available in the finite element version.
On the other hand, the finite element version uses PETSc's ``SNES``
Newton-method routines for solving the nonlinear problem, which leads to a 
robust convergence criterion independent of the number of processors.  
SSA inversion in PISM is based **only** on the finite-element implementation.

PISM treats the SSA as if it applies to the entire grid domain, even in 
ice-free locations.  Each grid point can be either icy or ice-free,
and either grounded or ocean, for a total of four states.  A point
is ice-free if the ice thickness :math:`H` falls below a 
small threshold :cfg:`mask_icefree_thickness_standard`.  
The distinction between
ground and ocean is made by computing what the surface elevation would be at that location for grounded ice and for floating ice; 
the maximum elevation determines the state.

In regions where :math:`H` is zero, the term 
:math:`\nu H` in equation :eq:`SSA` vanishes, and the ellipticity of this
equation is preserved only by the regularizing constant  :math:`\epsilon_{\rm SSA}`.  PISM has a second mechanism for maintaining the ellipticity
of this equation by by enforcing a minimum value for :math:`\nu H`.  
If the ice thickness falls below
a threshold :math:`H_{\rm ext}=` :cfg:`min_thickness_strength_extension_ssa`, 
then :math:`[\nu H+\epsilon_{\mathrm{SSA}}]` is replaced with 
:math:`\nu_{\mathrm{ext}} H_{\mathrm{ext}}` where 
:math:`\nu_{\rm ext}=` :cfg:`constant_nu_strength_extension_ssa`.  Note that this replacement
depends on the thickness :math:`H_{\rm ext}`, not the icy/ice-free mask 
condition.

The value of :math:`\tau_b` is also adjusted based on the ice/ice-free grounded/ocean status of a grid point.  For floating locations, the
value is set to 0, and for ice-free ground it is a 
large constant (:cfg:`beta_ice_free_bedrock`).  Consequently, :math:`\tau_b`
depends on the effective yield stress :math:`\tau_c` only for
grounded ice.
