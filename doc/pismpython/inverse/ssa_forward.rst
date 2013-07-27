.. _SSAForward:

SSA Forward Problems
====================

PISM supports two forward problems for the SSA that can be passed to inversion 
algorithms.  In each case, all of the parameters appearing in the 
SSA equations described in :ref:`PISM_SSA` are held constant, except 
for one of

  1. effective yield stress :math:`\tau_c`, or
  2. vertically averaged ice hardness :math:`B`,

which are the design variables. The corresponding state variable is 
the SSA velocity :math:`\vU`.  Roughly speaking, the forward problem 
is then "given a vector of  :math:`\tau_c` (or :math:`B`) at the grid 
points, determine :math:`\vU` at the grid points".  

There are two clarifications that need to be made to this description, 
however.  For concreteness, suppose :math:`\tau_c` is the design variable,
and let :math:`\calF` be the map taking :math:`\tau_c` to :math:`\vU`,

.. math::
  \vU=\calF_{\rm SSA}(\tau_c).

This forward problem is modified by allowing for a parameterization of the 
design variable.  For example, it is frequently the case in the 
glaciology literature that :math:`\tau_c` is written as the square 
of some other parameter to 
ensure its positivity.  Other parameterization 
choices are possible, and we assume
that :math:`\tau_c` is a function of a second variable, :math:`\zeta` (:ncvar:`zeta_inv`), which is the true design variable.  So there is a function :math:`g` such that

.. math::
  \tau_c = g(\zeta)

and the forward problem is

.. math::
  \vU=\calF(\zeta)=\calF_{\rm SSA}(g(\zeta)).
  
PISM supports a number of parameterizations, described in detail in :ref:`DesignParam`.  Sometimes we will speak loosely of :math:`\tau_c`
or :math:`B` as being the design variables, but strictly speaking
it is always :math:`\zeta` that corresponds to the design variables
:math:`\vd` from :ref:`inverse-background`.  If we need to make
a distinction we will call :math:`\tau_c` (or :math:`B`) the physical
design variable and :math:`\zeta` the parameterized design variable.

The second clarification is that there are locations where :math:`\tau_c` 
is not really free to change in PISM.  At grid points where ice is floating,
or where land is ice-free, the design variable :math:`\tau_c` is ignored and
replaced in computations with different constants.  This effectively reduces
the dimension of the space of design variables. The inversion algorithms use 
a mask, :ncvar:`zeta_fixed_mask`, to indicate grid points where :math:`\zeta` 
(and hence :math:`\tau_c`) is to be held constant.  The design 
variable :math:`\zeta` is free to change only at points where 
:ncvar:`zeta_fixed_mask` :math:`=0`;  at points where
:ncvar:`zeta_fixed_mask` :math:`=1`, :math:`\zeta` is maintained at the value
it had at the start of the inversion.  By default, :cfg:`zeta_fixed_mask`
is computed automatically, though it can be provided explicitly, if desired.

.. _DesignParam:

Design Variable Parameterizations
---------------------------------

For concreteness, we suppose that :math:`\tau_c` is the
physical design variable.  PISM uses a scale parameter
:math:`\tau_{c,\scale}=`\ :cfg:`design_param_tauc_scale` to keep
:math:`\zeta` of order 1 for typical values of :math:`\tau_c`.
The following four parameterizations are supported:

  * :cfg:`ident`\ : The most straightforward transformation between :math:`\zeta` and :math:`\tau_c`:

    .. math::
       \tau_c = \tau_{c,\scale}\; \zeta.
  
    Positivity is not enforced.
  
  * :cfg:`square`: A common choice in the glaciology literature to enforce positivity:

    .. math::
      \tau_c = \tau_{c,\scale}\;\zeta^2.
      
  * :cfg:`exp`: An alternative choice for enforcing positivity:
   
    .. math::
      \tau_c = \tau_{c,\scale}\;\exp(\zeta).
  
  * :cfg:`trunc`: A kind of truncated identity map that enforces positivity.  
    For large values of :math:`\zeta`, 
    :math:`\tau_c\approx \tau_{c,0}\zeta`, and :math:`\tau_c\ra 0` as
    :math:`\zeta\ra -\infty`.  Specifically:
             
      .. math::
        \tau_c = \tau_{c,\scale}\; \frac{\zeta+\sqrt{\zeta^2+4\zeta_0^2}}{2}
    
    where :math:`\zeta_0=` :cfg:`design_param_trunc_tauc0` / :cfg:`design_param_tauc_scale`.  The parameter :cfg:`design_param_trunc_tauc0`
    is the approximate point where the linear relationship between
    :math:`\tau_c` and :math:`\zeta` begins.

The same parameterizations can be used when :math:`B` is the physical design 
variable, in which case :cfg:`hardav` replaces :cfg:`tauc` in the 
configuration variable names.

