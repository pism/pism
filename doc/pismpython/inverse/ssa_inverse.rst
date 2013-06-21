.. _SSAInverse:

SSA Inverse Problems
====================


The following data are needed for an SSA inverse problem in PISM

  * The choice of physical design variable (:math:`\tau_c` or :math:`B`).
  * The choice of :ref:`design variable parameterization <DesignParam>`.
  * An initial best estimate for the physical design variable.
  * Choices of state and design functionals :math:`J_S` and :math:`J_D`
    as described in 
    the :ref:`inverse problems refresher <inverse-background>`.
  * The choice of inverse strategy and algorithm.
  * Regularization constants (i.e. :math:`\delta` 
    for :ref:`gradient algorithms <InverseGradient>` or :math:`\eta`
    for :ref:`Tikhonov algorithms <Tikhonov>`).

The :ref:`pismi` command documentation describes how to specify these
to run an inversion.  All of these inputs have been described previously,
except for the specific choices of functionals and
algorithms, which we turn to now.

.. _statefunc:

State Functionals
-----------------

State functionals are specified with the :cfg:`-inv_state_func` flag.

.. _meansquare:

*  **Weighted Mean Square** (:cfg:`meansquare`) [Default]

   A common choice of misfit functional has the form
 
   .. math::
     J_{\misfit}(\zeta) = \int_{\Omega} \left|\calF(\zeta)-\vU_\obs\right|^2,
  
   where :math:`\vU_\obs` are the observed SSA velocities.  This 
   corresponds to a state functional of the form
 
   .. math::
     J_S(Y) = \int_{\Omega} |Y|^2.

   For more flexibility, one might want to consider a weighted version of 
   this,
 
   .. math::
     J_S(Y) = \int_{\Omega} w |Y|^2.

   where :math:`w\ge 0` is a provided weight function.  The discrete analog of 
   this functional used in PISM is
 
   .. math::
      J_S(\tY) = c_N \sum_{i,j} \tw_{i,j} |\tY_{i,j}|^2

   where the sum is taken over all grid points :math:`i,j` and where:
     * :math:`c_N` is a normalization constant, and
     * :math:`\tw` is a weight vector provided in :ncvar:`vel_misfit_weight`.

   The normalization constant :math:`c_N` is selected so that if
   :math:`|\tY_{i,j}|=`\ :cfg:`inv_ssa_velocity_scale` 
   at all grid points, then :math:`J_S(\tY)=1`.
    
* **Log-Ratio** (:cfg:`log_ratio`)
  This is a functional similar to one appearing in :cite:`Morlighemetal2010`:
  
  .. math::
    J_S(\tY) = c_N \sum_i  \tw_{i,j} \left[
      \log\left( 
            \frac{|\tY_i+U_i|^2+\epsilon^2}{|U_{i}|^2+\epsilon^2}
         \right)
    \right]^2

  where

  * :math:`\tU_{(i,j)}` are the observed SSA velocities,
  * :math:`c_N` is a normalization constant, and
  * :math:`\tw` is a weight vector provided in :ncvar:`vel_misfit_weight`.
  
  The normalization constant :math:`c_N` is selected so that if
  :math:`|\tY_{i,j}|=s|\tU_{i,j}|` at all grid points, 
  then :math:`J_S(\tY)=1`, where :math:`s=`\ :cfg:`log_ratio_scale`. 

* **Log-Relative** (:cfg:`log_relative`)

 This is an experimental functional of the form

 .. math::
   J_S(\tY) = c_N \sum_{i,j} 
          \log\left( 
              1 + \tw_{i,j}\frac{|\tY_{i,j}|^2}{|\tU_{i,j}|^2+\epsilon^2}
               \right)

 where 
   * :math:`\tU_{(i,j)}` are the observed SSA velocities,
   * :math:`c_N` is a normalization constant, and
   * :math:`\tw` is a weight vector provided in :ncvar:`vel_misfit_weight`.

   The normalization constant :math:`c_N` is selected so that if
   :math:`|\tY_{i,j}|=`\ :cfg:`inv_ssa_velocity_scale` 
   at all grid points, then :math:`J_S(\tY)=1`.



Note that all these functionals supports grid points without SSA velocity
observations by setting the weight function :math:`\tw=0` at such points.

.. _designfunc:

Design Functionals
------------------

Design functionals are specified with the :cfg:`-inv_design_func` 
flag.

* **Sobolev** :math:`H^1` (:cfg:`sobolevH1`) [Default]

  The primary design functional has the form
  
  .. math::
    J_D(Z) = \frac{1}{|\Omega|} \int_\Omega \ell^2 c_{H^1} |\nabla Z|^2 + c_{L^2} Z^2

  where
  
  * :math:`|\Omega|` is the area of the rectangular grid domain,
  * :math:`\ell=` :cfg:`inv_ssa_length_scale`,
  * :math:`c_{H^1}=` :cfg:`inv_design_cH1`, and
  * :math:`c_{L^2}=` :cfg:`inv_design_cL2`.
  
  Integration is done with 
  numerical quadrature of finite element functions.
  
  Typical values for :math:`c_{H^1}` and :math:`c_{L^2}` range between
  0 and 1, and can be specified with the option flags
  :cfg:`-inv_design_cH1` and :cfg:`-inv_design_cL2`. 
  Setting either (but not both!) of these equal to zero is acceptable.  Note 
  that :math:`\zeta` is scaled to have typical values of 1, and hence typical
  values of :math:`J_D` are expected to be on the order of 1 as well.
  
  The purpose of the design functional is to determine distances
  between values of :math:`\zeta` and the original best estimate
  :math:`\zeta_0` via
  
  .. math::
    \zeta \mapsto J_D(\zeta-\zeta_0).
    
  Setting :math:`c_{H^1}` to a non-zero value penalizes wiggles and sharp   
  derivatives in the difference :math:`\zeta-\zeta_0`.  If the initial 
  estimate :math:`\zeta_0` is smooth, then :math:`\zeta` recovered by
  inversion will tend to be smooth when :math:`c_{H^1}\neq 0`.  
  Conversely, if the initial estimate :math:`\zeta_0` contains sharp features,
  :math:`\zeta` recovered by inversion will tend to keep those same sharp 
  features because
  
  .. math::
    \zeta = \zeta_0 + (\zeta-\zeta_0)
    
  and hence :math:`\zeta` is a smooth perturbation of the original estimate.
  
  This behavior is generally desirable, but can be problematic at the boundary   
  between grounded ice and floating ice when :math:`\tau_c` is 
  the physical design variable.  At this boundary there will be a 
  jump in :math:`\zeta_0`, and a jump in the inverted value of :math:`\zeta`,
  but there is not a good reason to try to enforce that there will be 
  exactly the same jump.  To avoid such artifacts, use the flag 
  :cfg:`-inv_ssa_grounded_ice_tauc`.  When this flag is set, the
  integral omits any floating or ice-free regions and therefore does not
  artificially penalize jumps in :math:`\zeta` at these boundaries. 
  (Specifically, finite elements are omitted from the integral if any
  of the nodes is ice-free or floating).

  Setting :math:`c_{H^1}=0` results in an :math:`L^2`-type functional.

* **Pseudo Total Variation** (:cfg:`tv`)

  The pseudo total variation functional has the form

  .. math::
    J_D(Z) = \frac{(\ell)^q}{|\Omega|} 
    \int_\Omega (\epsilon^2+|\nabla Z|^2)^{q/2}

  where 

    * :math:`|\Omega|` is the rectangular grid area,
    * :math:`\ell=` :cfg:`inv_ssa_length_scale`,
    * :math:`q=`\ :cfg:`inv_ssa_tv_exponent`,
    * :math:`\epsilon` is either specified directly with
      :cfg:`-inv_ssa_tv_epsilon` or has a default value of
      :math:`1/`\ :cfg:`Schoof_regularizing_length`.

  With these parameters, assuming that :math:`Z` is dimensionless, the
  functional is dimensionless.
    
  Strictly speaking, the total-variational functional corresponds to the case
  :math:`q=1` and :math:`\epsilon=0`.  
  Such functionals have the nice property that they do not
  penalize jumps across curves, but do penalize spikes and similar noisy
  singularities.  But the case :math:`q=1`, :math:`\epsilon=0` also causes
  numerical difficulties due to its lack of differentiability, and either of
  these parameters can be adjusted to help with this.  Note that if
  :math:`q=2` and :math:`\epsilon=0`, this is exactly the same functional
  as the Sobolev :math:`H^1` functional with :math:`c_{H^1}=1` 
  and :math:`c_{L^2}=0`.

Algorithm Selection
-------------------

.. _InvGradAlg:

Iterative Gradient Algorithms
'''''''''''''''''''''''''''''

PISM uses the ``siple`` python library :cite:`siple-web-page` to implement 
gradient algorithms.  All these algorithms approximately minimize
a sum-of-squares misfit functional

.. math::
  J_{\misfit}(\zeta) = J_S(\calF(\zeta)-\vU_\obs)

in an iterative fashion, terminating at the first iteration where
the misfit descends below a specified value.  The only 
sum-of-squares state functional currently supported by PISM is the 
:ref:`weighted mean square <meansquare>` functional,

.. math::
  J_S(\tY) = c_N \sum_{i,j} \tw_{i,j} |\tY_{i,j}|^2,

where the normalization constant :math:`c_N` is chosen so :math:`J_S=1` if
:math:`|\tY|=Y_\scale` everywhere, where :math:`Y_\scale=` 
:cfg:`inv_ssa_velocity_scale`.  The
functional therefore effectively has units of 
:math:`Y_\scale^2`.

.. _InvGradStop:

The stopping criterion is provided by a parameter 
:math:`\delta=` :cfg:`-inv_target_misfit` in 
units of :math:`m/a`, and iterations are stopped when

.. math::
  J_{\misfit}(\zeta) < \left(\frac{\delta}{Y_\scale}\right)^2.

There are three choices for the iterative algorithm for approximately 
minimizing :math:`J_\misfit`, which are specified by the command-line flag
:cfg:`-inv_ssa_method`.

*  **Steepest Descent** (:cfg:`sd`)

  This is a standard, robust, choice in the glaciology literature.
  It is also the slowest and least powerful of the minimization methods,
  and can fail to reduce the functional to the desired misfit.

*  **Nonlinear Conjugate Gradients** (:cfg:`nlcg`)

  This is a variation of the steepest descent method with superior
  speed.

* **Incomplete Gauss-Newton** (:cfg:`ign`)

  An iterative algorithm that solves a model linear inverse problem at each 
  minimization iteration using a Morozov discrepany principle stopping
  criterion for the model problems :cite:`siple-web-page`
  :cite:`Habermannetal2012`.  In many cases it is the fastest of 
  the three methods, but it can also sometimes generate solutions with more
  artifacts.

.. _TikhonovAlg:

Tikhonov Algorithms
'''''''''''''''''''
Tikhonov algorithms exactly minimize functionals of the form

.. math::
  J\Tik(\zeta) = \eta J_\misfit(\zeta) + J_D(\zeta-\zeta_0)
  
where :math:`\eta` is specified using :cfg:`-tikhonov_penalty`. Any
of the :ref:`misfit <statefunc>` and :ref:`design <designfunc>` functionals
described previously can be used.  There are choices to be made in the
algorithm used for minimizing the functional, and PISM relies on the
TAO optimization library :cite:`tao-user-ref` for much of this.  The
:cfg:`-inv_ssa_method` is used to indicate Tikhonov regularization with
a specified minimization approach.

* **TAO Limited Memory Variable Metric** (:cfg:`tikhonov_lmvm`)

  A large-scale unconstrained optimization algorithm requiring only
  function and gradient evaluations.  Hessians are approximated using
  the so-called BFGS update formula.
  
* **TAO Nonlinear Conjugate Gradient** (:cfg:`tikhonov_cg`)

  The nonlinear conjugate gradient method is used to exactly minimize the
  Tikhonov functional.  In general this is a slower algorithm than
  :cfg:`tikhonov_lmvm`.

The following algorithms are also available, but are still works in progress.

* **TAO Linearly Constrained Lagrangian** (:cfg:`tikhonov_lcl`)

* **TAO Bound Constraint Limited Memory Variable Metric** (:cfg:`tikhonov_blmvm`)

  A variation of :cfg:`tikhonov_lmvm` that enforces the constraint :math:`\zeta\ge 0`.  It is intended to be used only with 
  :cfg:`-inv_design_param ident`.
  
* **Gauss Newton** (:cfg:`tikhonov_gn`)

.. _TikConverge:

Tikhonov Convergence
''''''''''''''''''''

TAO minimization routines detect convergence based on parameters set by
flags :cfg:`-tao_fatol`, :cfg:`-tao-frtol` and some others.  See the TAO
User's Manual :cite:`tao-user-ref` for details.  In addition to these stopping criteria, PISM adds an additional convergence check.

The Tikhonov functional has the form

.. math::
  J_\Tik(\zeta) = \eta J_\misfit(\zeta) + J_D(\zeta-\zeta_0)

and at a minimizer :math:`\zeta_\reg`, 
:math:`\nabla J_\Tik(\zeta_\reg)=0`.  Hence

.. math::
  \nabla J_D(\zeta_\reg-\zeta_0) = -\eta \nabla J_\misfit(\zeta_\reg).

So convergence occurs when :math:`\nabla J_D(\zeta_\reg-\zeta_0)` and
:math:`\nabla J_\misfit(\zeta_\reg)` point in opposite directions (and
have the correct relative lengths determined by :math:`\eta`).  This leads to
the condition

.. math::
  |\nabla J_\Tik(\zeta)| < \epsilon \max(|\nabla J_D(\zeta-\zeta_0)|, |\nabla J_\misfit(\zeta)|)

where :math:`\epsilon` is specified by :cfg:`-tikhonov_rtol`.



