.. _inverse-background:

============================
Inverse Problems Background
============================

Suppose :math:`\calF` is a map from variables
:math:`\vd` to variables :math:`\vu`, i.e.

.. math:: \vu = \calF(\vd).

The corresponding forward problem is:

  Given :math:`\vd`, compute :math:`\calF(\vd)` to determine :math:`\vu`.

The inverse problem goes the other way:

  Given :math:`\vu`, find :math:`\vd` such that 
  :math:`\vu=\calF(\vd)`.

.. note::
  We will call :math:`\vd` the design variables and :math:`\vu` the 
  state variables, the point being that one controls the design to 
  determine the state.  As a concrete example, in
  PISM's inverse problems :math:`\vu` are 
  SSA velocities, and :math:`\vd` are other parameters occurring 
  the SSA equations. 

Frequently the inverse problem has undesirable features:

 #. Solutions may not exist.
 #. Multiple solutions may exists.
 #. The problem may be unstable.

A problem with any one of these properties is called ill-posed, and
the term *inverse problem* usually means an ill-posed problem.

The most important of these undesirable features is the third, which 
is sometimes called *discontinuity with respect to the data*.  It
refers to problems where it is difficult or impossible
to bound the error in :math:`\vu` from errors that may be 
present in :math:`\vd`.  

Instability
-----------

Consider the forward problem of finding a function
:math:`u` on :math:`[0,2\pi]` such that 

.. math:: 
  u'' = d

  u(0)=u(2\pi)=0

where :math:`d` is given data. The solution can be determined by means 
of a Fourier sine series

.. math::
  u = \sum_{k=1}^\infty c_k \sin(kx)

where

.. math::
  c_k = -\frac{1}{\pi k^2}\int_0^{2\pi} \sin(kx) d(x) \;dx.

The inverse problem is to determine :math:`d` from :math:`u`, which
is nothing more than taking two derivatives.  However,
taking a derivative is an unstable operation. 
Suppose :math:`u` is contaminated with small errors

.. math::
  u_{\epsilon}(x)  = u(x) + \epsilon \sin(2\pi Kx)
  
for some integer :math:`K`.  Note that :math:`\epsilon`, not :math:`K`
determines the size of the errors. The solution of the inverse problem 
for :math:`u_{\epsilon}` is 

.. math::
  d_{\epsilon}(x) = d(x) + \epsilon (2\pi K)^2 \sin(2\pi K x)
  
which has an error of order :math:`K^2` relative to :math:`\epsilon`.  
Consequently, we cannot control 
the size of the error in :math:`d_{\epsilon}` from the size of the error
in :math:`u_\epsilon`.

Regularization
--------------

The inverse problems considered in PISM suffer from instability (see, e.g., 
the demonstration in :cite:`Habermannetal2012`), and it is important to
understand techniques used to mitigate this problem.  

We assume that instead of knowing the true values :math:`\vu` of the
state variables, we only know an approximation :math:`\vu_{\mathrm{obs}}`.
Imagine that these are SSA velocities measured from observations.
There are errors present in them because

  1. they possess measurement errors, and
  2. they possess model errors (i.e. the SSA is only a simplified model 
     and nature does not conform exactly to that model).

The idea behind regularization is to **not** solve 

.. math::
  \vu_{\mathrm{obs}} = \calF(\vd)
  
exactly, since the errors in :math:`\vuobs` will lead to unacceptable
contamination of :math:`\vd` with errors.  Rather, we solve
an alternative, *regualarized* problem to find :math:`\vd_{\mathrm{reg}}` such that

.. math::
  \vu_{\mathrm{obs}} \approx \calF(\vd_{\mathrm{reg}})

and such that :math:`\vd_{\mathrm{reg}}` has been constrained on only recover
features that are justified by the error levels in :math:`\vuobs` 
and by additional information we may have about the solution.

PISM contains algorithms for two regularization strategies:

  1. :ref:`InverseGradient`
  2. :ref:`Tikhonov`

.. _InverseGradient:

Iterative Gradient Algorithms
"""""""""""""""""""""""""""""

For these algorithms, we start by imposing *sum-of-squares* type
functionals :math:`J_S` and :math:`J_D` for measuring the size of 
errors in the state and design spaces.  For example

.. math::
  J_S(\vu) = \int_{\Omega} |\vu|^2

or

.. math::
  J_S(\vu) = \int_{\Omega} |\nabla \vu|^2

would be legitimate state-space functionals for SSA problems, where :math:`\Omega` is the domain where the SSA is to be solved. 

We then introduce a functional to minimize:

.. math::
  J_\misfit(d) = J_S(\vu_{\obs}-\calF(d)).
  :label: J-misfit

If there exists a solution of :math:`\calF(d)=\vu_{\obs}`, it will be a
minimizer of :eq:`J-misfit`.  But we don't want to find this exact solution
(i.e. the exact minimizer) because it will be contaminated with errors.  

.. note::
  We call :math:`J_\misfit` the misfit functional, but we use the 
  word *misfit* for :math:`\sqrt{J_\misfit}` for sum-of-squares misfit
  functionals.  See also the :ref:`remarks below <misfit-clarify>`.

The regularized problem requires two pieces of *a-priori* data:

  1. An estimate for the size :math:`\delta` of the errors in
     :math:`\vu_{\obs}`.  That is, we  need an estimate :math:`\delta` 
     such that

       .. math::  J_S(\vu-\vu_{\obs}) < \delta^2.

  2. An initial estimate :math:`\vd_0` for the design variables that
     contains our best *a-priori* approximation of those values.

Regularization is achieved by minimizing
:math:`J_\misfit` (using any of a number of iterative techniques)  
starting from :math:`\vd_0` and by stopping at the first iteration 
where :math:`J_\misfit(d)<\delta^2`.  The technique of stopping early
based on error estimates is known as the Morozov discrepancy principle,
and is rigorously justified for linear problems :cite:`Hanke`.

The approximate minimization determines a solution :math:`\vd_{\reg}`
that is consistent with the errors in :math:`\vu_{\obs}` and
is not too far from the initial estimate :math:`\vd_0` for the design
variables. 

The design functional :math:`J_D` has not appeared explicitly in the 
discussion so far, and its role is subtle.  During the minimization,
typically steps are taken in directions based on the direction of
steepest descent of :math:`J`; these steps depend both on
:math:`J_S` and on :math:`J_D`.   When we say that the regularized
solution is not too far from the initial estimate :math:`\vd_0`,
this is measured with respect to :math:`J_D`.  Hence the choice of
:math:`J_D` is in some sense part of the *a-priori* data:
we believe that the solution is near :math:`\vd_0` as measured by :math:`J_D`.

A more complete description of iterative gradient algorithms
can be found in the ``siple`` documentation :cite:`siple-web-page`.  See 
also :cite:`Habermannetal2012`.

.. note::
  Inverse gradient algorithms are generally referred to as control theory
  algorithms in the glaciology literature, where they were introduced in
  :cite:`MacAyealtutorial`.  Traditionally, slow minimization methods (e.g. 
  steepest descent) have been used to minimize :math:`J_\misfit`, 
  which obscured the need for a careful consideration of the stopping
  criterion.  The use of the Morozov discrepancy principle stopping 
  criterion was first used for glaciology inverse problems 
  in :cite:`Maxwelletal2008`.

.. _Tikhonov:

Tikhonov Minimization
"""""""""""""""""""""""""""""

As with inverse gradient algorithms, we provide 
functionals :math:`J_S` and :math:`J_D` for measuring the size of 
errors in the state and design spaces, but now the functionals
need not be of sum-of-squares form.  For example,
the total-variation functional

.. math::
  J_D(\vd) = \int_{\Omega} |\nabla \vd|

is an acceptable functional that is not in sum-of-squares form.

.. _misfit-clarify:

Again, as with inverse gradient algorithms, we consider the 
misfit functional

.. math::
  J_\misfit(\vd) = J_S(\vu_{\obs}-\calF(\vd)),

which is minimized by an exact solution of the inverse problem.
Since :math:`J_S` need not be of sum-of-squares type,
it need not make sense to call :math:`\sqrt{J_\misfit}`.
We leave the word *misfit* undefined in general; though it 
will typically be a root of the misfit functional.
In the remainder of this section we will only speak of
"the value of the misfit functional", which is
unambiguous, and equals the square of the misfit for
sum-of-squares functionals.

For Tikhonov regularization, we provide a best-guess
initial estimate :math:`\vd_0` for the design parameters
and a penalty parameter :math:`\eta`.  We then
*exactly* minimize the Tikhonov functional

.. math::
  J_{\Tik}(\vd) = J_\misfit(\vd) + \frac{1}{\eta} J_D(\vd-\vd_0).

To explain the meaning of the parameter :math:`\eta`, suppose
that :math:`\vd_{\reg}` is a minimizer of the Tikhonov functional,
and let :math:`M=J_\misfit(\vd_\reg)` be the 
value of the misfit functional at the minimizer.
At the minimizer, the derivative of :math:`J_\Tik` vanishes and we have

.. math::
  \frac{\delta}{\delta \vd} J_D(\vd-\vd_0) = - \eta \frac{\delta}{\delta \vd}
  J_\misfit(\vd).
  :label: tik-min

Equation :eq:`tik-min` together with the equation 
:math:`J_\misfit(\vd)=M` is exactly the equation
to solve for minimizing  
:math:`\vd \mapsto J_D(\vd-\vd_0)` subject to the constraint
that :math:`J_\misfit(\vd)=M`, where :math:`-\eta` plays the
role of a Lagrange multiplier.  Moreover, for linear forward problems,
it can be shown that this is the same as minimizing 
:math:`\vd \mapsto J_D(\vd-\vd_0)` subject to the constraint  :math:`J_\misfit(\vd)\le M` :cite:`Tarantola`.

Selection of the penalty parameter :math:`\eta` therefore 
indirectly determines the acceptable value :math:`M` of the
misfit functional.
The inverse problem is regularized by seeking the design parameters 
:math:`\vd` that are closest (as measured by :math:`J_D`) 
to the initial estimate :math:`\vd_0` and such that the 
value of the misfit functional is consistent with the value 
specified indirectly by :math:`\eta`.  

.. note::
  Usually the Tikhonov functional is expressed in the form
  
  .. math::
    J_{\Tik}(\vd) = J_\misfit(\vd) + \mu J_D(\vd-\vd_0).
    :label: tikhonov_dual
  
  for some penalty parameter :math:`\mu`.  The two functionals
  are the same after identifying :math:`\mu=1/\eta`, so there is no 
  real difference between the formulations. Minimizing
  :math:`J_{\Tik}(\vd)` is equivalent to minimizing
  
  .. math::
    \vd \mapsto \eta J_\misfit(\vd) + J_D(\vd-\vd_0),
  
  which has the interpretation of minimizing the :math:`J_D` term subject
  to a constraint given by :math:`J_\misfit` (with :math:`\eta` being
  the Lagrange multiplier).  The dual interpretation of 
  equation :eq:`tikhonov_dual` is that we are minimizing :math:`J_\misfit`
  with a constraint given by the :math:`J_D` term (and :math:`\mu` being the
  Lagrange multiplier).  The constraint given by the misfit term is 
  more natural, and to emphasize this we specify the penalty weight 
  with :math:`\eta`.  For larger values of :math:`\eta`, you can expect
  that a minimization algorithm will take a greater
  number of iterations to terminate with a smaller resulting value 
  of :math:`J_\misfit`.

.. _InvAlgCompare:

Comparison of the Regularization Techniques
-------------------------------------------

Both approaches discussed above minimize the distance between the
solution design parameters :math:`\vd` and the initial estimate :math:`\vd_0`
subject to a condition that :math:`\calF(\vd)` needs to be close to
the observed state variables :math:`\vu_\obs`.

 * For gradient algorithms, :math:`\vd` is kept 
   close to :math:`\vd_0` by stoping a minimization of :math:`J_\misfit` 
   early according to the Morozov discrepancy principle.
 * For Tikhonov algorithms, :math:`\vd` is kept close to :math:`\vd_0`
   by exactly solving a constrained minimization problem encoded by
   :math:`J_\Tik`.

For gradient algorithms, the closeness is imposed heuristically 
and the solution of the inverse problem 
depends on the specific algorithm being used to approximately 
minimize :math:`J_\misfit`.  Conversely, because Tikhonov functionals 
are minimized exactly, the solution is independent of the minimization 
algorithm.

Since gradient algorithms do not exactly minimize anything, they can be
expected to terminate faster than corresponding Tikhonov algorithms.
This comes at the tradeoff of a lower quality solution: it is frequently
the case that the regularized solution :math:`\vd_\reg` is not the closest
possible to :math:`\vd_0` and contains unnecessary features
in the solution.  On the other hand, it is often the case that such 
solutions are 'good enough', which accounts for the wide-spread use of 
control theory methods in glaciology.

The condition that ":math:`\calF(\vd)` needs to be close to
the observed state variables :math:`\vu_\obs`" is specified differently
in the two approaches.  

  * For gradient algorithms, an estimate for the final value
    :math:`M` of the misfit functional 
    directly specified via :math:`\delta` according to the formula 
    :math:`M\le \delta^2`. 
        
  * For Tikhonov minimization, the final value of the misfit functional
    is specified indirectly via
    :math:`\eta`.  Once the Tikhonov functional has been minimized
    to for a fixed parameter :math:`\eta` to determine :math:`d_\reg`, 
    let :math:`M=J_\misfit(d_\reg)`.  The regularized
    solution is then the closest solution to :math:`d_0` 
    such that :math:`J_\misfit=M`. Adjusting the value of of 
    :math:`\eta` can then
    be used to adjust :math:`M`.  Increasing :math:`\eta`
    leads to a greater emphasis on the misfit in the Tikhonov functional
    and therefore a decrease in the misfit at the minimizer.

Tikhonov algorithms allow for greater flexibility in the choice of functionals
:math:`J_D` and :math:`J_S`.  For gradient methods, these must be of 
sum-of-squares type, a condition that is not needed for Tikhonov methods.

Final Remarks
-------------

There is no right answer to the solution of an ill-posed inverse problem.
Many, many approximate solutions will be consistent with the errors in the
state variables, and they will look qualitatively different from each other.
Some will be highly oscillatory, and some will not.  Some will have sharp
jumps, and some will have smooth ramps.  They will all approximately solve
the equation

.. math::
  \vu_\obs = \calF(\vd)

and none of them has true claim to being the solution of the inverse problem.

Choices made in regularizing the inverse problem affect the solution; they
pick out one from the myriad others that could have been selected.  These
choices include

  * state and design functionals
  * regularization strategy
  * initial estimates :math:`\vd_0`
  * Final values :math:`M` of the misfit functional (as determined by  choices 
    of :math:`\delta` and :math:`\eta` for gradient and Tikhonov algorithms).

There is *no preferred way* to make these choices.  
Estimates for :math:`M` are especially critical and especially hard to make.
Selecting an :math:`M` that is too large leads to loss of resolution; one that
is too small leads to spurious features in the design variables that cannot be 
trusted. One approach to estimating these parameters via L-curves can be 
found in :cite:`Habermannetal2013`.  Other methods are possible; see e.g. the 
texts :cite:`Tarantola` or :cite:`Vogel`.
