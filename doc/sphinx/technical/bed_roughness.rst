.. include:: ../global.txt

.. only:: html

   .. in LaTeX I can just \usepackage{txfonts}...
   
   `\require{cancel} \newcommand{\fint}{\cancel{\phantom{x}}\kern{-1.1em}\int}`
   
.. _sec-bed-roughness:

Using Schoof's parameterized bed roughness technique in PISM
============================================================

.. contents::

An explanation
--------------

If the bed elevation ``topg`` is smoothed by preprocessing then we observe a reduction in
the peak values of the SIA diffusivities. From such smoothing there is (generically) also
a reduction in the peak magnitudes of horizontal velocities from both the SIA and SSA
models. The major consequence of these reductions, through the adaptive time-stepping
mechanism, is that PISM can take longer time steps and thus that it can complete model
runs in shorter time.

Large peak diffusivities coming from bed roughness are located (generically) at margin
locations where the ice is on, or has flowed onto, fjord-like bed topography. At coarser
resolutions (e.g. 20km and up), it appears that the effect of increasing bed roughness is
not as severe as at finer resolutions (e.g. 10km, 5km and finer). Of course it is true
that the shallow models we use, namely the SIA and SSA models, are missing significant
stress gradients at the same margin locations which have large bed slopes.

Here we are emphasizing the performance "hit" which the whole model experiences if some
small part of the ice sheet is on a rough bed. That part therefore is not well-modeled
anyway, compared to the majority of the ice sheet. (Switching to full Stokes or Blatter
higher order models without major spatial adaptivity would probably imply a gain in the
balanced stress components *and* a loss of the ability to model the ice sheet at high
resolution. There is a tradeoff between completeness of the continuum model and usable
resolution needed to resolve the features of the real ice sheet.)

There exists a theory which addresses exactly this situation for the SIA model, and
explains rigorously that one should use a smoothed bed in that model, but with an
associated reduction in diffusivity. This theory explains how to improve the SIA model to
handle bed roughness more correctly, because it parameterizes the effects of
"higher-order" stresses which act on the ice as it flows over bed topography. Specifically
it shows the way to a double performance boost for PISM:

- smoothed beds give longer time steps directly, and
- the parameterized effect of the local bed roughness is to further reduce the
  diffusivity, giving even longer time-steps.


Theory
------

The theory is in Christian Schoof's (2003) *The effect of basal topography on ice sheet
dynamics* :cite:`Schoofbasaltopg2003`. His mathematical technique is to expand the
Stokes equations in two levels of horizontal scales, one for the entire ice sheet (denoted
`[L]`) and one for the horizontal scale (wavelength) of bed topography (`[S]`).
The "inner" scaling assumes that the typical ice sheet thickness `[D]` is small
compared to `[S]`, while the "outer" scaling assumes that `[S]` is small compared
to `[L]`:

.. math::

   \nu = \frac{[D]}{[S]} \ll 1, \qquad \delta = \frac{[S]}{[L]} \ll 1.  

Specifically, there is an "inner" horizontal variable `x` describing the local
topography on scales comparable to `[S]` or smaller, and an "outer" horizontal
variable `X` describing the large scale bed topography and ice sheet flow on scales
larger than `[S]`.

In order to describe the Schoof scheme using PISM notation, we start by recalling the mass
continuity equation which is fundamental to any shallow ice theory:

.. math::
   
    \frac{\partial H}{\partial t} = (M - S) - \nabla\cdot \mathbf{q}. 

Within PISM this equation is handled by GeometryEvolution. Recall that `M-S` is the mass
balance added to the ice column per time. (It plays no further role here.) In the SIA case
with zero basal sliding, the horizontal mass flux is

.. math::
   
    \mathbf{q} = - D_{SIA} \nabla h, 

where `D_{SIA}\ge 0` is given next. Thus the mass continuity equation is *diffusive*. The
diffusivity `D_{SIA}` is a function of the ice geometry and the ice flow law. In the
isothermal Glen power law (power `= n`) case we recall

.. math::
   :label: eq-siadiffusivity
   
   D_{SIA} = \Gamma H^{n+2} |\nabla h|^{n-1}

where `\Gamma = 2 A (\rho g)^n / (n+2)` (c.f. details in :cite:`BLKCB`\).

Consider now the "original" bed topography `b_0(x_1,x_2)`, which we assume for the moment
is independent of time. (Time-independence is not actually critical, and such a
restriction can be removed.) We will use `x_1,x_2` to denote the horizontal model
coordinates, though they are denoted `x,y` elsewhere in these PISM docs. Suppose a
locally-smoothed bed is computed by averaging `b_0` over a rectangular region of sides
`2\lambda_1` by `2\lambda_2`, namely:

.. math::
   
    b_s(x_1,x_2) = \fint b_0(x_1+\xi_1,x_2+\xi_2)\,d\xi_1\,d\xi_2

where the slashed integral symbol is defined as

.. math::
   
    \fint F(\xi_1,\xi_2)\,d\xi_1\,d\xi_2 = \frac{1}{(2\lambda_1)(2\lambda_2)} \int_{-\lambda_1}^{\lambda_1} \int_{-\lambda_2}^{\lambda_2} F(\xi_1,\xi_2)\,d\xi_1\,d\xi_2. 

Consider also the "local bed topography"

.. math::
   
    \tilde b(x_1,x_2,\xi_1,\xi_2) = b_0(x_1+\xi_1,x_2+\xi_2) - b_s(x_1,x_2). 

As a function of the local coordinates `\xi_1,\xi_2`, the local bed topography `\tilde b`
is the amount by which the bed deviates from the "local average" `b_s(x_1,x_2)`. Generally
we will use `-\lambda_1 \le \xi_1 \le \lambda_1`, `-\lambda_2 \le \xi_2 \le \lambda_2` as
the smoothing domain, but these specific ranges are not required by the formulas above.
Note that the average of the local bed topgraphy is zero by definition:

.. math::
   
    \fint \tilde b(x_1,x_2,\xi_1,\xi_2)\,d\xi_1\,d\xi_2 = 0. 

The result of Schoof's scaling arguments (:cite:`Schoofbasaltopg2003`, equation (49)) is to
modify the diffusivity by multiplying by a factor `0 \le \theta \le 1`:

.. math::
   
    D = \theta(h(x_1,x_2),x_1,x_2) \, D_{SIA}. 

where `D_{SIA}` is defined by :eq:`eq-siadiffusivity` earlier, and

.. math::
   :label: eq-thetadefn
   
    \theta(h,x_1,x_2) = \left[ \fint \left(1 - \frac{\tilde b(x_1,x_2,\xi_1,\xi_2)}{H}\right)^{-(n+2)/n}\,d\xi_1\,d\xi_2 \right]^{-n}

Here the ice thickness and ice surface elevation are related to the *smoothed* bed
topography, so that in PISM notation

.. math::
   
    H(t,x_1,x_2) = h(t,x_1,x_2) - b_s(x_1,x_2). 

This can be treated as the *definition* of the ice thickness `H` in the above formula for
`\theta`.

The formula for `\theta` has additional terms if there is basal sliding, but we consider
only the non-sliding SIA here.

The very important fact that `0 \le \theta \le 1` is proven in appendix A of
:cite:`Schoofbasaltopg2003` by a Jensen's inequality argument. (See also the convexity
argument at the bottom of this page.)


Practical application, and Taylor approximation
-----------------------------------------------

The above formulas already reflect the recommendations Schoof gives on how to apply his
formulas (:cite:`Schoofbasaltopg2003`, subsection 4.2). The rest of this page is devoted to
how the class ``stressbalance::BedSmoother`` implements a practical version of this
theory, based on these recommendations plus some additional approximation.

The averages appearing in his scaling arguments are over an infinite domain, e.g.

.. math::
   
    f_s(x) = \lim_{R\to\infty} \frac{1}{2R} \int_{-R}^R f(\xi,x)\,d\xi. 

For practical modeling use, Schoof specifically recommends averaging over some finite
length scale which should be "considerably greater than the grid spacing, but much smaller
than the size of the ice sheet." Furthermore he recommends that, because of the typical
aspect ratio of ice sheets, "Bed topography on much larger length scales than 10 km should
then be resolved explicitly through the smoothed bed height `b_s` rather than the
correction factor `\theta`." Thus in PISM we use `\lambda_1 = \lambda_2 = 5` km as the
default (set :config:`stress_balance.sia.bed_smoother.range` to change this value).

It is, of course, possible to have bed roughness of significant magnitude at essentially
any wavelength. We make no claim that PISM results are good models of ice flow over
*arbitrary* geometry; clearly the current models cannot come close to the non-shallow
solution (Stokes) in such cases. Rather, the goal right now is to improve on the existing
shallow models, the diffusive SIA specifically, while maintaining or increasing
high-resolution performance and comprehensive model quality, which necessarily includes
many other modeled physical processes like ice thermal state, basal lubrication, and so
on. The desirable properties of the Schoof scheme are accepted not because the resulting
model is perfect, but because we gain both a physical modeling improvement *and* a
computational performance improvement from its use.

How do we actually compute expression :eq:`eq-thetadefn` quickly? Schoof has this suggestion,
which we follow: "To find `\theta` for values of [surface elevation for which `\theta` has
not already been computed], some interpolation scheme should then be used. `\theta` is
then represented at each grid point by some locally-defined interpolating function [of the
surface elevation]."

We need a "locally-defined interpolating function". As with any approximation scheme,
higher accuracy is achieved by doing "more work", which here is an increase in memory used
for storing spatially-dependent coefficients. Pade rational approximations, for example,
were considered but are excluded because of the appearance of uncontrolled poles in the
domain of approximation. The 4th degree Taylor polynomial is chosen here because it shares
the same convexity as the rational function it approximates; this is proven below.

Use of Taylor polynomial `P_4(x)` only requires the storage of three fields (below),
but it has been demonstrated to be reasonably accurate by trying beds of various
roughnesses in Matlab/Octave scripts. The derivation of the Taylor polynomial is most
easily understood by starting with an abstract rational function like the one which
appears in :eq:`eq-thetadefn`, as follows.

The fourth-order Taylor polynomial of the function `F(s)=(1-s)^{-k}` around `s=0` is

.. math::
   
    P_4(s) = 1 + k s + \frac{k(k+1)}{2} s^2 + \dots + \frac{k(k+1)(k+2)(k+3)}{4!} s^4, 

so `F(s) = P_4(s) + O(s^5)`.  Let

.. math::
   
    \omega(f,z) = \fint (1 - f(\xi) z)^{-k}\,d\xi 

where `f` is some function and `z` a scalar.  Then

.. math::

   \omega(f,z) &= \fint (1 - f(\xi) z)^{-k}\,d\xi = \fint P_4(f(\xi) z)\,d\xi + O((\max |f(\xi)|)^5 |z|^5)
   
   &\approx 1 + k z \fint f(\xi)\,d\xi + \frac{k(k+1)}{2} z^2 \fint f(\xi)^2\,d\xi + \dots + \frac{k(k+1)(k+2)(k+3)}{4!} z^5 \fint f(\xi)^4\,d\xi.


Now, `\theta` can be written

.. math::
   
    \theta(h,x_1,x_2) = \left[ \omega(\tilde b(x_1,x_2,\cdot,\cdot), H^{-1}) \right]^{-n}.  

So our strategy should be clear, to approximate `\omega(\tilde b(x_1,x2,\cdot,\cdot),
H^{-1})` by the Taylor polynomial as a function of `H^{-1}`, whose the coefficients depend
on `x_1,x_2`. We thereby get a rapidly-computable approximation for `\theta` using stored
coefficients which depend on `x_1,x_2`. In fact, let `f(\xi) = \tilde
b(x_1,x_2,\xi_1,\xi_2)` for fixed `x_1,x_2`, and let `z=H^{-1}`. Recall that the mean of
this `f(\xi)` is zero, so that the first-order term drops out in the above expansion of
`\omega`. We have the following approximation of `\theta`:

.. math::
   :label: eq-thetaapprox

    \theta(h,x_1,x_2) \approx \left[ 1 + C_2(x_1,x_2) H^{-2} + C_3(x_1,x_2) H^{-3} + C_4(x_1,x_2) H^{-4} \right]^{-n}

where

.. math::
   
    C_q(x_1,x_2) = \frac{k(k+1)\dots(k+q-1)}{q!} \fint \tilde b(x_1,x_2,\xi_1,\xi_2)^q\,d\xi_1\,d\xi_2 

for `q=2,3,4` and `k = (n+2)/n`.

Note that the coefficients `C_q` depend only on the bed topography, and not on the ice
geometry. Thus we will pre-process the original bed elevations `b_0` to compute and store
the fields `b_s,C_2,C_3,C_4`. The computation of equation :eq:`eq-thetaapprox` is fast and
easily-parallelized if these fields are pre-stored. The computation of the coefficients
`C_q` and the smoothed bed `b_s` at the pre-processing stage is more involved, especially
when done in parallel.

The parameters `\lambda_1,\lambda_2` must be set, but as noted above we use a default
value of 5 km based on Schoof's recommendation. This physical distance may be less than or
more than the grid spacing. In the case that the grid spacing is 1 km, for example, we see
that there is a large smoothing domain in terms of the number of grid points. Generally,
the ghosting width (in PETSc sense) is unbounded. We therefore move the unsmoothed
topography to processor zero and do the smoothing and the coefficient-computing there. The
class ``stressbalance::BedSmoother`` implements these details.

Convexity of `P_4`
^^^^^^^^^^^^^^^^^^

The approximation :eq:`eq-thetaapprox` given above relates to the Jensen's inequality argument
used by Schoof in Appendix A of :cite:`Schoofbasaltopg2003`. For the nonsliding case, his
argument shows that because `F(s)=(1-s)^{-k}` is convex on `[-1,1]` for `k>0`, therefore
`0 \le \theta \le 1`.

Thus it is desirable for the approximation `P_4(z)` to be convex on the same interval,
and this is true. In fact,

.. math::
   
    P_4''(s) =  k(k+1) \left[1 + (k+2) s + \frac{(k+2)(k+3)}{2} s^2\right], 

and this function turns out to be positive for all `s`. In fact we will show that the
minimum of `P_4''(s)` is positive. That minimum occurs where `P_4'''(s)=0` or
`s=s_{min}=-1/(k+3)`. It is a minimum because `P_4^{(4)}(s)` is a positive
constant. And

.. math::
   
    P_4''(s_{min}) =  \frac{k(k+1)(k+4)}{2(k+3)} > 0. 



