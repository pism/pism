/*! \page bedrough Using Schoof's (2003) parameterized bed roughness technique in PISM

\section justify An explanation

If the bed elevation \c topg is smoothed by preprocessing then we observe a reduction in the peak values of the SIA diffusivities.  From such smoothing there is (generically) also a reduction in the peak magnitudes of horizontal velocities from both the SIA and SSA models.  The major consequence of these reductions, through the adaptive time-stepping mechanism, is that PISM can take longer time steps and thus that it can complete model runs in shorter time.

Large peak diffusivities coming from bed roughness are located (generically) at margin locations where the ice is on, or has flowed onto, fjord-like bed topography.  At coarser resolutions (%e.g. 20km and up), it appears that the effect of increasing bed roughness is not as severe as at finer resolutions (%e.g. 10km, 5km and finer).  Of course it is true that the shallow models we use, namely the SIA and SSA models, are missing significant stress gradients at the same margin locations which have large bed slopes.

Here we are emphasizing the performance "hit" which the whole model experiences if some small part of the ice sheet is on a rough bed.  That part therefore is not well-modeled anyway, compared to the majority of the ice sheet.  (Switching to full Stokes or Blatter higher order models without major spatial adaptivity would probably imply a gain in the balanced stress components \e and a loss of the ability to model the ice sheet at high resolution.  There is a tradeoff between completeness of the continuum model and usable resolution needed to resolve the features of the real ice sheet.)

There exists a theory which addresses exactly this situation for the SIA model, and explains rigorously that one should use a smoothed bed in that model.  But with an associated reduction in diffusivity.  This theory explains how to improve the SIA model to handle bed roughness more correctly, because it parameterizes the effects of "higher-order" stresses which act on the ice as it flows over bed topography.  Specifically it shows the way to a double performance boost for PISM:
-# smoothed beds give longer time steps directly, and
-# the parameterized effect of the local bed roughness is to further reduce the diffusivity, giving even longer time-steps.


\section schoofstheory Theory

The theory is in Christian Schoof's (2003) <i>The effect of basal topography on ice sheet dynamics</i> [\ref Schoofbasaltopg2003].  His mathematical technique is to expand the Stokes equations in two levels of horizontal scales, one for the entire ice sheet (denoted \f$[L]\f$) and one for the horizontal scale (wavelength) of bed topography (\f$[S]\f$).  The "inner" scaling assumes that the typical ice sheet thickness \f$[D]\f$ is small compared to \f$[S]\f$, while the "outer" scaling assumes that \f$[S]\f$ is small compared to \f$[L]\f$:
	\f[\nu = \frac{[D]}{[S]} \ll 1, \qquad \delta = \frac{[S]}{[L]} \ll 1.  \f]
Specifically, there is an "inner" horizontal variable \f$x\f$ describing the local topography on scales comparable to \f$[S]\f$ or smaller, and an "outer" horizontal variable \f$X\f$ describing the large scale bed topography and ice sheet flow on scales larger than \f$[S]\f$.

In order to describe the Schoof scheme using PISM notation, we start by recalling the mass continuity equation which is fundamental to any shallow ice theory:
  \f[ \frac{\partial H}{\partial t} = (M - S) - \nabla\cdot \mathbf{q}. \f]
Within PISM this equation is handled by IceModel::massContExplicitStep().  Recall that \f$M-S\f$ is the mass balance added to the ice column per time.  (It plays no further role here.)  In the SIA case with zero basal sliding, the horizontal mass flux is
  \f[ \mathbf{q} = - D_{SIA} \nabla h, \f]
where \f$D_{SIA}\ge 0\f$ is given next.  Thus the mass continuity equation is \e diffusive.  The diffusivity \f$D_{SIA}\f$ is a function of the ice geometry and the ice flow law.  In the isothermal Glen power law (power \f$= n\f$) case we recall
  \f[ D_{SIA} = \Gamma H^{n+2} |\nabla h|^{n-1} \tag{siadiffusivity} \f]
where \f$\Gamma = 2 A (\rho g)^n / (n+2)\f$ [%c.f. details in \ref BLKCB].

Consider now the "original" bed topography \f$b_0(x_1,x_2)\f$, which we assume for the moment is independent of time.  (Time-independence is not actually critical, and such a restriction can be removed.)  We will use \f$x_1,x_2\f$ to denote the horizontal model coordinates, though they are denoted \f$x,y\f$ elsewhere in these PISM docs.  Suppose a locally-smoothed bed is computed by averaging \f$b_0\f$ over a rectangular region of sides \f$2\lambda_1\f$ by \f$2\lambda_2\f$, namely:
  \f[ b_s(x_1,x_2) = \fint b_0(x_1+\xi_1,x_2+\xi_2)\,d\xi_1\,d\xi_2  \f]
where the slashed integral symbol is defined as
  \f[ \fint F(\xi_1,\xi_2)\,d\xi_1\,d\xi_2 = \frac{1}{(2\lambda_1)(2\lambda_2)} \int_{-\lambda_1}^{\lambda_1} \int_{-\lambda_2}^{\lambda_2} F(\xi_1,\xi_2)\,d\xi_1\,d\xi_2. \f]
Consider also the "local bed topography"
  \f[ \tilde b(x_1,x_2,\xi_1,\xi_2) = b_0(x_1+\xi_1,x_2+\xi_2) - b_s(x_1,x_2). \f]
As a function of the local coordinates \f$\xi_1,\xi_2\f$, the local bed topography \f$\tilde b\f$ is the amount by which the bed deviates from the "local average" \f$b_s(x_1,x_2)\f$.  Generally we will use \f$-\lambda_1 \le \xi_1 \le \lambda_1\f$, \f$-\lambda_2 \le \xi_2 \le \lambda_2\f$ as the smoothing domain, but these specific ranges are not required by the formulas above.  Note that the average of the local bed topgraphy is zero by definition:
  \f[ \fint \tilde b(x_1,x_2,\xi_1,\xi_2)\,d\xi_1\,d\xi_2 = 0. \f]

The result of Schoof's scaling arguments ([\ref Schoofbasaltopg2003], equation (49)] is to modify the diffusivity by multiplying by a factor \f$0 \le \theta \le 1\f$:
  \f[ D = \theta(h(x_1,x_2),x_1,x_2) \, D_{SIA}. \f]
where \f$D_{SIA}\f$ is defined by (siadiffusivity) earlier, and
  \f[ \theta(h,x_1,x_2) = \left[ \fint \left(1 - \frac{\tilde b(x_1,x_2,\xi_1,\xi_2)}{H}\right)^{-(n+2)/n}\,d\xi_1\,d\xi_2 \right]^{-n} \tag{thetadefn} \f]
Here the ice thickness and ice surface elevation are related to the \e smoothed bed topography, so that in PISM notation
  \f[ H(t,x_1,x_2) = h(t,x_1,x_2) - b_s(x_1,x_2). \f]
This can be treated as the \e definition of the ice thickness \f$H\f$ in the above formula for \f$\theta\f$.

The formula for \f$\theta\f$ has additional terms if there is basal sliding, but we consider only the non-sliding SIA here.

The very important fact that \f$0 \le \theta \le 1\f$ is proven in appendix A of [\ref Schoofbasaltopg2003] by a Jensen's inequality argument.  (See also the convexity argument at the bottom of this page.)


\section application Practical application, and Taylor approximation

The above formulas already reflect the recommendations Schoof gives on how to apply his formulas ([\ref Schoofbasaltopg2003], subsection 4.2).  The rest of this page is devoted to how the class PISMBedSmoother implements a practical version of this theory, based on these recommendations plus some additional approximation.

The averages appearing in his scaling arguments are over an infinite domain, %e.g.
  \f[ f_s(x) = \lim_{R\to\infty} \frac{1}{2R} \int_{-R}^R f(\xi,x)\,d\xi. \f]

For practical modeling use, Schoof specifically recommends averaging over some finite length scale which should be "considerably greater than the grid spacing, but much smaller than the size of the ice sheet."  Furthermore he recommends that, because of the typical aspect ratio of ice sheets, ""Bed topography on much larger length scales than 10 km should then be resolved explicitly through the smoothed bed height \f$b_s\f$ rather than the correction factor \f$\theta\f$.""  Thus in PISM we use \f$\lambda_1 = \lambda_2 = 5\f$ km as the default.  Naturally the values are configurable also.

It is, of course, possible to have bed roughness of significant magnitude at essentially any wavelength.  We make no claim that PISM results are good models of ice flow over \e arbitrary geometry; clearly the current models cannot come close to the non-shallow solution (Stokes) in such cases.  Rather, the goal right now is to improve on the existing shallow models, the diffusive SIA specifically, while maintaining or increasing high-resolution performance and comprehensive model quality, which necessarily includes many other modeled physical processes like ice thermal state, basal lubrication, and so on.  The desirable properties of the Schoof scheme are accepted not because the resulting model is perfect, but because we gain both a physical modeling improvement \e and a computational performance improvement from its use.

How do we actually compute expression (thetadefn) quickly?  Schoof has this suggestion, which we follow: ""To find \f$\theta\f$ for values of [surface elevation for which \f$\theta\f$ has not already been computed], some interpolation scheme should then be used.  \f$\theta\f$ is then represented at each grid point by some locally-defined interpolating function [of the surface elevation].""  

We need a "locally-defined interpolating function".  As with any approximation scheme, higher accuracy is achieved by doing "more work", which here is an increase in memory used for storing spatially-dependent coefficients.  Pade rational approximations, for example, were considered but are excluded because of the appearance of uncontrolled poles in the domain of approximation.  The 4th degree Taylor polynomial is chosen here because it shares the same convexity as the rational function it approximates; this is proven below.

Use of Taylor polynomial \f$P_4(x)\f$ only requires the storage of three fields (below), but it has been demonstrated to be reasonably accurate by trying beds of various roughnesses in Matlab/Octave scripts.  The derivation of the Taylor polynomial is most easily understood by starting with an abstract rational function like the one which appears in (thetadefn), as follows.

The fourth-order Taylor polynomial of the function \f$F(s)=(1-s)^{-k}\f$ around \f$s=0\f$ is
    \f[ P_4(s) = 1 + k s + \frac{k(k+1)}{2} s^2 + \dots + \frac{k(k+1)(k+2)(k+3)}{4!} s^4, \f]
so \f$F(s) = P_4(s) + O(s^5)\f$.  Let
    \f[ \omega(f,z) = \fint (1 - f(\xi) z)^{-k}\,d\xi \f]
where \f$f\f$ is some function and \f$z\f$ a scalar.  Then
\f{align*}
\omega(f,z) &= \fint (1 - f(\xi) z)^{-k}\,d\xi = \fint P_4(f(\xi) z)\,d\xi + O((\max |f(\xi)|)^5 |z|^5) \\
   &\approx 1 + k z \fint f(\xi)\,d\xi + \frac{k(k+1)}{2} z^2 \fint f(\xi)^2\,d\xi + \dots + \frac{k(k+1)(k+2)(k+3)}{4!} z^5 \fint f(\xi)^4\,d\xi.
\f}

Now, \f$\theta\f$ can be written
  \f[ \theta(h,x_1,x_2) = \left[ \omega(\tilde b(x_1,x_2,\cdot,\cdot), H^{-1}) \right]^{-n}.  \f]
So our strategy should be clear, to approximate \f$\omega(\tilde b(x_1,x2,\cdot,\cdot), H^{-1})\f$ by the Taylor polynomial as a function of \f$H^{-1}\f$, whose the coefficients depend on \f$x_1,x_2\f$.  We thereby get a rapidly-computable approximation for \f$\theta\f$ using stored coefficients which depend on \f$x_1,x_2\f$.  In fact, let \f$f(\xi) = \tilde b(x_1,x_2,\xi_1,\xi_2)\f$ for fixed \f$x_1,x_2\f$, and let \f$z=H^{-1}\f$.  Recall that the mean of this \f$f(\xi)\f$ is zero, so that the first-order term drops out in the above expansion of \f$\omega\f$.  We have the following approximation of \f$\theta\f$:
   \f[ \theta(h,x_1,x_2) \approx \left[ 1 + C_2(x_1,x_2) H^{-2} + C_3(x_1,x_2) H^{-3} + C_4(x_1,x_2) H^{-4} \right]^{-n} \tag{thetaapprox} \f]
where
  \f[ C_q(x_1,x_2) = \frac{k(k+1)\dots(k+q-1)}{q!} \fint \tilde b(x_1,x_2,\xi_1,\xi_2)^q\,d\xi_1\,d\xi_2 \f]
for \f$q=2,3,4\f$ and \f$k = (n+2)/n\f$.

Note that the coefficients \f$C_q\f$ depend only on the bed topography, and not on the ice geometry.  Thus we will pre-process the original bed elevations \f$b_0\f$ to compute and store the fields \f$b_s,C_2,C_3,C_4\f$.  The computation of equation (thetaapprox) is fast and easily-parallelized if these fields are pre-stored.  The computation of the coefficients \f$C_q\f$ and the smoothed bed \f$b_s\f$ at the pre-processing stage is more involved, especially when done in parallel.

The parameters \f$\lambda_1,\lambda_2\f$ must be set, but as noted above we use a default value of 5 km based on Schoof's recommendation.  This physical distance may be less than or more than the grid spacing.  In the case that the grid spacing is 1 km, for example, we see that there is a large smoothing domain in terms of the number of grid points.  Generally, the ghosting width (in PETSc sense) is unbounded.  Therefore move the unsmoothed topography to processor zero and do the smoothing and the coefficient-computing there.  The class PISMBedSmoother implements these details.

\subsection proveconvex Convexity of P_4

The approximation (thetaapprox) given above relates to the Jensen's inequality argument used by Schoof in Appendix A of [\ref Schoofbasaltopg2003].  For the nonsliding case, his argument shows that because \f$F(s)=(1-s)^{-k}\f$ is convex on \f$[-1,1]\f$ for \f$k>0\f$, therefore \f$0 \le \theta \le 1\f$.

Thus it is desirable for the approximation \f$P_4(z)\f$ to be convex on the same interval, and this is true.  In fact,
   \f[ P_4''(s) =  k(k+1) \left[1 + (k+2) s + \frac{(k+2)(k+3)}{2} s^2\right], \f]
and this function turns out to be positive for all \f$s\f$.  In fact we will show that the minimum of \f$P_4''(s)\f$ is positive.  That minimum occurs where \f$P_4'''(s)=0\f$ or \f$s=s_{min}=-1/(k+3)\f$.  It is a minimum because \f$P_4^{(4)}(s)\f$ is a positive constant.  And 
   \f[ P_4''(s_{min}) =  \frac{k(k+1)(k+4)}{2(k+3)} > 0. \f]


*/  /* <--- IF THIS END-OF-COMMENT IS LOST, CONFUSION RESULTS ... */

