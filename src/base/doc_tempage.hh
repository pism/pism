

/*! \fn PetscErrorCode IceModel::temperatureStep()
    \brief Takes a semi-implicit time-step for the temperature equation.

In summary, the conservation of energy equation is
    \f[ \rho c_p \frac{dT}{dt} = k \frac{\partial^2 T}{\partial z^2} + \Sigma,\f] 
where \f$T(t,x,y,z)\f$ is the temperature of the ice.  This equation is the shallow approximation
of the full three-dimensional conservation of energy.  Note \f$dT/dt\f$ stands for the material
derivative, so advection is included.  Here \f$\rho\f$ is the density of ice, 
\f$c_p\f$ is its specific heat, and \f$k\f$ is its conductivity.  Also \f$\Sigma\f$ is the volume
strain heating.

In summary, the numerical method is first-order upwind for advection and centered-differences with
semi-implicitness for the vertical conduction term.  We work from the bottom 
of the column upward in building the system to solve (in the semi-implicit time-stepping scheme).
The excess energy above pressure melting is converted to melt-water, and that a fraction 
of this melt water is transported to the base according to the scheme in excessToFromBasalMeltLayer().

The method uses equally-spaced calculation but the methods getValColumn(), setValColumn() interpolate 
back and forth from this equally-spaced calculational grid to the (usually) non-equally space storage 
grid.

In this procedure four scalar fields are modified: vHmelt, vbasalMeltRate, Tb3, and Tnew3.
But vHmelt, vbasalMeltRate and Tb3 will never need to communicate ghosted values (i.e. horizontal 
stencil neighbors.  The ghosted values for T3 are updated from the values in Tnew3 in the
communication done by temperatureAgeStep().
 
@cond REFMAN
Here is a more complete discussion and derivation.

Consider a column of a slowly flowing and heat conducting material as shown in the left side 
of the next figure.  (The left side shows a general column of flowing and heat conduction 
material showing a small segment \f$V\f$.  The right side shows a more specific column of ice 
flowing and sliding over bedrock.)  This is an \em Eulerian view so the material flows through 
a column which remains fixed (and is notional).  The column is vertical.  We will 
assume when needed that it is rectangular in cross-section with cross-sectional area 
\f$\Delta x\Delta y\f$.

\image latex earlycols.png "Left: a general column of material.  Right: ice over bedrock." width=3in

[FIXME: CONTINUE TO MINE eqns3D.tex FOR MORE]

The application of the geothermal flux at the base of a column is a special case for which 
we give a finite difference argument.  This scheme follows the equation (2.114) in 
Morton and Mayers (2005).  We have the boundary condition
	\f[  -k \frac{\partial T}{\partial z} = G(t,x,y) \f]
where \f$G(t,x,y)\f$ is the applied geothermal flux, and it is applied at level \f$z=-B_0\f$ 
in the bedrock (which is the only case considered here).  We <em> add a virtual lower grid 
point </em> \f$z_{-1} = z_0 - \Delta  z\f$ and we approximate the above boundary condition 
at $z_0$ by the centered-difference
	\f[  -k \frac{T_{1} - T_{-1}}{2 \Delta z} = G. \f]
Here \f$T_k = T_{ijk}^{l+1}\f$.  We also apply the discretized conduction equation at $z_0$.  These two combined equations 
yield a simplified form
	\f[(1 + 2 KR) T_0 - 2 K T_1 = T_0 + \frac{2\Delta t}{\rho c_p \Delta z} G \f]
where \f$K = k \Delta t (\rho c \Delta z^2)^{-1}\f$.

@endcond
 */
 
 

/*! \page bombproof Appendix:  BOMBPROOF, a proposed numerical scheme for temperature and age
\latexonly

One of the essential goals for a thermomechanically coupled numerical ice sheet model 
is a completely bombproof numerical\footnote{Here, finite difference, of course.} 
scheme for the advection-conduction-reaction problem for the temperature within the ice.  
``Bombproof'' means, roughly, as stable as possible in as many ways as possible, as 
we explain below.  Accuracy is also a goal, but it is necessarily a second goal,
partly achieved here.

The scheme in this section supercedes the numerical scheme for the temperature equation 
described in \cite{BBL}.  That scheme is robust but has a CFL condition on vertical
velocity, an avoidable time step restriction.

The region in which the temperature equation 
needs to be solved, namely within the ice sheet and attached ice shelves,
is changing in time.  This is an essential complicating factor in ice sheet modeling.
Furthermore, according to user options, there may be a model
for the temperature within a layer of bedrock below the ice sheet.

Equally essential is the fact that the velocity field has an
extraordinarily complicated and subtle provenance.  It comes from a variety of stress
balance equations, the nature of which is highly variable and poorly understood; 
the stress balance could be the SIA, the SSA, a mix of these, or a less shallow model.
In particular, we want to make no assumptions about the regularity of the velocity 
field in space or time, nor do we want the numerical
scheme for advection to need any information about the velocity except its value at the
beginning of the time step.  Thus we would not accept a scheme which required the Jacobian
of the velocity field with respect to changes in temperature, for example.  We cannot
implement a fully implicit scheme without some blind iteration (e.g.~with
no guarantee of convergence of the iteration), while the scheme we propose 
involves no iteration.

As described elsewhere in this manual, the temperature equation is
\begin{equation}\label{basicConserve}
\rho c_p(T) \frac{dT}{dt} = k \frac{\partial^2 T}{\partial z^2} + \Sigma(T),
\end{equation}
where $T(t,x,y,z)$ is the temperature of the ice.  This equation is the shallow approximation
of the full three-dimensional conservation of energy, and it does not include horizontal
conduction \cite{Fowler}.  The units of this equation are energy per time, thus J/s in MKS.
Note $dT/dt$ stands for the material derivative, so advection 
is included.  The function $\Sigma(T)$ is the strain-heating term, 
the details of which will not affect the derivation here.  Note that $\Sigma(T)$ may be
described as a (self-)reaction term, in analogy with standard language for fluids 
in which chemical reactions are occurring \cite{Fowler}.

We will be focussing on the dimension in which the temperature is changing fastest, 
namely the vertical coordinate $z$.  Rewriting equation \eqref{basicConserve} 
to emphasize the vertical terms we have
\begin{align}\label{vertProblem}
\rho c_p(T)\left(\frac{\partial T}{\partial t} + w \frac{\partial T}{\partial z}\right) 
     &= k \frac{\partial^2 T}{\partial z^2} + \Phi
\end{align}
where
\begin{equation*}
\Phi = \Sigma(T) - \rho c_p(T) \left(u \frac{\partial T}{\partial x}
                         + v \frac{\partial T}{\partial y}\right)
\end{equation*}

For the discussion of the numerical scheme below, let $T_{ijk}^n$ be our approximation to 
the exact temperature $T$ at the grid point with coordinates $(x_i,y_j,z_k)$ at 
time $t_n$.  When $i,j$ are uninteresting we will suppress them and write $T_k^n$.  
We will use similar notation for numerical approximations to the other quantities, like
velocity components $u,v,w$ and also $\Sigma$ or $\Phi$ above.

\newcommand{\Up}[2]{\operatorname{Up}\left(#1\big|#2\right)}

We include the horizontal advection terms explicitly because (semi-)implicit treatment of 
these terms would require a coupled system distributed across processors, and that kind of 
system can be avoided given the other constraints on accuracy and implementation.
In fact, as follows, the scheme for these terms 
is explicit first-order upwinding.  There is a CFL condition for the scheme to be stable, 
based on the magnitude of the horizontal velocity components.
Let $\Up{f_{\bullet}}{\alpha_i} = \alpha_i(f_i-f_{i-1})$ 
if $\alpha_i\ge 0$ while $\Up{f_{\bullet}}{\alpha_i}=\alpha_i(f_{i+1}-f_i)$ if $\alpha_i<0$.
The horizontal advection terms, part of $\Phi$, are approximated
	$$\Phi_{ijk}^n = \Sigma(T_{ijk}^n) - \rho c_p(T_{ijk}^n) 
	                   \left( \frac{\Up{T_{\bullet jk}^n}{u_{ijk}^n}}{\Delta x}
	                          + \frac{\Up{T_{i\bullet k}^n}{v_{ijk}^n}}{\Delta y} \right).$$
The CFL condition for the stability of this scheme is
\begin{equation}\label{CFL}
    \Delta t \,\left( \left|\frac{u_{ijk}^n}{\Delta x}\right|
                           + \left|\frac{v_{ijk}^n}{\Delta y}\right| \right) \le 1.
\end{equation}
The routine \endlatexonly IceModel::computeMax3DVelocities() \latexonly computes the 
maximum of velocity magnitudes.  This produces a time step restriction based on 
the above CFL condition.  Note \endlatexonly IceModel::determineTimeStep() \latexonly 
implements adaptive time-stepping based on this and other stability criteria.

We will assume an equally-spaced grid $z_0,\dots,z_{M_z}$ with $\Delta z = z_{k+1} - z_k$.
In fact PISM has an optional remapping scheme (options \texttt{-quadZ} and \texttt{-chebZ}), 
where the temperatures/ages
in a column of ice are stored on an unequally-spaced vertical grid, but mapped to a fine
equally spaced grid when we do a temperature/age time step.  See the 
\endlatexonly IceModelVec3 \latexonly class.

Finally we get to the heart of the matter.  The $z$ derivative terms in \eqref{vertProblem}
will be approximated implicitly.  Let $\lambda$ be in the interval $0 \le \lambda \le 1$.  The scheme BOMBPROOF is
\begin{align}\label{bombone}
\rho c_p(T_k^n) &\left( 
       \frac{T_k^{n+1} - T_k^n}{\Delta t}
       + \lambda w_k^n \frac{T_{k+1}^{n+1} - T_{k-1}^{n+1}}{2 \Delta z}
       + (1-\lambda) \frac{\Up{T_{\bullet}^{n+1}}{w_k^{n}}}{\Delta z} \right) \\
     &= k\, \frac{T_{k+1}^{n+1} - 2 T_{k}^{n+1} + T_{k-1}^{n+1}}{\Delta z^2} + \Phi_k^n \notag
\end{align}
with the value for $\lambda$ determined by equation ??? below.  Here, as promised, 
we have suppressed indices $i,j$.

Equation \eqref{bombone} has two approximations of the vertical advection, 
combined using nonnegative coefficient which sum to one.  Both approximations
are (semi-)implicit, but the centered formula has higher accuracy,
	$$w_k^n \frac{T_{k+1}^{n+1} - T_{k-1}^{n+1}}{2 \Delta z}
	   = w \frac{\partial T}{\partial z} + O(\Delta t,\Delta z^2),$$
while the first order upwind formula has less,
	$$\frac{\Up{T_{\bullet}^{n+1}}{w_k^{n+1}}}{\Delta z}
	   = w \frac{\partial T}{\partial z} + O(\Delta t,\Delta z).$$
We prefer for accuracy reasons to use the centered formula.  We will only include 
upwinding when it is needed for its added stability benefits, as explained next
in the context of two definitions of stability.

Let
	$$R = \frac{k \Delta t}{\Delta z^2}$$
and let
	$$\nu = \frac{\Delta t}{\Delta z}.$$

There are two unconditionally stable ways


\endlatexonly
 */ 


