

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
is a completely bombproof numerical 
scheme for the advection-conduction-reaction problem for the temperature within the ice.  
``Bombproof'' means, roughly, as stable as possible in as many ways as possible, as 
we explain below.  Accuracy is also a goal, but it is necessarily a second goal,
partly achieved here.

The finite difference scheme in this section supercedes 
the scheme for the temperature equation 
described in \cite{BBL}.  That scheme is robust but has a CFL condition on vertical
advection, an avoidable time step restriction.  Note the CFL condition on horizontal
remains in the new scheme.

The region in which the temperature equation 
needs to be solved, namely within the ice sheet and attached ice shelves,
is changing in time.  This is an essential complicating factor in ice sheet modeling.
According to user options, there may also be a model
for the temperature within a layer of bedrock below the ice sheet.

Equally essential is the fact that the velocity field has an
complicated and subtle provenance.  It comes from a variety of stress
balance equations, the nature of which is highly variable and poorly understood; 
the stress balance could be the SIA, the SSA, a mix of these in different parts of the
ice sheet/shelf system, or a less shallow model.
In particular, we want to make no assumptions about the regularity of the velocity 
field in space or time.  Nor do we want the numerical
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
of the full three-dimensional conservation of energy, so it does not include horizontal
conduction \cite{Fowler}.  The units of this equation are energy per time, thus J/s in MKS.
Note $dT/dt$ stands for the material derivative, so advection 
is included.  The function $\Sigma(T)$ is the strain-heating term, 
the details of which will not affect the derivation here.\footnote{The units of 
$\Sigma(T)$ in equation \eqref{basicConserve}, and in the rest of this Appendix,
are, however, different from those of 
$\Sigma$ elsewhere in this note.  FIXME.}  Note that $\Sigma(T)$ may be
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
The following scheme for these terms is explicit first-order upwinding.  
There is a CFL condition for the scheme to be stable, 
based on the magnitude of the horizontal velocity components.

Let\footnote{This is a slightly different definition of ``$\Up{f_{\bullet}}{\alpha_i}$''
from that used in the appendices of \cite{BBL}, though with exactly the same intent.}
$\Up{f_{\bullet}}{\alpha_i} = f_i-f_{i-1}$ 
if $\alpha_i\ge 0$ while $\Up{f_{\bullet}}{\alpha_i}=f_{i+1}-f_i$ if $\alpha_i<0$.

The horizontal advection terms, thus $\Phi$, are approximated
	$$\Phi_{ijk}^n = \Sigma(T_{ijk}^n) - \rho c_p(T_{ijk}^n) 
	                   \left( u_{ijk}^n\,\frac{\Up{T_{\bullet jk}^n}{u_{ijk}^n}}{\Delta x}
	                          + v_{ijk}^n\,\frac{\Up{T_{i\bullet k}^n}{v_{ijk}^n}}{\Delta y} \right).$$
Thus $\Phi$ is approximated in a completely explicit way, at time $t_n$ at the beginning 
of the time step.  The CFL stability condition implied by this part of the scheme is
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
in a column of ice are stored on an unequally-spaced vertical grid, but mapped to a fine,
equally-spaced grid for the temperature/age time step.  See the 
\endlatexonly IceModelVec3 \latexonly class and procedures \endlatexonly 
IceModel::temperatureStep(), \latexonly \endlatexonly IceModel::ageStep(), \latexonly
and \endlatexonly IceModel::getVertLevsForTempAge(). \latexonly

Finally we get to the heart of the matter.  The $z$ derivative terms in \eqref{vertProblem}
will be approximated implicitly.  Let $\lambda$ be in the interval $0 \le \lambda \le 1$.  The scheme BOMBPROOF is
\begin{align}\label{bombone}
\rho c_p(T_k^n) &\left( 
       \frac{T_k^{n+1} - T_k^n}{\Delta t}
       + \lambda w_k^n \frac{T_{k+1}^{n+1} - T_{k-1}^{n+1}}{2 \Delta z}
       + (1-\lambda) w_k^{n} \frac{\Up{T_{\bullet}^{n+1}}{w_k^{n}}}{\Delta z} \right) \\
     &= k\, \frac{T_{k+1}^{n+1} - 2 T_{k}^{n+1} + T_{k-1}^{n+1}}{\Delta z^2} + \Phi_k^n \notag
\end{align}
with the value for $\lambda$ determined by \eqref{lambdachoice} below.  As promised, 
from now on we suppress indices $i,j$.

Equation \eqref{bombone} has two approximations of vertical advection which are 
combined using nonnegative coefficients which sum to one.  Both approximations
are (semi-)implicit, but the centered formula has higher accuracy,
	$$w_k^n \frac{T_{k+1}^{n+1} - T_{k-1}^{n+1}}{2 \Delta z}
	   = w \frac{\partial T}{\partial z} + O(\Delta t,\Delta z^2),$$
while the first order upwind formula has less,
	$$w_k^{n} \frac{\Up{T_{\bullet}^{n+1}}{w_k^{n}}}{\Delta z}
	   = w \frac{\partial T}{\partial z} + O(\Delta t,\Delta z).$$
We prefer for accuracy reasons to use the centered formula.  We will only include 
upwinding when it is needed for its added stability benefits (below).

Let
	$$C_k^n = \rho c_p(T_k^n), \qquad\quad \nu = \frac{\Delta t}{\Delta z},
	 \qquad \text{and} \qquad R_k^n = \frac{k \Delta t}{C_k^n \Delta z^2}.$$
The following equation is the same as \eqref{bombone} but rewritten for computational 
purposes as a line in a system of equations for the unknowns $\{T_k^{n+1}\}$, with 
coefficients scaled so that the diagonal entries of the matrix have limit one as
$\Delta t\to 0$.  That is, we take equation \eqref{bombone} and multiply it by
$\Delta t$ divide it by $C_k^n$, and rearrange:
\newcommand{\uppair}[2]{\left\{\begin{matrix} #1 \\ #2 \end{matrix}\right\}}
\begin{align}\label{bombtwo}
&\left(-R_k^n - \nu w_k^n \uppair{1-\lambda/2}{\lambda/2}\right) T_{k-1}^{n+1}  \\
&\qquad\qquad + \left(1 + 2 R_k^n + \nu w_k^n (1-\lambda) \uppair{+1}{-1}\right) T_k^{n+1} \notag \\
&\qquad\qquad + \left(-R_k^n + \nu w_k^n \uppair{\lambda/2}{1-\lambda/2} \right) T_{k+1}^{n+1} \notag \\
&\quad = T_k^n + (C_k^n)^{-1}\Phi_k^n \notag
\end{align}
Here $\uppair{a}{b} = a$ when $w_k^n\ge 0$ and $\uppair{a}{b} = b$ when $w_k^n < 0$.

Equation \eqref{bombtwo} combines two unconditionally stable ways of approximating advection,
namely implicit centered and implicit first-order upwinding as note above.  We want to
be precise about the phrase ``unconditionally stable,'' however.  To do so we assume particular
realistic boundary conditions to problem \eqref{vertProblem}, namely a surface temperature
condition at the top of the ice and a geothermal flux condition at the bottom of the ice:
\begin{equation}\label{vPbcs}
T(t,z=H) = T_s(t,x,y), \qquad -k \frac{\partial T}{\partial z}\Big|_{z=0} = G(t,x,y).
\end{equation}

\newtheorem{theorem}{Theorem}

\begin{theorem}  Assume (for the precise but limited assertion of this theorem) 
that the surface temperature $T_s$ and the geothermal flux $G$ are constant 
in time and space.  Assume also that the entire source function $\Phi$ is identically zero 
(but see comments below).

Fix an equally-spaced vertical grid $z_0=0 < z_1 < \dots < z_N=H$, so that 
the upper grid point coincides with the surface of the ice.

If
\begin{equation}\label{lambdachoice}
\lambda = \lambda^n = \min\left\{1, \quad 
	               \min_{k=0,\dots,N}\left\{\frac{2 k}{w_k^n C_k^n \Delta z}\right\}
	               \quad \right\}
\end{equation}
then scheme \eqref{bombone}/\eqref{bombtwo} is unconditionally stable 
in the following two senses
\begin{enumerate}
\item  A maximum principle applies under the minimal additional assumption that 
$c_p(T)$ is bounded below by a positive constant.  
\item  Suppose we freeze the coefficients of the problem to make it constant coefficient.
(Concretely, we assume that the specific heat is temperature independent, that
$C_k^n = C_0$ and $R_k^n=R_0$ are positive constants, that 
$\lambda^n=\lambda_0$ is chosen independent of time step $n$, and that $\Delta t$ is 
the same for each time step.  We also consider a spatially periodic or 
unbounded version of our problem, with no boundary conditions.)
Then a von Neumann  analysis of the constant coefficient problem yields 
a growth factor less than one for all modes on the grid.
\end{enumerate}
\end{theorem}

\emph{Remark}.  The phrases \emph{maximum principle} and \emph{von Neumann analysis} 
will be precisely illustrated in the following proof.  Both approaches are well-explained in 
\cite{MortonMayers}.  There is additional information on the von Neumann analysis
of implicit methods for advection in \cite{Strikwerda}.

\begin{proof}  \end{proof}



\endlatexonly
 */ 


