
/*! @cond REFMAN  \page bombproof Appendix:  BOMBPROOF, a numerical scheme for temperature and age
\latexonly
\newtheorem{theorem}{Theorem}
\newcommand{\Up}[2]{\operatorname{Up}\left(#1\big|#2\right)}
\newcommand{\uppair}[2]{\left\{\begin{matrix} #1 \\ #2 \end{matrix}\right\}}

\index{BOMBPROOF!description}
One of the essential goals for any thermomechanically coupled numerical ice sheet model 
is a completely bombproof numerical 
scheme for the advection-conduction-reaction problem for the temperature within the ice.  
``Bombproof'' means, roughly, as stable as possible in as many realistic modeling
contexts as possible.  PISM illustrates this idea by including a scheme which is 
observed to be highly robust in practice.  In this section we prove stability 
properties of the scheme. 

Accuracy is also a goal, but it is necessarily a second goal.  It is partly achieved 
here because our scheme has truncation error $O(\Delta z^2)$
in many circumstances, though it reverts to lower order when it ``detects trouble''
in the form of large vertical velocities.  More completely stated, the scheme has 
order $O(\Delta t,\Delta x,\Delta y, \Delta z^2)$ in circumstances where the vertical
ice flow velocity is small enough.

The finite difference scheme in this section supercedes the PISM scheme for the 
temperature equation described in the appendices of \cite{BBL}.  That scheme 
is robust but it has a CFL condition on vertical advection, an avoidable time step 
restriction as a practical matter (because the scheme described here has been 
implemented successfully!).  Note the CFL condition on horizontal
remains in the new scheme, however.  The scheme here is conditionally stable 
as an adaptive time-stepping technique for the 3D temperature equation.
The length of the time step is limited only by the maximum magnitude of the 
horizontal velocities within the ice.

The region in which the temperature equation needs to be solved, namely within 
the ice sheet and attached ice shelves, is changing in time.  This is an 
essential complicating factor in ice sheet modeling.  According to user options, 
there may also be a model for the temperature within a layer of bedrock 
below the ice sheet.  The implementation must include this possibility, 
although it is not addressed here.

Equally essential is the fact that the velocity field has a
complicated and subtle provenance.  It comes from a variety of stress
balance equations, the nature of which is highly variable and poorly understood.
The stress balance could be the SIA, the SSA, a mix of these in different parts of the
ice sheet/shelf system, or a less shallow model.
In particular, we want to make no assumptions about the regularity of the velocity 
field in space or time.  Nor do we want the numerical
scheme for advection to need any information about the velocity except its value at the
beginning of the time step.  Thus we would not accept a scheme which required the Jacobian
of the velocity field with respect to changes in temperature, for example.
At very least we cannot implement a fully implicit scheme without some 
blind iteration (e.g.~with no guarantee of convergence of the iteration), while 
the scheme we propose involves no iteration.

As described elsewhere in this manual, the temperature equation is
\begin{equation}\label{basicConserve}
\rho c_p(T) \frac{dT}{dt} = k \frac{\partial^2 T}{\partial z^2} + \Sigma(T),
\end{equation}
where $T(t,x,y,z)$ is the temperature of the ice.  This equation is the shallow approximation
of the full three-dimensional conservation of energy, so it does not include horizontal
conduction \cite{Fowler}.  The units of this equation are energy per time, thus J/s in MKS.
We will assume, without exception, that $\rho,c_p(T),k$ are positive.
Note $dT/dt$ stands for the material derivative, so advection 
is present.  The function $\Sigma(T)$ is the strain-heating term, 
the details of which will not affect the derivation here.  The SI units of 
$\Sigma(T)$ are J$/(\text{s} \text{m}^3)$.  Note that $\Sigma(T)$ is naturally
described as a (self-)reaction term \cite{Fowler}.

We will be focusing on the dimension in which the temperature is changing fastest, 
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
The boundary conditions to problem \eqref{vertProblem} are a surface temperature
condition at the top of the ice and a geothermal flux condition at the bottom of the ice:
\begin{equation}\label{vPbcs}
T(t,z=H) = T_s(t,x,y), \qquad -k \frac{\partial T}{\partial z}\Big|_{z=0} = G(t,x,y).
\end{equation}
The case where there is a thermally-modeled layer of bedrock below the ice is not
considered here, though that case is implemented in 
\endlatexonly IceModel::temperatureStep() \latexonly .

For the discussion of the numerical scheme below, let $T_{ijk}^n$ be our approximation to 
the exact temperature $T$ at the grid point with coordinates $(x_i,y_j,z_k)$ at 
time $t_n$.  When $i,j$ are uninteresting we will suppress them and write $T_k^n$.  
We will use similar notation for numerical approximations to the other quantities, like
velocity components $u,v,w$ and $\Phi$.

We put the horizontal advection terms in the source term $\Phi$ because
(semi-)implicit treatment of horizontal advection would require a 
coupled system distributed across processors.  That kind of numerical difficulty
can be avoided, even given the other constraints on accuracy and implementation.
The following scheme for horizontal advection is explicit first-order upwinding.  
There is a CFL condition for the scheme to be stable, 
based on the magnitude of the horizontal velocity components.

Let\footnote{This is a slightly different definition of ``$\Up{f_{\bullet}}{\alpha_i}$''
from that used in the appendices of \cite{BBL}, though with exactly the same intent.}
   $$\Up{f_{\bullet}}{\alpha_i} = \begin{cases} f_i-f_{i-1}, & \alpha_i\ge 0, \\
                                                f_{i+1}-f_i, & \alpha_i<0.\end{cases}$$

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
See procedures \endlatexonly IceModel::temperatureStep() \latexonly and
\endlatexonly IceModel::ageStep(). \latexonly
In fact PISM has an optional remapping scheme (options \texttt{-quadZ} and 
\texttt{-chebZ}), where the temperatures/ages
in a column of ice are stored on an unequally-spaced vertical grid, but mapped to a fine,
equally-spaced grid for the temperature/age computation.  See the 
\endlatexonly IceModelVec3 \latexonly class and \endlatexonly 
IceModel::getVertLevsForTempAge(). \latexonly

\index{BOMBPROOF!statement of the numerical scheme}
Finally we get to the heart of the matter.  The $z$ derivative terms in \eqref{vertProblem}
will be approximated implicitly.  Let $\lambda$ be in the interval $0 \le \lambda \le 1$.  
The scheme BOMBPROOF is
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
while the first order upwind formula has lower accuracy,
	$$w_k^{n} \frac{\Up{T_{\bullet}^{n+1}}{w_k^{n}}}{\Delta z}
	   = w \frac{\partial T}{\partial z} + O(\Delta t,\Delta z).$$
We prefer for accuracy reasons to use the centered formula.  We will only include 
implicit upwinding when it is needed for its added stability benefits (below).

Let
	$$C_k^n = \rho c_p(T_k^n), \qquad\quad \nu = \frac{\Delta t}{\Delta z},
	 \qquad \text{and} \qquad R_k^n = \frac{k \Delta t}{C_k^n \Delta z^2}.$$
We now rewrite \eqref{bombone} for computational purposes as one of a system of equations 
for the unknowns $\{T_k^{n+1}\}$.  (In this system the coefficients are
scaled so that the diagonal entries of the matrix have limit one as
$\Delta t\to 0$.)  That is, we multiply equation \eqref{bombone} by
$\Delta t$, divide it by $C_k^n$, and rearrange:
\begin{align}
&\left(-R_k^n - \nu w_k^n \uppair{1-\lambda/2}{\lambda/2}\right) T_{k-1}^{n+1}  
   + \left(1 + 2 R_k^n + \nu w_k^n (1-\lambda) \uppair{+1}{-1}\right) T_k^{n+1} \notag \\
&\qquad\qquad + \left(-R_k^n + \nu w_k^n \uppair{\lambda/2}{1-\lambda/2} \right) T_{k+1}^{n+1}  = T_k^n + (C_k^n)^{-1}\Phi_k^n \label{bombtwo}
\end{align}
Here $\uppair{a}{b} = a$ when $w_k^n\ge 0$ and $\uppair{a}{b} = b$ when $w_k^n < 0$.

Equations \eqref{bombone} and \eqref{bombtwo} combine two unconditionally 
stable ways of approximating advection, namely implicit centered difference 
($\lambda = 1$) and implicit first-order upwinding ($\lambda=0$),
as note above.  The combination is ``convex'' because the coefficients are
nonnegative and sum to one.

The pure implicit centered difference scheme, namely
\begin{align}
&\left(-R_k^n - \nu w_k^n/2\right) T_{k-1}^{n+1} + \left(1 + 2 R_k^n\right) T_k^{n+1} \notag \\
&\qquad\qquad + \left(-R_k^n + \nu w_k^n/2\right) T_{k+1}^{n+1} 
              = T_k^n + (C_k^n)^{-1}\Phi_k^n \label{centered}
\end{align}
is \emph{less stable} than implicit first-order upwinding.  It is less stable 
in the same sense that Crank-Nicolson is a less stable scheme than 
backwards Euler for the simplest
heat equation $u_t = u_{xx}$ because, although oscillatory modes 
cannot grow exponentially they \emph{can} appear when
none are present already \cite{MortonMayers}, even in the homogeneous case $\Phi_k^n=0$.
One way of stating the greater stability of first-order upwinding is to say it
satisfies a ``maximum principle.''  On the other hand, as noted above, the centered difference
formula is desirable because it has higher order truncation error.

An example of a maximum principle for this kind of finite difference 
scheme is that if $U_{k-1}^n,U_k^n,U_{k+1}^n$ are adjacent gridded values at time step
$t_n$, and if the next value satisfies the scheme
	$$U_k^{n+1} = C_{-1} U_{k-1}^n + C_0 U_k^n + C_{+1} U_{k+1}^n$$
for \emph{nonnegative} coefficients $C_i$ summing to one ($C_{-1} + C_0 + C_{+1} = 1$), 
then it follows by the real triangle inequality that
	$$\min\{|U_{k-1}^n|, |U_k^n|, |U_{k+1}^n|\} 
	       \le |U_k^{n+1}| \le \max\{|U_{k-1}^n|, |U_k^n|, |U_{k+1}^n|\}.$$
Thus a ``wiggle'' cannot appear.  The proof below illustrates how this is handled 
for scheme \eqref{bombtwo}.

Thus we want to be precise about the phrase ``unconditionally stable'' for BOMBPROOF.
To do so we consider somewhat simplified cases which are amenable to analysis, and
we prove both particular stability properties we are claiming as the advantages of 
BOMBPROOF.

\index{BOMBPROOF!proof of stability properties}
\begin{theorem}  Assume (for the precise but limited assertion of this theorem) 
that the surface temperature $T_s$ and the geothermal flux $G$ are constant 
in time and space.  Assume also that the entire source function $\Phi$ is identically zero 
(but see comments below).

Fix an equally-spaced vertical grid $z_0=0 < z_1 < \dots < z_N=H$, so that 
the upper grid point coincides with the surface of the ice.

Suppose $c_p(T)$ is bounded below by a positive constant so that equation \eqref{bombtwo}
is a valid rewriting of equation \eqref{bombone}.  If
\begin{equation}\label{lambdachoice}
\lambda = \lambda^n = \min\left\{1, \quad 
	               \min_{k=0,\dots,N}\left\{\frac{2 k}{|w_k^n| C_k^n \Delta z}\right\}
	               \quad \right\}
\end{equation}
then scheme \eqref{bombone}/\eqref{bombtwo} is unconditionally stable 
in the following two senses:
\renewcommand{\labelenumi}{(\Roman{enumi})}\begin{enumerate}

\item  A maximum principle applies without further assumptions.
 
\item  Suppose we freeze the coefficients of the problem to have constant values.
(Concretely, we assume that the specific heat is temperature independent, so
$C_k^n = C_0$ and $R_k^n=R_0$ are positive constants, that 
$\lambda^n=\lambda_0$ is chosen independent of time step $n$, and that $\Delta t$ is 
the same for each time step.  We assume constant vertical velocity $w_k^n=w_0$.
We also consider a spatially periodic or unbounded version of our problem, 
with no boundary conditions.)
Then a von Neumann  analysis of the constant coefficient problem yields 
a growth factor less than one for all modes on the grid.
\end{enumerate}

These statements also apply in case $k=0$, in which case \eqref{lambdachoice} implies
$\lambda=0$, and the method reduces to implicit first-order upwinding.
\end{theorem}

\emph{Remark}.  The phrases \emph{maximum principle} and \emph{von Neumann analysis} 
will be precisely illustrated in the following proof.  Both approaches are explained 
from the beginning in \cite{MortonMayers}.  There is additional information on the 
von Neumann analysis of implicit methods for advection in \cite{Strikwerda}.

\begin{proof} \emph{Case (I)}.  In the case considered for the maximum principle, with $\Phi_k^n=0$, 
we can rewrite \eqref{bombtwo} as
\begin{align}\label{formformax}
&\left(1 + 2 R_k^n + \nu w_k^n (1-\lambda) \uppair{+1}{-1}\right) T_k^{n+1} \\
&\qquad = T_k^n + \left(R_k^n + \nu w_k^n \uppair{1-\lambda/2}{\lambda/2}\right) T_{k-1}^{n+1}
                + \left(R_k^n - \nu w_k^n \uppair{\lambda/2}{1-\lambda/2}\right) T_{k+1}^{n+1}.
                \notag
\end{align}

We claim that with choice \eqref{lambdachoice} for $0 \le \lambda \le 1$, all 
coefficients in \eqref{formformax} are nonnegative.  At one extreme, in 
the upwinding case ($\lambda=0$), all the coefficients are nonnegative.  Otherwise, note that
$\nu w_k^n (1-\lambda) \uppair{+1}{-1}$ is nonnegative for any valid value 
of $\lambda$ and for any value of $w_k^n$, noting the meaning of the ``$\uppair{+1}{-1}$'' 
symbol.  Thus the coefficient on the left is always nonnegative.  The coefficient of 
$T_{k-1}^{n+1}$ is clearly nonnegative for any valid value of $\lambda$ if $w_k^n \ge 0$.
The coefficient of $T_{k+1}^{n+1}$ is clearly nonnegative for any valid value of $\lambda$ if 
$w_k^n \le 0$.

Therefore the only concerns are for the coefficient of $T_{k-1}^{n+1}$ when $w_k^n\le 0$ and the 
coefficient of $T_{k+1}^{n+1}$ when $w_k^n\ge 0$.  But if $\lambda$ is smaller than 
$2k/(|w_k^n|C_k^n \Delta z)$ then 
\begin{align*}
R_k^n - \nu |w_k^n| (\lambda/2) &= \frac{k \Delta t}{C_k^n \Delta z^2}
	                                       - \frac{\Delta t |w_k^n|}{\Delta z} \frac{\lambda}{2}
	        &\ge \frac{k \Delta t}{C_k^n \Delta z^2}
	            - \frac{\Delta t |w_k^n|}{\Delta z} \frac{k}{|w_k^n|C_k^n \Delta z} = 0.
\end{align*}

Thus all the coefficients in \eqref{formformax} are positive, and a maximum principle applies.  
In particular, if we define the pointwise numerical error 
\cite{MortonMayers} $e_k^n = T_k^n - T(t_n,x_i,y_j,z_k)$, where $T(\dots)$ is the 
exact solution, then \eqref{formformax} implies an equality of the form
	$$A e_k^{n+1} = e_k^n + B_- e_{k-1}^{n+1} + B_+ e_{k+1}^{n+1} + \Delta t\, \tau_k^n$$
where $\tau_k^n$ is the truncation error of the scheme and $A,B_\pm$ are nonnegative 
coefficients.  Note that
	$$B_- + B_+ = 2 R_k^n + \nu w_k^n (1-\lambda) \uppair{+1}{-1}.$$

Letting $E^n = \max_k |e_k^n|$ we have, because of the positivity of coefficients,
\begin{align}
&\left(1 + 2 R_k^n + \nu w_k^n (1-\lambda) \uppair{+1}{-1}\right)|e_k^{n+1}|  \label{prebound}\\
& \qquad\qquad \le E^n + \left(2 R_k^n + \nu w_k^n (1-\lambda) \uppair{+1}{-1}\right)E^{n+1}
     + \Delta t\,\bar\tau^n \notag
\end{align}
for all $k$, where $\bar\tau^n = \max_k |\tau_k^n|$.  Now let $k$ be the index for 
which $|e_k^{n+1}| =E^{n+1}$.  For that $k$ we can replace ``$|e_k^{n+1}|$'' in 
equation \eqref{prebound} with ``$E^{n+1}$''.  Subtracting the same quantity from 
each side of the resulting inequality gives
\begin{equation*}
E^{n+1} \le E^n + \Delta t\,\bar\tau^n,
\end{equation*}
Thus have a maximum principle for BOMBPROOF, and have used it to prove convergence 
in the standard way \cite{MortonMayers}.

\noindent \emph{Case (II)}.  As a von Neumann analysis is much more restrictive 
than the analysis above, we will be brief here.  
Let's assume the velocity is downward, $w_0<0$; the other case 
is very similar.  Equation \eqref{bombtwo} becomes 
\begin{align}
&\left(-R_0 - \nu w_0 (\lambda/2)\right) T_{k-1}^{n+1}  
   + \left(1 + 2 R_0 - \nu w_0 (1-\lambda)\right) T_k^{n+1} \label{prevonN}\\
&\qquad\qquad + \left(-R_0 + \nu w_0 (1-\lambda/2) \right) T_{k+1}^{n+1}  = T_k^n.\notag
\end{align}

The heart of the von Neumann analysis is the substitution of a growing or decaying
(with $n$) oscillatory mode on the grid of spatial wave number $\mu$:
	$$T_k^n = \sigma^n e^{i\mu\,(k\Delta z)}.$$
Note $k\Delta z = z_k$ is a grid point.  Such a mode is a solution to \eqref{prevonN} 
if and only if
\begin{align*}
\sigma\Big[  &(-R_0-\nu w_0(\lambda/2)) e^{-i\mu\Delta z}
	              + (1+2R_0 - \nu w_0 (1-\lambda)) \\
	           &\quad   + (-R_0  + \nu w_0 (\lambda/2)) e^{+i\mu\Delta z}
	                    + \nu w_0 (1-\lambda) e^{+i\mu\Delta z} \Big] = 1.
\end{align*}
This equation reduces by standard manipulations to
	$$\sigma = \frac{1}{1 + \left(4 R_0 - 2 \nu w_0 (1-\lambda)\right)\cos^2(\mu \Delta z/2)
	                      + i\,\nu w_0 (1-\lambda/2)\sin(\mu\Delta z)}.$$
Note $4 R_0 - 2 \nu w_0 (1-\lambda) \ge 0$ without restrictions on 
numerical parameters $\Delta t$, $\Delta z$, because $w_0<0$ in the 
case under consideration.  Therefore
	$$|\sigma|^2 = \frac{1}{\left[1 + \left(4 R_0 - 2 \nu w_0 (1-\lambda)\right)
	                              \cos^2(\mu \Delta z/2)\right]^2
	                      + \left[\nu w_0 (1-\lambda/2)\sin(\mu\Delta z)\right]^2} < 1,$$
and it follows that all modes decay exponentially.
\end{proof}

Note that $\lambda$ is carefully chosen in \eqref{lambdachoice} so that the 
maximum principle \emph{(I)} applies.  On the other hand, both the 
implicit first-order upwind and
the implicit centered difference formulas have unconditional stability in the von Neumann
sense.  The proof of case \emph{(II)} above is thus a formality, just showing that a
convex combination of unconditionally stable (von Neumann sense) schemes is still
unconditionally stable in the same sense.

Recall we assumed in Theorem 1 that the entire ``source'' $\Phi_k^n$ was identically zero.
Of course this is not realistic.  What we understand is provable, however, is that if a
numerical scheme for a linear advection/conduction equation
	$$u_t + A u_x = B u_{xx}$$
is stable in the general sense of numerical schemes for partial differential equations
(e.g.~as defined in subsection 5.5 of \cite{MortonMayers}) then the same scheme is stable 
in the same general sense when applied to the equation with (linear) lower order terms:
	$$u_t + A u_x = B u_{xx} + C u + D.$$
A precise statement of this general fact is hard to find in the literature, 
to put it mildly, but theorem 2.2.3 of \cite{Strikwerda} is one interesting case
($B=0$ and $D=0$).

\endlatexonly
@endcond
 */ 

