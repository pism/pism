// Copyright (C) 2008--2010 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <petscsnes.h>
#include "iceModelVec.hh"
#include "iceModel.hh"


/* THIS SOURCE FILE CONTAINS THE TWO LEAST TRIVIAL ROUTINES IN THE INVERSE
   MODEL COMPONENTS OF IceModel, NAMELY THE TWO WHICH BUILD PETSc Mat AND SNES
   OBJECTS.  SEE src/base/iMinverse.cc FOR OTHER INVERSE MODEL-RELATED ROUTINES.  */


//! Compute basal shear stress from a given sliding velocity using SSA equations.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

We assume the basal sliding velocity is known, though this ``basal'' velocity
is usually the vertically-average horizontal velocity computed by subtracting the
nonsliding SIA part of an observed surface velocity.

Using ice thickness, ice surface elevation, and an ice temperature field, 
we use SSA equations on the to compute the basal shear stress.  There is no 
smoothing in this implementation.  In particular:
\latexonly
\def\ddt#1{\ensuremath{\frac{\partial #1}{\partial t}}}
\def\ddx#1{\ensuremath{\frac{\partial #1}{\partial x}}}
\def\ddy#1{\ensuremath{\frac{\partial #1}{\partial y}}}
\begin{align*}
  - 2 \ddx{}\left[\nu H \left(2 \ddx{u} + \ddy{v}\right)\right]
        - \ddy{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
        + \tau_{(b)x}  &=  - \rho g H \ddx{h}, \\
    - \ddx{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
      - 2 \ddy{}\left[\nu H \left(\ddx{u} + 2 \ddy{v}\right)\right]
        + \tau_{(b)y}  &=  - \rho g H \ddy{h}, 
\end{align*}
\endlatexonly
Note \f$\nu\f$ is temperature-dependent.

In this case the known field is \f$(u,v)\f$.  We are solving for 
\f$\tau_b = (\tau_{(b)x},\tau_{(b)y})\f$.  The equations are solved everywhere, 
regardless of whether there were observed surface velocities or even positive
ice thickness.  Because the SSA differential operator only reads in a 9 point
stencil, the condition (invMask == 2) determines whether the result of this
routine could be meaningful.

This procedure calls
- assembleSSAMatrix(),
- assembleSSARhs(),
- computeEffectiveViscosity().

It saves temporary versions of various flags modified by those routines, to
avoid changing the state of IceModel.
 */
PetscErrorCode IceModel::computeBasalShearFromSSA() {
  PetscErrorCode ierr;

  // effective viscosity for ub_in, vb_in: save flags before changing them
  const bool leaveNuHAloneSSA_save = leaveNuHAloneSSA, 
    useConstantNuHForSSA_save = config.get_flag("use_constant_nuh_for_ssa");
  leaveNuHAloneSSA = PETSC_FALSE;
  config.set_flag("use_constant_nuh_for_ssa", false);
  IceModelVec2S myvNuH[2] = {vWork2d[0], vWork2d[1]}; // already allocated space
  // eps=0.0 in bdd-below regularization; Schoof type regularization does occur;
  ierr = computeEffectiveViscosity(myvNuH, 0.0); CHKERRQ(ierr); // uses ssavel
  leaveNuHAloneSSA = leaveNuHAloneSSA_save;
  config.set_flag("use_constant_nuh_for_ssa", useConstantNuHForSSA_save);

  // allocate Mat and Vecs; compare setup in IceModel::createVecs() and
  //   IceModel::velocitySSA()
  Mat A;
  Vec x, result, rhs;
  const PetscInt M = 2 * grid.Mx * grid.My;
  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE, M, M,
                         13, PETSC_NULL, 13, PETSC_NULL,
                         &A); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &result); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &rhs); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &x); CHKERRQ(ierr);

  // Build approximation of SSA equations (A x = rhs)
  //   Generally we want basal shear as part of matrix:
  //                L[U] + beta U = (driving),
  //   so  A = L + beta I. In this case we don't want basal shear terms 
  //   as part of the matrix:
  //                tau_b := L[U] - (driving),
  //   so  A = L.
  ierr = assembleSSAMatrix(false, myvNuH, A); CHKERRQ(ierr);
 
  // Note rhs contains driving terms  - \rho g H \grad h
  ierr = assembleSSARhs(rhs); CHKERRQ(ierr);

  // Set x = [u v]'  (interleaved).
  const PetscInt  twoMy = 2 * grid.My;
  ierr = VecSet(x, 0.0); CHKERRQ(ierr);

  ierr = ssavel.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = VecSetValue(x, i*twoMy + 2*j,   ssavel(i,j).u, INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(x, i*twoMy + 2*j+1, ssavel(i,j).v, INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = ssavel.end_access(); CHKERRQ(ierr);

  // assemble matrix and vec before multiplication; communicate!
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

  // major action: With new A and x, compute  result = Ax - rhs = L[u] - (driving)
  ierr = MatMult(A,x,result); CHKERRQ(ierr);  
  
  ierr = VecAXPY(result,-1.0,rhs); CHKERRQ(ierr);  // note result = 0 if MASK_SHEET
  
  // transfer result to taub output vectors; compare code in
  //   IceModel::moveVelocityToDAVectors(); only transfer if full 9 pt
  //   neighborhood had valid velocities
  PetscScalar     **tbx, **tby, *res;
  Vec             resultLoc = SSAXLocal;
  ierr = VecScatterBegin(SSAScatterGlobalToLocal, result, resultLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(SSAScatterGlobalToLocal, result, resultLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGetArray(resultLoc, &res); CHKERRQ(ierr);
  ierr = inv.taubxComputed->get_array(tbx); CHKERRQ(ierr);
  ierr = inv.taubyComputed->get_array(tby); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // note that computed taub should only be used to get new phi at
      //   points where inv.taubValidMask is 1.0
      tbx[i][j] = res[i*twoMy + 2*j];
      tby[i][j] = res[i*twoMy + 2*j+1];
    }
  }
  ierr = inv.taubxComputed->end_access(); CHKERRQ(ierr);
  ierr = inv.taubyComputed->end_access(); CHKERRQ(ierr);
  ierr = VecRestoreArray(resultLoc, &res); CHKERRQ(ierr);

  // de-allocate
  ierr = MatDestroy(A); CHKERRQ(ierr);
  ierr = VecDestroy(x); CHKERRQ(ierr);
  ierr = VecDestroy(result); CHKERRQ(ierr);
  ierr = VecDestroy(rhs); CHKERRQ(ierr);

  return 0;
}


extern PetscErrorCode RegPoissonFunctionLocal(
                DALocalInfo *info, PetscScalar **x, PetscScalar **F,
                RegPoissonCtx *user);


//! Compute the till friction angle for the given basal shear stress by minimizing a regularized objective.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

A more complete but still brief description might be:  Compute \f$\mu = \tan\phi\f$,
where \f$\phi\f$ is the till friction angle, from the basal shear stress map 
and the basal velocity (both vector-valued) which come from the observed surface 
velocity.  The resulting \f$\mu\f$ map minimizes a combination
of a smoothness norm and the norm of the difference between the basal
shear stress from the observed velocity and that from the pseudo-plastic 
model using \f$\mu\f$. The smoothness norm is the square of the gradient of 
the difference between \f$\mu\f$ and the value of \f$\tan\phi=\mu_i\f$
for the IceModel state at the time the inverse model is invoked.

This procedure should not be called in the purely plastic case.  It requires 
IceBasalResistancePlasticLaw::pseudo_plastic to be \c TRUE, which occurs if the 
<tt>-pseudo_plastic_q</tt> option is used with a positive value.  Recall that 
in this case the method IceBasalResistancePlasticLaw::drag() computes the basal shear stress as
    \f[ \tau_b = - \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U} \f]
where \f$\tau_b=(\tau_{(b)x},\tau_{(b)y})\f$, \f$U=(u,v)\f$,
\f$q=\f$ <tt>pseudo_q</tt>, and \f$U_{\mathtt{th}}=\f$ <tt>pseudo_u_threshold</tt>.
Typical values for the constants are \f$q=0.25\f$ and \f$U_{\mathtt{th}} = 100\f$
m/a.  The linearly viscous till case \f$q=1\f$ is allowed, in which case 
\f$\beta = \tau_c/U_{\mathtt{th}}\f$.

The rest of this comment is a sketch of the procedure, including the 
regularization which leads to a Poisson-like problem and the use of 
a PETSc solver.

Recall the Mohr-Coulomb criterion for till yield stress, where, assuming zero 
till cohesion,
	\f[\tau_c = (\tan\phi) N = \mu N.\f]
Here \f$N\f$ is the effective pressure on the till which is known.  It is a
function of overburden pressure \f$\rho g H\f$ and of modeled pore water pressure.

Let
	\f[ \mathbf{V} = \frac{\mathbf{U}}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q}. \f]
Pseudo-plastic till satisfies \f$\tau_b = - \tau_c \mathbf{V}\f$.  In the current
procedure, \f$\mathbf{V}\f$ is a known vector field (in the map-plane).

We want to determine \f$\mu\f$ and \f$\tau_c\f$ given \f$\tau_b\f$, 
\f$\mathbf{U}\f$, and \f$N\f$.  But according to the above relations 
\f$\tau_b\f$ and \f$\mathbf{V}\f$ are supposed to be proportional vectors, while in this 
context they may not even be parallel.  (Otherwise we could merely divide in the 
equation for pseudo-plastic till.  Which is essentially what is done by
computeTFAFromBasalShearNoReg().)  Thus we first refine the goal to the 
minimization of the functional
    \f[ I[\mu] = \frac{1}{2} \int_\Omega \left|\tau_b + \tau_c \mathbf{V}\right|^2\,dx\,dy 
            = \frac{1}{2} \int_\Omega \left|\tau_b + \mu N \mathbf{V}\right|^2\,dx\,dy.  \f]
There is no regularization yet, but still we do not expect a minimum of zero because
the vectors \f$\tau_b\f$ and \f$\mathbf{V}\f$ are not generally parallel.

The minimization of \f$I[\mu]\f$ is presumably too strong of a goal because the
input \f$\tau_b\f$ will have significant artifacts from the original observed
surface velocity and from the application of the SSA differential operator.  Also there
are large gaps in observed surface velocity and it follows that there will be large
gaps in \f$\tau_b\f$ as well.  Furthermore the run which preceded the use of 
this inverse model involved a map of till friction angle \f$\phi\f$, and thus we may
assume there is an initial and existing field of values \f$\mu_i = \tan\phi_i\f$.

Let \f$\Omega_O\f$ be the \e open set on which observed velocities are available,
and thus \f$\tau_b\f$ is valid.  (The idea that observations are given
an open set is an artifact of the mathematical formulation below.  The point 
is that smooth functions defined on an open set can be differentiated on the
same open set.  Here this amounts to the condition invMask == 2 because of the way
invMask is constructed by readObservedSurfVels().)

In regions where there are \e no observed velocities
we want the result \f$\mu\f$ to deviate from \f$\mu_i\f$ only as necessary to 
smoothly transition to meet the goal in \f$\Omega_O\f$.  Let \f$\chi_O\f$ be the
function which is one on \f$\Omega_O\f$ and zero otherwise so \f$\chi_O\f$ 
is the ``characteristic function'' of the open set \f$\Omega_O\f$.  Assume the 
initial field \f$\mu_i\f$ is smooth.  We seek a smooth function \f$\mu\f$ which 
minimizes the regularized functional
	\f[ I_\epsilon[\mu] = \frac{1}{2} \int_{\Omega_O} \left|\tau_b
              + \mu N \mathbf{V}\right|^2\,dx\,dy 
           + \frac{1}{2} \int_{\Omega} \epsilon|\nabla(\mu - \mu_i)|^2\,dx\,dy  \f]
This functional is quadratic in \f$\mu\f$.  It is positive definite on the space 
of functions \f$\mu\f$ which have a square-integrable gradient and which go to
\f$\mu_i\f$ at the boundary (i.e. \f$\mu-\mu_i \in H_0^1(\Omega)\f$).

Let \f$\omega=\mu-\mu_i\f$.  The minimum condition for the functional \f$I_\epsilon[\mu]\f$
is a linear equation for \f$\omega\f$.  It is the Poisson-like PDE
	\f[ - \epsilon \triangle \omega + N |\mathbf{V}|^2 \chi_O\, \omega
	       = - N (\tau_b \cdot \mathbf{V} - |\mathbf{V}|^2\,\mu_i) \chi_O. \f]
We solve this PDE numerically, for \f$\omega\f$, and recover the till
friction angle by \f$\phi = \arctan\left(\omega + \mu_i\right)\f$.

To solve the PDE numerically we use a PETSc \c SNES object.  We write the equation
in slightly abstracted form as 
	\f[-\epsilon \triangle \omega + f(x,y) \omega = g(x,y).\f]
This is a \e linear equation but nonetheless we use the general, nonlinear PETSc
SNES tools because they are convenient and tunable.  An example code in the PISM
source, src/trypetsc/poisson.c, gives a standalone example of solving this kind
of Poisson-like equation.  It is a verification code for our solution method.

We call fillRegPoissonData() to fill in the coefficients \f$f(x,y)\f$ and \f$g(x,y)\f$
from known values \f$\tau_b\f$, \f$\mathbf{V}\f$, \f$N\f$ and the mask which corresponds
to \f$O\f$.  The starting till friction angle \f$\phi_i\f$ gives \f$\mu_i\f$.  Note
that the fields in IceModel which determine the effective pressure \f$N\f$, namely \c vH
and \c vHmelt, are needed to compute the coefficients; see getEffectivePressureForInverse().
Also note that \f$f(x,y)\f$ and \f$g(x,y)\f$ are identically zero outside of \f$\Omega_O\f$.
Thus \f$f(x,y)\f$ and \f$g(x,y)\f$ are generally discontinuous even if the other data
of the problem are smooth.

Regarding the size of the regularization constant \f$\epsilon>0\f$, which is input
argument \c invRegEps, we obviously expect smoother \f$\mu\f$ for larger values of
\f$\epsilon\f$.  In the functional (objective) above, the terms \f$\left|\tau_b + 
\mu N \mathbf{V}\right|^2\f$ and \f$\epsilon|\nabla(\mu - \mu_i)|^2\f$ are compared.
Since \f$|\mathbf{V}|\f$ is of size one at locations where there is significant sliding
but also some deformation, i.e. when \f$|\mathbf{V}| \sim U_{\mathtt{th}}\f$, a
small amount of regularization is a value for \f$\epsilon\f$ satisfying
	\f[ \epsilon|\nabla(\mu - \mu_i)|^2 \ll \max\left\{|\tau_b|^2,\mu^2 N^2\right\}. \f]
A scale for spatial variation of \f$\omega = \mu-\mu_i\f$ might be a change of one degree 
in \f$\phi\f$ in one kilometer, or a change of \f$\delta \omega \approx 
(\pi/180) \delta \phi = \pi/180 \approx 0.02\f$ in one kilometer,
	\f[ \frac{\delta \omega}{\delta x} = 2\times 10^{-5}\, \text{m}^{-1}. \f]
If \f$10^5\f$ Pa = 1 bar is a scale for \f$\tau_b\f$ then, by the above standard, small
regularization might be
        \f[ \epsilon|\nabla \omega|^2 \ll 10^{10}\, \text{Pa}^2, \f]
Using the above estimate of \f$\delta \omega/\delta x\f$ as a size for the
gradient of \f$\omega\f$, we get
	\f[ \epsilon \ll  \frac{10^{10}}{4 \times 10^{-10}} 
	                  \approx 10^{19}\,\text{Pa}^2\,\text{m}^2. \f]

It may also be fair to relate values for \f$\epsilon\f$ to the steady state
of a notional diffusion (smoothing) process for \f$\omega = \mu-\mu_i\f$, but this idea
has not been pursued, yet.

The above is the whole story if the ice is grounded and there is positive
ice thickness.  If a point is marked \c MASK_FLOATING or \c MASK_FLOATING_OCEAN0, 
however, or if a point is grounded but the ice thickness is zero, then 
\f$\mu=\mu_i\f$.

Finally we enforce the limits \c phi_min, \c phi_max.
 */
PetscErrorCode IceModel::computeTFAFromBasalShear(
                const PetscScalar phi_low, const PetscScalar phi_high,
                const PetscScalar invRegEps, const char *invfieldsfilename) {               
  PetscErrorCode ierr;
  
  // choose regularization; see guesses in comment at top
  RegPoissonCtx  user;        // user-defined work context; see iceModel.hh
  user.epsilon = invRegEps; 
  user.da = grid.da2;
  user.grid = &grid;

  // space for f(x,y) and g(x,y)
  user.f = new IceModelVec2S;
  ierr = user.f->create(grid, "f_invmodel", false); CHKERRQ(ierr);
  user.g = new IceModelVec2S;
  ierr = user.g->create(grid, "g_invmodel", false); CHKERRQ(ierr);

  // main added content relative to ex5: fill nontrivial coeffs f,g for reg version
  ierr = fillRegPoissonData(user); CHKERRQ(ierr);

  // write f,g to file, if a filename is given
  if (strlen(invfieldsfilename) > 0) {
    ierr = user.f->set_attrs(
       "inverse_output", "f(x,y), coeff in Poisson-like eqn for regularizing mu",
       "", ""); CHKERRQ(ierr);
    ierr = user.f->write(invfieldsfilename, NC_FLOAT); CHKERRQ(ierr);
    ierr = user.g->set_attrs(
       "inverse_output", "g(x,y), coeff in Poisson-like eqn for regularizing mu",
       "", ""); CHKERRQ(ierr);
    ierr = user.g->write(invfieldsfilename, NC_FLOAT); CHKERRQ(ierr);
  }

  // create appropriate SNES, Mat
  SNES               snes;        /* nonlinear solver */
  Vec                x,r;         /* solution, residual vectors */
  PetscInt           its;         /* iterations for convergence */

  ierr = SNESCreate(grid.com,&snes);CHKERRQ(ierr);
  ierr = DACreateGlobalVector(user.da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
  
  // see src/trypetsc/poisson.c for ideas here
  ierr = SNESSetFunction(snes,r,SNESDAFormFunction,&user);CHKERRQ(ierr);

  ierr = DASetLocalFunction(user.da,(DALocalFunction1)RegPoissonFunctionLocal);CHKERRQ(ierr);
  ierr = PetscOptionsSetValue("-snes_mf", PETSC_NULL); CHKERRQ(ierr);  // FIXME?: no Jacobian!!

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  ierr = VecSet(x,0.0); CHKERRQ(ierr); // initial guess \omega = 0 is \mu = \mu_i

  // solve the linear problem as though it is nonlinear
  ierr = SNESSolve(snes, PETSC_NULL, x);CHKERRQ(ierr);

  // some feedback
  PetscReal resnorm;
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr); 
  ierr = VecNorm(r,NORM_INFINITY,&resnorm); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
            "    done; number of Newton iterations = %d;  |residual|_infty = %9.3e\n",
            its, resnorm); CHKERRQ(ierr);

  // convert solution in x = \omega = \mu - \mu_i into phi 
  PetscScalar **tillphi, **oldphi, **result;
  ierr = vtillphi.get_array(tillphi);  CHKERRQ(ierr);
  ierr = inv.oldtillphi->get_array(oldphi); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, x, &result); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar muinitial = tan((pi/180.0) * oldphi[i][j]);
      tillphi[i][j] = (180.0/pi) * atan(result[i][j] + muinitial);
      if (tillphi[i][j] > phi_high)  tillphi[i][j] = phi_high;
      if (tillphi[i][j] < phi_low)   tillphi[i][j] = phi_low;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, x, &result); CHKERRQ(ierr);
  ierr = vtillphi.end_access();  CHKERRQ(ierr);
  ierr = inv.oldtillphi->end_access(); CHKERRQ(ierr);

  // de-allocate SNES stuff
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);      
  delete user.f;  user.f = PETSC_NULL;
  delete user.g;  user.g = PETSC_NULL;
  ierr = SNESDestroy(snes);CHKERRQ(ierr);

  return 0;
}


//! Compute non-constant coefficients in Poisson-like regularized problem for yield stress in inverse modeling.
PetscErrorCode IceModel::fillRegPoissonData(RegPoissonCtx &user) {
  PetscErrorCode ierr;

  PetscScalar **oldphi, **imask, **N, **taubx, **tauby,
              **ff, **gg;
  PetscScalar magVsqr, V_x, V_y;

  ierr = inv.oldtillphi->get_array(oldphi);  CHKERRQ(ierr);
  ierr = inv.invMask->get_array(imask);  CHKERRQ(ierr);
  ierr = inv.effPressureN->get_array(N);  CHKERRQ(ierr);
  ierr = inv.taubxComputed->get_array(taubx); CHKERRQ(ierr);
  ierr = inv.taubyComputed->get_array(tauby); CHKERRQ(ierr);
  ierr = ssavel.begin_access(); CHKERRQ(ierr);
  ierr = user.f->get_array(ff); CHKERRQ(ierr);
  ierr = user.g->get_array(gg); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (imask[i][j] < 1.5) { // coeffs f,g have factor \chi_O; stencil width needed to be in O
        ff[i][j] = 0.0;
        gg[i][j] = 0.0;
      } else {
        ierr = getVfromUforInverse(ssavel(i,j).u,
				   ssavel(i,j).v, V_x, V_y, magVsqr); CHKERRQ(ierr);
        const PetscScalar 
              taubdotV = taubx[i][j] * V_x + tauby[i][j] * V_y,
              muinitial = tan((pi/180.0) * oldphi[i][j]);
        // f(x,y) = + N |V|^2 \chi_O
        ff[i][j] = N[i][j] * magVsqr;
        // g(x,y) = - N (tau_b . V - |V|^2 \mu_i) \chi_O
        gg[i][j] = - N[i][j] * (taubdotV - magVsqr * muinitial);
      }
    }
  }
  ierr = user.f->end_access(); CHKERRQ(ierr);
  ierr = user.g->end_access(); CHKERRQ(ierr);
  ierr = inv.oldtillphi->end_access();  CHKERRQ(ierr);
  ierr = inv.invMask->end_access();  CHKERRQ(ierr);
  ierr = inv.effPressureN->end_access();  CHKERRQ(ierr);
  ierr = inv.taubxComputed->end_access(); CHKERRQ(ierr);
  ierr = inv.taubyComputed->end_access(); CHKERRQ(ierr);
  ierr = ssavel.end_access(); CHKERRQ(ierr);

  return 0; 
}


//! Evaluates function F(x) for SNES for a Poisson-like regularized problem in inverse modeling.
/*!
This method is used in IceModel::computeTFAFromBasalShear().  It is handed to the SNES.

This is one of several routines called by IceModel::invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

Compare examples/trypetsc/poisson.c.
 */
PetscErrorCode RegPoissonFunctionLocal(
                  DALocalInfo *info, PetscScalar **x, PetscScalar **F,
                  RegPoissonCtx *user) {
  PetscErrorCode ierr;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      dx,dy,sc,scxx,scyy;
  PetscScalar    u,neguxx,neguyy;
  PetscScalar    **ff, **gg;

  PetscFunctionBegin;

  // use transpose as in IceModel
  Mx = info->my; My = info->mx;
  xs = info->ys; ys = info->xs;
  xm = info->ym; ym = info->xm;

  // scaling constants
  dx   = (user->grid)->dx;
  dy   = (user->grid)->dy;
  sc   = dx * dy;
  scxx = sc / (dx * dx);
  scyy = sc / (dy * dy);

  ierr = user->f->get_array(ff); CHKERRQ(ierr);
  ierr = user->g->get_array(gg); CHKERRQ(ierr);
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        F[i][j] = x[i][j];
      } else {
        u       = x[i][j];
        neguxx  = (2.0*u - x[i-1][j] - x[i+1][j]) * scxx;
        neguyy  = (2.0*u - x[i][j-1] - x[i][j+1]) * scyy;
        F[i][j] = user->epsilon * (neguxx + neguyy) + sc * ff[i][j] * u - sc * gg[i][j];
      }
    }
  }
  ierr = user->f->end_access(); CHKERRQ(ierr);
  ierr = user->g->end_access(); CHKERRQ(ierr);

  PetscFunctionReturn(0); 
} 

