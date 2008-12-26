// Copyright (C) 2008 Ed Bueler
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


/* THIS SOURCE FILE CONTAINS THE TWO MOST NONTRIVIAL ROUTINES IN THE INVERSE
   MODEL COMPONENTS OF IceModel, NAMELY THE TWO WHICH BUILD Mat and SNES OBJECTS  */

//! Compute basal shear stress from a given sliding velocity using SSA equations.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

We assume the basal sliding velocity is known.  Using ice thickness, 
ice surface elevation, and an ice temperature field, we use 
SSA equations to compute the basal shear stress.

There is no smoothing in this initial implementation.

In particular:
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
\f$\tau_b = (\tau_{(b)x},\tau_{(b)y})\f$.  The equations are solved everywhere 
that there is positive thickness; at other points this procedure returns 
\f$\tau_b = 0\f$.

This procedure calls

- assembleSSAMatrix(),

- assembleSSARhs(),

- computeEffectiveViscosity().

It saves temporary versions of various flags modified by those routines.

The result in IceModelVec2 taubx_out, tauby_out must be allocated (created) 
before calling this routine.
 */
PetscErrorCode IceModel::computeBasalShearFromSSA(
                          IceModelVec2 ub_in, IceModelVec2 vb_in,
			  IceModelVec2 &taubx_out, IceModelVec2 &tauby_out) {
  PetscErrorCode ierr;

  // effective viscosity for ub_in, vb_in: save flags before changing them
  const PetscTruth leaveNuHAloneSSA_save = leaveNuHAloneSSA, 
                   useConstantNuHForSSA_save = useConstantNuHForSSA;
  leaveNuHAloneSSA = PETSC_FALSE;
  useConstantNuHForSSA = PETSC_FALSE;
  IceModelVec2 myvNuH[2] = {vWork2d[0], vWork2d[1]}; // already allocated space
  IceModelVec2 vubarSSAold = vWork2d[2], 
               vvbarSSAold = vWork2d[3];
  ierr = vubarSSA.copy_to(vubarSSAold); CHKERRQ(ierr);
  ierr = vvbarSSA.copy_to(vvbarSSAold); CHKERRQ(ierr);
  ierr = ub_in.copy_to(vubarSSA); CHKERRQ(ierr);
  ierr = vb_in.copy_to(vvbarSSA); CHKERRQ(ierr);
  // eps=0.0 in bdd-below regularization; Schoof type regularization does occur;
  ierr = computeEffectiveViscosity(myvNuH, 0.0); CHKERRQ(ierr);
  ierr = vubarSSAold.copy_to(vubarSSA); CHKERRQ(ierr);
  ierr = vvbarSSAold.copy_to(vvbarSSA); CHKERRQ(ierr);
  leaveNuHAloneSSA = leaveNuHAloneSSA_save;
  useConstantNuHForSSA = useConstantNuHForSSA_save;

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
  ierr = assembleSSARhs(false, rhs); CHKERRQ(ierr);

  // Set x = [u v]'  (interleaved).
  PetscScalar **ub, **vb;
  const PetscInt  twoMy = 2 * grid.My;
  ierr = VecSet(x, 0.0); CHKERRQ(ierr);
  ierr = ub_in.get_array(ub); CHKERRQ(ierr);
  ierr = vb_in.get_array(vb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = VecSetValue(x, i*twoMy + 2*j, ub[i][j], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(x, i*twoMy + 2*j+1, vb[i][j], INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = ub_in.end_access(); CHKERRQ(ierr);
  ierr = vb_in.end_access(); CHKERRQ(ierr);

  // assemble matrix and vec before multiplication; communicate!
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

  // With new A and x, compute  result = Ax - rhs = L[u] - (driving)
  // HERE IS INVERSE MODEL!!
  
  ierr = MatMult(A,x,result); CHKERRQ(ierr);  
  
  ierr = VecAXPY(result,-1.0,rhs); CHKERRQ(ierr);  // note result = 0 if MASK_SHEET
  
  // transfer result to taub output vectors; compare code in
  //   IceModel::moveVelocityToDAVectors()
  PetscScalar     **tbx, **tby, *res;
  Vec             resultLoc = SSAXLocal;
  ierr = VecScatterBegin(SSAScatterGlobalToLocal, result, resultLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(SSAScatterGlobalToLocal, result, resultLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGetArray(resultLoc, &res); CHKERRQ(ierr);
  ierr = taubx_out.get_array(tbx); CHKERRQ(ierr);
  ierr = tauby_out.get_array(tby); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      tbx[i][j] = res[i*twoMy + 2*j];
      tby[i][j] = res[i*twoMy + 2*j+1];
    }
  }
  ierr = taubx_out.end_access(); CHKERRQ(ierr);
  ierr = tauby_out.end_access(); CHKERRQ(ierr);
  ierr = VecRestoreArray(resultLoc, &res); CHKERRQ(ierr);

  // de-allocate
  ierr = MatDestroy(A); CHKERRQ(ierr);
  ierr = VecDestroy(x); CHKERRQ(ierr);
  ierr = VecDestroy(result); CHKERRQ(ierr);
  ierr = VecDestroy(rhs); CHKERRQ(ierr);

  return 0;
}


//! Compute the yield stress for the given basal shear stress using a regularized objective.
/*!
The more complete brief description might be: "Compute till yield stress
from saved basal shear stress map and saved basal velocity 
(both vector-valued) using pseudo-plastic till constitutive relation, and so
that the resulting yield stress map minimizes a combination of a smoothness
measure and the difference between yield stress in the pseudo-plastic model
and the basal shear stress."

This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

This procedure should not be called in the purely plastic case, that is, 
if <tt>PlasticBasalType::pseudo_plastic</tt> is \c FALSE.  Recall that 
when <tt>pseudo_plastic</tt> is \c TRUE, the method PlasticBasalType::drag() 
will compute the basal shear stress as
    \f[ \tau_{(b)}
           = - \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U}    \f]
where \f$\tau_{(b)}=(\tau_{(b)x},\tau_{(b)y})\f$, \f$U=(u,v)\f$,
\f$q=\f$ <tt>pseudo_q</tt>, and 
\f$U_{\mathtt{th}}=\f$ <tt>pseudo_u_threshold</tt>.
Typical values for the constants are \f$q=0.25\f$ and \f$U_{\mathtt{th}} = 100\f$
m/a.  The linearly viscous till case \f$q=1\f$ is allowed, in which case 
\f$\beta = \tau_c/U_{\mathtt{th}}\f$.

Let
	\f[ \mathbf{V} = \frac{\mathbf{U}}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q}, \f]
so that, in this procedure, \f$\mathbf{V}\f$ is a known vector field (in the map-plane).

We want to determine \f$\tau_c\f$ given \f$\tau_{(b)}\f$ and \f$\mathbf{U}\f$ in the
above equation.  The two latter quantities are vectors, however, and they may not be
parallel.  Thus we refine the goal to the minimization of the functional
	\f[ I[\tau_c] = \int_\Omega \left|\tau_{(b)}
           + \tau_c \mathbf{V}\right|^2\,dx\,dy  \f]
and we do not expect a minimum of zero.

Even the minimization of \f$I[\tau_c\f$ is probably too strong of a goal because the
input \f$\tau_{(b)}\f$ will have significant artifacts from the original observed
surface velocity and from the application of the SSA differential operator.
Thus we seek a smoother function \f$\tau_c\f$ which minimizes the regularized functional
	\f[ I_\epsilon[\tau_c] = \int_\Omega \left|\tau_{(b)}
           + \tau_c \mathbf{V}\right|^2 + \epsilon|\nabla\tau_c|^2\,dx\,dy  \f]

This functional is quadratic in \f$\tau_c\f$ so its minimum condition is a linear
equation for \f$\tau_c\f$.  It is this Poisson-like problem
	\f[ - \epsilon \triangle \tau_c + |\mathbf{V}|^2 \tau_c
	             = -\tau_{(b)} \cdot \mathbf{V}. \f]
Here we have assumed, and use as a boundary value when solving, that \f$\tau_c = 0\f$
at the edge of the computational domain, which is compatible with the whole ice sheet
nature of PISM, and with the mathematical analysis in \lo\cite{SchoofStream}\elo.
That is, we assume the ice sheet is surrounded by ice shelves or open ocean.

So we use a PETSc \c SNES object to set up and solve this Poisson-like equation.
We write the equation in slightly abstracted form as 
\f$-\epsilon \triangle \tau_c + f(x,y) \tau_c = g(x,y)\f$.
This is a \e linear equation but nonetheless we just use the general PETSc \c SNES
nonlinear methods.  We call fillRegPoissonTaucData() to fill in the
coefficients \f$f(x,y)\f$ and \f$g(x,y)\f$ from known values \f$\tau_{(b)}\f$ 
and \f$\mathbf{U}\f$ and parameters <tt>pseudo_q</tt>, <tt>pseudo_u_threshold</tt>
in PlasticBasalType.

An example code in the PISM source, examples/trypetsc/poisson.c, gives
a standalone example of solving this kind of Poisson-like equation.  It also
verifies the method we use here.

Regarding the size of the regularization constant \f$\epsilon\f$, we expect a
smoother result \f$\tau_c\f$ for larger (positive) values of \f$\epsilon\f$.
In the functional (objective) above, the terms \f$\left|\tau_{(b)} + 
\tau_c \mathbf{V}\right|^2\f$ and \f$\epsilon|\nabla\tau_c|^2\f$ are compared.
We assert \f$|\mathbf{V}|\f$ is of size one, so we assert that a reasonable 
value for \f$\epsilon\f$ follows from the heuristic
	\f[ \epsilon|\nabla\tau_c|^2 \ll \max\left\{|\tau_{(b)}|^2,\tau_c^2\right\}. \f]
We take this value as the scale for spatial variation in \f$\tau_c\f$:
	\f[ \frac{\Delta \tau_c}{\Delta x} = \frac{10^5\,\text{Pa}}{10^4\,\text{m}} 
	            = 10 \frac{\text{Pa}}{\text{m}}. \f]
Supposing we want
        \f[ \epsilon|\nabla\tau_c|^2 = (10^4\, \text{Pa})^2 = 10^8 \,\text{Pa}^2, \f]
we conclude with a value:
	\f[ \epsilon = 10^6\,\text{m}^2. \f]
This seems to be comparable to the area of a grid cell, i.e. related to the steady state
of a notional diffusion (smoothing) process for \f$\tau_c\f$.  It may be that 
\f$\epsilon = 10^8\,\text{m}^2\f$ is a good value to start with.

The above is the whole story if the ice is grounded and there is positive
ice thickness.  If a point is marked \c MASK_FLOATING or \c MASK_FLOATING_OCEAN0, 
however, then the yield stress is set to zero.  If a point is grounded but the ice
thickness is zero then the yield stress is set to a high value of
1000 kPa = 10 bar.

This procedure reads N??? fields from the current model state, namely ice 
thickness \c vH, the mask \c vMask, ???
 */
PetscErrorCode IceModel::computeYieldStressFromBasalShearUsingPseudoPlastic(
                IceModelVec2 ub_in, IceModelVec2 vb_in,
	        IceModelVec2 taubx_in, IceModelVec2 tauby_in, 
                IceModelVec2 &tauc_out) {               
  PetscErrorCode ierr;
  
  if (doPseudoPlasticTill == PETSC_FALSE) {
    ierr = verbPrintf(1, grid.com, 
       "WARNING: computeTFAFromBasalShearStress() should only be called with q > 0.0\n"
       "  in pseudo-plastic model;  here is PlasticBasalType::printInfo() output:\n");
       CHKERRQ(ierr);
    ierr = basal->printInfo(1,grid.com); CHKERRQ(ierr);
    ierr = verbPrintf(1, grid.com, "  CONTINUING.  May crash.\n"); CHKERRQ(ierr);
  }

  RegPoissonTaucCtx  user;        // user-defined work context; see iceModel.hh

  // choose regularization; see guesses in comment at top
  user.epsilon = 1.0e8; 
  ierr = PetscOptionsGetReal(PETSC_NULL, "-inverse_reg_eps", &user.epsilon, PETSC_NULL);
             CHKERRQ(ierr);

  user.da = grid.da2;
  user.grid = &grid;

  SNES               snes;        /* nonlinear solver */
  Vec                x,r;         /* solution, residual vectors */
  Mat                J;           /* Jacobian matrix */
  PetscInt           its;         /* iterations for convergence */

  ierr = SNESCreate(grid.com,&snes);CHKERRQ(ierr);

  ierr = DACreateGlobalVector(user.da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.f);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.g);CHKERRQ(ierr);

  // main added content relative to ex5: fill nontrivial coeffs
  ierr = fillRegPoissonTaucData(ub_in,vb_in,taubx_in,tauby_in,user); CHKERRQ(ierr);
  
  ierr = DAGetMatrix(user.da,MATAIJ,&J);CHKERRQ(ierr);
  
  // use default method of Jacobian eval (use SNESDAComputeJacobian option and
  //   DASetLocalJacobian() below); also preconditioner is same as Jacobian;
  //   compare different approaches here in ex5.c
  ierr = SNESSetJacobian(snes,J,J,SNESDAComputeJacobian,&user);CHKERRQ(ierr); // default

  ierr = DASetLocalFunction(user.da,(DALocalFunction1)RegPoissonTaucFunctionLocal);CHKERRQ(ierr);
  ierr = DASetLocalJacobian(user.da,(DALocalFunction1)RegPoissonTaucJacobianLocal);CHKERRQ(ierr); 

  ierr = SNESSetFunction(snes,r,SNESDAFormFunction,&user);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  ierr = VecSet(x,0.0); CHKERRQ(ierr); // a linear problem, so initial guess is easy ...

  // solve the linear problem as though it is nonlinear!
  ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr);

  // some feedback appropriate, but FIXME
  PetscReal resnorm;
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr); 
  ierr = VecNorm(r,NORM_INFINITY,&resnorm); CHKERRQ(ierr);
  ierr = PetscPrintf(grid.com,
            "number of Newton iterations = %d;  |residual|_infty = %9.3e\n",
            its, resnorm); CHKERRQ(ierr);

  // copy solution in x into tauc_out, zeroing the yield stress where floating 
  PetscScalar **tauc, **mask, **result;
  ierr = tauc_out.get_array(tauc);  CHKERRQ(ierr);
  ierr =    vMask.get_array(mask);  CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, x, &result); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        tauc[i][j] = 0.0;
      } else {
        tauc[i][j] = result[i][j];
      }
    }
  }
  ierr = tauc_out.end_access(); CHKERRQ(ierr);
  ierr =    vMask.end_access(); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, x, &result); CHKERRQ(ierr);

  // de-allocate SNES stuff
  ierr = MatDestroy(J);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);      
  ierr = VecDestroy(user.f);CHKERRQ(ierr);      
  ierr = VecDestroy(user.g);CHKERRQ(ierr);      
  ierr = SNESDestroy(snes);CHKERRQ(ierr);

  return 0;
}


//! Compute non-constant coefficients in Poisson-like regularized problem for yield stress in inverse modeling.
PetscErrorCode IceModel::fillRegPoissonTaucData(
                    IceModelVec2 ub_in, IceModelVec2 vb_in,
	            IceModelVec2 taubx_in, IceModelVec2 tauby_in,
	            RegPoissonTaucCtx &user) {
  PetscErrorCode ierr;

  PetscScalar **ub, **vb, **taubx, **tauby, **ff, **gg;
  const PetscReal q   = basal->pseudo_q,
                  Uth = basal->pseudo_u_threshold;

  ierr = taubx_in.get_array(taubx); CHKERRQ(ierr);
  ierr = tauby_in.get_array(tauby); CHKERRQ(ierr);
  ierr = ub_in.get_array(ub); CHKERRQ(ierr);
  ierr = vb_in.get_array(vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(user.da, user.f, &ff); CHKERRQ(ierr);
  ierr = DAVecGetArray(user.da, user.g, &gg); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // compute                  U
      //              V = ------------------
      //                  |U|^{1-q} U_th^q
      const PetscScalar 
              U_x = ub[i][j],
              U_y = vb[i][j],
              denom = pow(PetscSqr(U_x)+PetscSqr(U_y),(1.0-q)/2.0) * pow(Uth,q),
              V_x = U_x / denom,
              V_y = U_y / denom;
      // and  f(x,y) = + |V|^2  and  g(x,y) = - tau_b . V  (dot product)
      ff[i][j] = PetscSqr(V_x) + PetscSqr(V_y);
      gg[i][j] = - taubx[i][j] * V_x - tauby[i][j] * V_y;
    }
  }
  ierr = DAVecRestoreArray(user.da, user.f, &ff); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user.da, user.g, &gg); CHKERRQ(ierr);
  ierr = ub_in.end_access(); CHKERRQ(ierr);
  ierr = vb_in.end_access(); CHKERRQ(ierr);
  ierr = taubx_in.end_access(); CHKERRQ(ierr);
  ierr = tauby_in.end_access(); CHKERRQ(ierr);

  return 0; 
}


//! Evaluates function F(x) for SNES, corresponding to Poisson-like regularized problem for yield stress in inverse modeling.
/*!
This method is used in computeYieldStressFromBasalShearUsingPseudoPlastic().  It is handed to
the SNES.

Compare examples/trypetsc/poisson.c.
 */
PetscErrorCode RegPoissonTaucFunctionLocal(
                  DALocalInfo *info, PetscScalar **x, PetscScalar **F,
                  RegPoissonTaucCtx *user) {
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

  ierr = DAVecGetArray(user->da, user->f, &ff); CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da, user->g, &gg); CHKERRQ(ierr);
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
  ierr = DAVecRestoreArray(user->da, user->f, &ff); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da, user->g, &gg); CHKERRQ(ierr);

  PetscFunctionReturn(0); 
} 


//! Evaluates Jacobian J=dF(x) for SNES, corresponding to Poisson-like regularized problem for yield stress in inverse modeling.
/*!
This method is used in computeYieldStressFromBasalShearUsingPseudoPlastic().  It is handed to
the SNES.

Compare examples/trypetsc/poisson.c.
 */
PetscErrorCode RegPoissonTaucJacobianLocal(
                  DALocalInfo *info, PetscScalar **x, Mat jac,
                  RegPoissonTaucCtx *user) {
  PetscErrorCode ierr;
  PetscInt       i,j;
  MatStencil     col[5],row;
  PetscScalar    v[5],dx,dy,sc,scxx,scyy,eps;
  PetscScalar    **ff;

  PetscFunctionBegin;

  dx   = (user->grid)->dx;
  dy   = (user->grid)->dy;
  sc   = dx * dy;
  scxx = sc/(dx * dx);
  scyy = sc/(dy * dy);
  eps  = user->epsilon;

  ierr = DAVecGetArray(user->da, user->f, &ff); CHKERRQ(ierr);
  for (i=info->xs; i<info->xs+info->xm; i++) {
    for (j=info->ys; j<info->ys+info->ym; j++) {
      row.j = j; row.i = i;
      /* boundary points */
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) {
        v[0] = 1.0;
        ierr = MatSetValuesStencil(jac,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
      } else {
      /* interior grid points */
        v[0] = -scyy;                                     col[0].j = j - 1; col[0].i = i;
        v[1] = -scxx;                                     col[1].j = j;     col[1].i = i-1;
        v[2] = eps * 2.0 * (scxx + scyy) + sc * ff[i][j]; col[2].j = row.j; col[2].i = row.i;
        v[3] = -scxx;                                     col[3].j = j;     col[3].i = i+1;
        v[4] = -scyy;                                     col[4].j = j + 1; col[4].i = i;
        ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = DAVecRestoreArray(user->da, user->f, &ff); CHKERRQ(ierr);

  /* Assemble matrix, using the 2-step process */
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  /* Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error. */
#if (PISM_HAVE_PETSC3)
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
#else
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

