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

#include <cstring>
#include <cstdlib>
#include <petscda.h>
#include <netcdf.h>
#include "nc_util.hh"
#include "iceModel.hh"

/* TO DO:
-- read user options "-cbar_for_basal foo.nc" and "-csurf_for_basal foo.nc"
59.  write out cbase as diagnostic, routinely
60.  write out till friction angle as diagnostic, routinely
61.  write out magnitude of basal shear stress as diagnostic, routinely

*/

/*
QUOTING FROM EMAIL TO TOLLY:

I am about to (re-)implement a scheme for reading a map of vertically
averaged horizontal velocity from a NetCDF file and using it, along with
the velocity field produced by the thermomechanically-coupled non-sliding
SIA, and with the MacAyeal equations, to determine a basal shear stress
and a basal sliding velocity.

That is, the shear from the SIA will be removed from the vertically
averaged velocity, the rest will be attributed to sliding, and this
sliding velocity field will be fed into the MacAyeal equations to estimate
a basal shear stress.  This is inverse modeling.  Done badly, compared to
the principled methods underway by Martin and Hilmar and Olga and Doug and
Ian etc.  This is the part where I am hoping interaction with Martin et al
will give better methods, but ones which are still computable in a
PISM-like framework at full ice sheet scale.

Spin up matters in all this because the removal of vertical shear using
the SIA depends fairly critically on the flow factor in the flow law, and
thus on the temperature field within the ice.

The "vertically averaged horizontal velocity" could be produced by any
other model, but I have in mind using the balance velocities from Bamber. 
He produced them from a lousy model for ice sheets, namely the underlying
idea for balance velocities(!), but at least that model produces high
velocities where they are justified by large catchments and substantial
surface accumulation.  In other words, I know how to, and am about to
(re)write the PISM code necessary to read in balance velocities and use
them to estimate sliding ``coefficients''.  This code was written in draft
for Antarctica runs, but I can make it part of standard PISM now.

If we have surface velocities, the modifications are simple.  I will
subtract SIA-predicted vertical shear from the surface velocities, not the
vertically-averaged velocities, to produce basal sliding, and then proceed
as before.

There will be new PISM options like "-cbar_for_basal foo.nc" and
"-csurf_for_basal foo.nc"; "c" means ice velocity (actually, speed) in a
number of PISM options.


*/

/*!
This procedure computes a z-independent horizontal ice velocity, conceptually
a basal sliding velocity, from either surface horizontal ice speed or 
a vertically-averaged horizontal speed by using the thermomechanically coupled 
SIA equations.

We could be given a surface horizontal speed \c csurf or a vertically-integrated 
horizontal speed \c cbar.  The argument \c CisSURF is \c PETSC_TRUE in the 
former case and \c PETSC_FALSE in the latter case.

It is assumed that the input speed gives a velocity pointing down the 
surface gradient.  That is,
  \f[ (u,v) = - c \frac{\nabla h}{|\nabla h|} \f]
at all points where \f$|\nabla h| \ge \mathtt{alphacrit}\f$.  At points where 
\f$|\nabla h| < \mathtt{alphacrit}\f$ the returned velocity is zero.  

Where there is no ice the returned velocity is zero.

This procedure must save temporary versions of everything modified by 
IceModel::velocitySIAStaggered() because it calls that routine to compute 
\c vuvbar[2].

The result in \c Vecs \c u and \c v is on the regular grid.
 */
PetscErrorCode IceModel::removeVerticalPlaneShearRateFromC(
                 const PetscTruth CisSURF, Vec myC, const PetscScalar alphacrit,
                 Vec &u, Vec &v) {
  PetscErrorCode ierr;
  SETERRQ(1,"NOT IMPLEMENTED");
  return 0;
}


//! Compute basal shear stress from a given z-independent horizontal ice velocity using SSA equations.
/*!
From the input velocity, and using ice thickness, ice surface elevation, and an 
ice temperature field, we use the SSA equations to compute the basal shear
stress.  In particular, the equations we are using are
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

In the case where the known field is \f$(u,v)\f$.  We are solving for 
\f$\tau_b = (\tau_{(b)x},\tau_{(b)y})\f$.

The equations are solved everywhere that there is positive thickness; at other points 
this procedure returns \f$\tau_b = 0\f$.

This procedure must save temporary versions of everything modified by 
IceModel::assembleSSAMatrix() and IceModel::assembleSSARhs() and
IceModel::computeEffectiveViscosity()
because it calls those routines.

The result in \c Vecs \c taubx and \c taubx are on the regular grid.
 */
PetscErrorCode IceModel::computeBasalShearFromSSA(
                 Vec myu, Vec myv, Vec &taubx, Vec &tauby) {
  PetscErrorCode ierr;
  SETERRQ(1,"NOT IMPLEMENTED");

#if 0
  Mat A;
  Vec x, result, rhs;
  PetscScalar **ub, **vb, **cbarin, **taub;

//FIXME: need to compute effective viscosity first

  // allocate Mat and Vecs; compare IceModel::velocitySSA()
  const PetscInt M = 2 * grid.Mx * grid.My;
  ierr = MatDuplicate(SSAStiffnessMatrix, &A); CHKERRQ(ierr);
//  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE, M, M,
//                         13, PETSC_NULL, 13, PETSC_NULL,
//                         &A); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &result); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &rhs); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &x); CHKERRQ(ierr);

  // build approximation of SSA equations (A x = rhs)
  // FIXME:  Issue here is about how assembleSSAMatrix() treats the mask.
  //         Generally we want basal shear as part of matrix:
  //                L[U] + beta U = (driving),
  //         so  A = L + beta I.
  //         In this case we don't want basal shear terms as part of the matrix:
  //                tau_b := L[U] - (driving),
  //         so  A = L.
// I think we should save vMask and rebuild it based only on H>0; if H>0 mark
// as FLOATING; if H==0 mark as SHEET.
  Vec vNu[2] = {vWork2d[0], vWork2d[1]}; // already allocated space
  ierr = assembleSSAMatrix(vNu, A); CHKERRQ(ierr);
  ierr = assembleSSARhs(false, rhs); CHKERRQ(ierr);

  // Note rhs contains driving terms  \rho g H \grad h  but  A  contains basal
  // drag term [betax*u betay*v]'.  It must be removed.
  // Also set x = [u v]'  (i.e. interleaved).
  ierr = VecSet(x, 0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vmagbalvel, &balvel); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt    J = 2*j;
      const PetscInt    rowU = i*2*grid.p->My + J;
      const PetscInt    rowV = i*2*grid.p->My + J+1;
      // remove old drag term
      if (intMask(mask[i][j]) == MASK_DRAGGING) {
        ierr = MatSetValue(A, rowU, rowU,
                           - basalDragx(beta, tauc, u, v, i, j), ADD_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(A, rowV, rowV,
                           - basalDragy(beta, tauc, u, v, i, j), ADD_VALUES); CHKERRQ(ierr);
      }
      // remove deformational from balance velocities; put in x
      const PetscScalar c = sqrt(PetscSqr(u[i][j]) + PetscSqr(v[i][j]));
      PetscScalar       residualFactor = (balvel[i][j] / c) - 1.0;
      // if (residualFactor < 0.0)   residualFactor = 0.0;  // avoid negative drag coeff here?  probably not
      ierr = VecSetValue(x, rowU, residualFactor * u[i][j], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(x, rowV, residualFactor * v[i][j], INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbeta, &beta); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vmagbalvel, &balvel); CHKERRQ(ierr);
  // communicate!
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

  // With new A and x, compute  result = Ax - rhs.
  ierr = MatMult(A,x,result); CHKERRQ(ierr);
  ierr = VecAXPY(result,-1.0,rhs); CHKERRQ(ierr);  // note result = 0 if MASK_SHEET
  
  // As result is [betax*u betay*v]', divide by velocities, but only where grounded!
  // where floating, report no drag coeff
  ierr = VecPointwiseDivide(dragxy,result,x); CHKERRQ(ierr);
  ierr = VecScale(dragxy,-1.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt    J = 2*j;
      const PetscInt    rowU = i*2*grid.p->My + J;
      const PetscInt    rowV = i*2*grid.p->My + J+1;
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        ierr = VecSetValue(dragxy, rowU, 0.0, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(dragxy, rowV, 0.0, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dragxy); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dragxy); CHKERRQ(ierr);
  
  ierr = MatDestroy(A); CHKERRQ(ierr);
  ierr = VecDestroy(x); CHKERRQ(ierr);
  ierr = VecDestroy(result); CHKERRQ(ierr);
  ierr = VecDestroy(rhs); CHKERRQ(ierr);

#endif
  return 0;
}


/*!
First, here is the context for this procedure:

This procedure should not be called if <tt>PlasticBasalType::pseudo_plastic</tt>
is \c FALSE, that is, in the purely plastic case.  (This case is validly
implemented in \c PlasticBasalType, however.)

If <tt>PlasticBasalType::pseudo_plastic</tt> is \c TRUE then
<tt>PlasticBasalType::drag()</tt> will compute the basal shear stress
as
    \f[ (\tau_{(b)x},\tau_{(b)y})
           = - \frac{\tau_c}{|U|^{1-q} |U_{\mathtt{thresh}}|^q} U\f]
where \f$U=(u,v)\f$, \f$q=\f$ <tt>pseudo_q</tt>, and
\f$U_{\mathtt{thresh}}=\f$ <tt>pseudo_u_threshold</tt>.
Typical values for the latter two constants are \f$q=0.25\f$ 
and 100 m/a for \f$U_{\mathtt{thresh}}\f$, but \f$q\f$ could even be
set to 1 for a linear till model, in which case
\f$\beta = \tau_c/|U_{\mathtt{thresh}}|\f$.

On the other hand, <tt>updateYieldStressFromHmelt()</tt> computes
the till yield stress \f$\tau_c\f$, a scalar, by 
    \f[ \tau_c = c_0 + \tan\left(\frac{\pi}{180} \phi\right)\, N \f]
where \f$\phi\f$ is the map of till friction angle, in degrees, stored in
\c vtillphi.  The till cohesion \f$c_0\f$ (= \c plastic_till_c_0) is usually zero.
Here \f$N\f$ is the effective pressure on the till, namely \f$N = \rho g H - p_w\f$.
The porewater pressure \f$p_w\f$ is estimated in terms of the stored basal water
\c vHmelt; see <tt>updateYieldStressFromHmelt()</tt>.

In this context, the current procedure uses a sliding velocity 
and a basal shear stress (sign choice so that it is stress applied to the ice)
to compute the till friction angle.  It is a three step process:
  - It is not assumed that the vectors \f$U\f$ and \f$\tau_{(b)}\f$
    at each grid point are pointing in the same direction 
    on input to this routine.  If the velocity
    is below 1 m/a in magnitude then \f$\tau_c\f$ is returned as 
    \f$|\tau_{(b)}|\f$ and then the till friction angle is computed
    by the above formula. If the velocity is larger in magnitude than
    1 m/a then the angle between \f$U\f$ and \f$\tau_{(b)}\f$ is checked; 
    if that angle exceeds 45 degrees then \f$\tau_c\f$ is returned as 
    \f$|\tau_{(b)}|\f$ and \f$\phi\f$ is computed as above.
  - Compute the yield stress through magnitudes:
    \f[ \tau_c = \frac{|\tau_{(b)}|\,|U_{\mathtt{thresh}}|^q}{|U|^q} \f]
  - Compute the till friction angle
    \f[ \phi = \frac{180}{\pi} \arctan\left(\frac{\tau_c - c_0}{N}\right) \f]

Values of both \f$\tau_c\f$ and \f$\phi\f$ are returned by this routine
but it is envisioned that the value of \f$\tau_c\f$ may evolve in a run
even if \f$\phi\f$ is treated as time-independent.
 */
PetscErrorCode IceModel::computeTFAFromBasalShearStressUsingPseudoPlastic(
                 const Vec myu, const Vec myv, const Vec mytaubx, const Vec mytauby, 
                 Vec &tauc_out, Vec &tfa_out) {
  PetscErrorCode ierr;
  SETERRQ(1,"NOT IMPLEMENTED");
  
  const PetscScalar sufficientSpeed = 1.0 / secpera,
                    sameDirAngle = pi/4.0;
  PetscScalar **tauc, **phi, **u, **v, **taubx, **tauby, **Hmelt, **H;
  ierr = DAVecGetArray(grid.da2, myu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, myv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, mytaubx, &taubx); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, mytauby, &tauby); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, tauc_out, &tauc); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, tfa_out, &phi); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar speed = sqrt(PetscSqr(u[i][j]) + PetscSqr(v[i][j])),
                        taumag = sqrt(PetscSqr(taubx[i][j]) + PetscSqr(tauby[i][j]));
      if (speed < sufficientSpeed) {
      } else {
/* FIXME:      
        const PetscScalar cosAngle = ;
        
        // this code is from updateYieldStressFromHmelt():
          const PetscScalar
                   overburdenP = ice->rho * grav * H[i][j],
                   pwP = plastic_till_pw_fraction * overburdenP,
                   lambda = Hmelt[i][j] / Hmelt_max,
                   N = overburdenP - lambda * pwP;  // effective pressure on till
*/
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, myu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, myv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, mytaubx, &taubx); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, mytauby, &tauby); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, tauc_out, &tauc); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, tfa_out, &phi); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  return 0;
}

