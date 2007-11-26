// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

#include <cmath>
#include <petscda.h>
#include "iceModel.hh"


//! Compute the surface gradient in advance of the SIA velocity computation.
/*! 
There are two methods for computing the surface gradient.  The default is to transform the 
thickness to something more regular and differentiate that.  In particular, as shown 
in (Calvo et al 2002) for the flat bed and \f$n=3\f$ case, if we define
	\f[\eta = H^{(2n+2)/n}\f]
then \f$\eta\f$ is more regular near the margin than \f$H\f$.  So the default method for computing
the surface gradient is to compute
   \f[\nabla h = \frac{n}{(2n+2)} \eta^{(-n-2)/(2n+2)} \nabla \eta + \nabla b,\f]
recalling that \f$h = H + b\f$.  This method is only applied when \f$\eta > 0\f$ at a given point;
otherwise \f$\nabla h = \nabla b\f$.

We are computing this gradient by finite differences onto a staggered grid.  We do so 
by centered differences using (roughly) the same method for \f$\eta\f$ and 
\f$b\f$ that (Mahaffy 1976) applies directly to the surface elevation \f$h\f$.

The optional method is to directly differentiate the surface elevation \f$h\f$ by the
(Mahaffy 1976) method.
 */
PetscErrorCode IceModel::surfaceGradientSIA() {
  PetscErrorCode  ierr;

  const PetscScalar   dx=grid.p->dx, dy=grid.p->dy;
  PetscScalar **h_x[2], **h_y[2];

  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);

  if (transformForSurfaceGradient == PETSC_TRUE) {
    PetscScalar **eta, **b, **H;
    const PetscScalar n = isothermalFlux_n_exponent; // presumably 3.0
    const PetscScalar etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
                      invpow  = 1.0 / etapow,
                      dinvpow = (- n - 2.0) / (2.0 * n + 2.0);
    // compute eta = H^{8/3}, which is more regular, on reg grid
    ierr = DAVecGetArray(grid.da2, vWork2d[4], &eta); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        eta[i][j] = pow(H[i][j], etapow);
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vWork2d[4], &eta); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    // communicate eta: other processors will need ghosted for d/dx and d/dy
    ierr = DALocalToLocalBegin(grid.da2, vWork2d[4], INSERT_VALUES, vWork2d[4]); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vWork2d[4], INSERT_VALUES, vWork2d[4]); CHKERRQ(ierr);
    // now use Mahaffy on eta to get grad h on staggered;
    // note   grad h = (3/8) eta^{-5/8} grad eta + grad b  because  h = H + b
    ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[4], &eta); CHKERRQ(ierr);
    for (PetscInt o=0; o<2; o++) {
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
        for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
          if (o==0) {     // If I-offset
            const PetscScalar mean_eta = 0.5 * (eta[i+1][j] + eta[i][j]);
            if (mean_eta > 0.0) {
              const PetscScalar factor = invpow * pow(mean_eta, dinvpow);
              h_x[o][i][j] = factor * (eta[i+1][j] - eta[i][j]) / dx;
              h_y[o][i][j] = factor * (+ eta[i+1][j+1] + eta[i][j+1]
                                     - eta[i+1][j-1] - eta[i][j-1]) / (4.0*dy);
            } else {
              h_x[o][i][j] = 0.0;
              h_y[o][i][j] = 0.0;
            }
            // now add bed slope to get actual h_x,h_y
            h_x[o][i][j] += (b[i+1][j] - b[i][j]) / dx;
            h_y[o][i][j] += (+ b[i+1][j+1] + b[i][j+1]
                             - b[i+1][j-1] - b[i][j-1]) / (4.0*dy);
          } else {        // J-offset
            const PetscScalar mean_eta = 0.5 * (eta[i][j+1] + eta[i][j]);
            if (mean_eta > 0.0) {
              const PetscScalar factor = invpow * pow(mean_eta, dinvpow);
              h_y[o][i][j] = factor * (eta[i][j+1] - eta[i][j]) / dy;
              h_x[o][i][j] = factor * (+ eta[i+1][j+1] + eta[i+1][j]
                                     - eta[i-1][j+1] - eta[i-1][j]) / (4.0*dx);
            } else {
              h_y[o][i][j] = 0.0;
              h_x[o][i][j] = 0.0;
            }
            // now add bed slope to get actual h_x,h_y
            h_y[o][i][j] += (b[i][j+1] - b[i][j]) / dy;
            h_x[o][i][j] += (+ b[i+1][j+1] + b[i+1][j]
                             - b[i-1][j+1] - b[i-1][j]) / (4.0*dx);
          }
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vWork2d[4], &eta); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  } else {  // if !transformForSurfaceGradient; the old way
    PetscScalar **h;
    ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
    for (PetscInt o=0; o<2; o++) {
      for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
        for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
          if (o==0) {     // If I-offset
            h_x[o][i][j] = (h[i+1][j] - h[i][j]) / dx;
            h_y[o][i][j] = (+ h[i+1][j+1] + h[i][j+1]
                            - h[i+1][j-1] - h[i][j-1]) / (4.0*dy);
          } else {        // J-offset
            h_y[o][i][j] = (h[i][j+1] - h[i][j]) / dy;
            h_x[o][i][j] = (+ h[i+1][j+1] + h[i+1][j]
                            - h[i-1][j+1] - h[i-1][j]) / (4.0*dx);
          }
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  } // end if (transformForSurfaceGradient)
  
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  return 0;
}


//!  Compute the vertically-averaged horizontal velocity according to the non-sliding SIA.
/*!
See the comment for massBalExplicitStep() before reading the rest of this comment.

Note that one may write 
  \f[ \mathbf{q} = \bar{\mathbf{U}} H = D \nabla h + \mathbf{U}_b \cdot H\f]
in shallow ice approximation (SIA) areas.  Here \f$h\f$ is the surface elevation of the ice
\f$\mathbf{U}_b\f$ is the basal sliding velocity, and \f$D\f$ is the diffusivity (which 
is computed in this method).

At the end of this routine the value of \f$D\f$ and of the <em>deformational part of</em> 
the vertically-averaged horizontal velocity, namely \f$D \nabla h\f$, is known at all staggered 
grid points.  It is stored in the \c Vec pair called \c vuvbar.

The scheme used for this is
the one first proposed in the context of ice sheets by Mahaffy (1976).  That is, the method 
is "type I" in the classification described in (Hindmarsh and Payne 1996).  Note that the 
surface slope \f$\nabla h\f$ is needed on the staggered grid although the surface 
elevation \f$h\f$ itself is known on the regular grid.  

This routine also computes the (ice volume but not basal) strain-heating.  In particular,
the staggered grid value of \f$\Sigma\f$ is computed using the formula appropriate to the SIA
case and is put in a workspace \c Vec.  See correctSigma().
 */
PetscErrorCode IceModel::velocitySIAStaggered() {
  PetscErrorCode  ierr;

  const PetscScalar   dz=grid.p->dz;
  PetscScalar *delta, *K, *J[2];
  delta = new PetscScalar[grid.p->Mz];
  K = new PetscScalar[grid.p->Mz];
  J[0] = new PetscScalar[grid.p->Mz];
  J[1] = new PetscScalar[grid.p->Mz];

  PetscScalar **h_x[2], **h_y[2], **H, **uvbar[2];
  PetscScalar ***I[2], ***Sigma[2], ***T, ***gs;

  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[2], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[3], &Sigma[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
  
  // staggered grid computation of: I, J, Sigma
  for (PetscInt o=0; o<2; o++) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        // staggered point: o==0 is right, o==1 is up
        const PetscInt      oi = 1-o, oj=o;

        PetscScalar  slope;
        if (o==0) {     // If I-offset
          slope = h_x[o][i][j];
        } else {        // J-offset
          slope = h_y[o][i][j];
        }
        const PetscScalar   thickness = 0.5 * (H[i][j] + H[i+oi][j+oj]);
 
        if (thickness > 0) { 
          const PetscInt      ks = static_cast<PetscInt>(floor(thickness/dz));
          const PetscScalar   alpha =
                  sqrt(PetscSqr(h_x[o][i][j]) + PetscSqr(h_y[o][i][j]));

          if (ks>grid.p->Mz) {
            ierr = PetscPrintf(grid.com,
                 "[[error LOCATION: i, j, ks, H = %5d %5d %5d %10.2f]]\n",i, j, ks, H[i][j]); 
            SETERRQ(1, "thickness overflow in SIA velocity: ks>Mz");
          }
          I[o][i][j][0] = 0;   J[o][0] = 0;   K[0] = 0;
          for (PetscInt k=0; k<=ks; ++k) {
            const PetscScalar   s = k * dz;
            const PetscScalar   pressure = ice.rho * grav * (thickness - s);
            delta[k] = (2 * pressure * enhancementFactor
                        * ice.flow(alpha * pressure,
                                   0.5 * (T[i][j][k] + T[i+oi][j+oj][k]), pressure,
                                   0.5 * (gs[i][j][k] + gs[i+oi][j+oj][k])));
            // for Sigma, ignor mask value and assume SHEET; will be overwritten
            // by correctSigma() in iMssa.cc
            Sigma[o][i][j][k] = delta[k] * PetscSqr(alpha) * pressure
                                  / (ice.rho * ice.c_p);
            if (k>0) { // trapezoid rule for I[][][][k] and K[k]
              I[o][i][j][k] = I[o][i][j][k-1] + 0.5 * dz * (delta[k-1] + delta[k]);
              K[k] = K[k-1] + 0.5 * dz * ((s-dz)*delta[k-1] + s*delta[k]);
              J[o][k] = s * I[o][i][j][k] - K[k];
            }
          }
          for (PetscInt k=ks+1; k<grid.p->Mz; ++k) { // above the ice
            Sigma[o][i][j][k] = 0.0;
            I[o][i][j][k] = I[o][i][j][ks];
            J[o][k] = k * dz * I[o][i][j][ks]; // = J[o][ks];
          }  

          // diffusivity for deformational flow (vs basal diffusivity, incorporated in ub,vb)
          const PetscScalar
            Dfoffset = J[o][ks] + (thickness - ks * dz) * I[o][i][j][ks];
          // vertically-averaged velocity; note uvbar[0][i][j] is  u  at right staggered
          // point (i+1/2,j) but uvbar[1][i][j] is  v  at up staggered point (i,j+1/2)
          // here we use stale (old) diffusivity if faststep
          uvbar[o][i][j] = - Dfoffset * slope / thickness;
         
        } else {  // zero thickness case
          uvbar[o][i][j] = 0;
          for (PetscInt k=0; k < grid.p->Mz; k++) {
            I[o][i][j][k] = 0;
            Sigma[o][i][j][k] = 0;
          }
        } 
      } // o
    } // j
  } // i

  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vgs, &gs); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[2], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[3], &Sigma[1]); CHKERRQ(ierr);

  delete [] delta;   delete [] K;   delete [] J[0];   delete [] J[1];
  return 0;
}


//! Compute the basal sliding and frictional heating if (where) SIA sliding rule is used.
/*! This routine calls the SIA-type sliding law, which may return zero in the frozen base
case.  I.e. basalVelocity().  The basal sliding velocity is computed for all SIA 
points.  This routine also computes the basal frictional heating.  

The basal velocity \c Vecs \c vub and \c vvb and the frictional heating \c Vec are all
fully over-written by this routine.  Where the ice is floating, they all have value zero.

See correctBasalFrictionalHeating().
 */
PetscErrorCode IceModel::basalSIA() {
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **ub, **vb, **Rb, **mask, **H;
  PetscScalar ***T;

  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vRb, &Rb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  
  for (PetscInt o=0; o<2; o++) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (modMask(mask[i][j]) == MASK_FLOATING) {
          ub[i][j] = 0.0;
          vb[i][j] = 0.0;
          Rb[i][j] = 0.0;
        } else { 
          // basal velocity
          const PetscScalar
                  myx = -grid.p->Lx + grid.p->dx * i, 
                  myy = -grid.p->Ly + grid.p->dy * j,
                  myhx = 0.25 * (  h_x[0][i][j] + h_x[0][i-1][j]
                                 + h_x[1][i][j] + h_x[1][i][j-1]),
                  myhy = 0.25 * (  h_y[0][i][j] + h_y[0][i-1][j]
                                 + h_y[1][i][j] + h_y[1][i][j-1]),
                  alpha = sqrt(PetscSqr(myhx) + PetscSqr(myhy)),
                  basalC = basalVelocity(myx, myy, H[i][j], T[i][j][0], alpha, muSliding);
          ub[i][j] = - basalC * myhx;
          vb[i][j] = - basalC * myhy;
        
          // basal frictional heating; note P * dh/dx is x comp. of basal shear stress
          // in ice streams this result will be *overwritten* by
          //   correctBasalFrictionalHeating() if useSSAVelocities==TRUE
          const PetscScalar P = ice.rho * grav * H[i][j];
          Rb[i][j] = - (P * myhx) * ub[i][j] - (P * myhy) * vb[i][j];
        }
      }
    }
  }
  
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vRb, &Rb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  return 0;
}


//! Average staggered-grid vertically-averaged horizontal velocity onto regular grid.
/*! 
At the end of velocitySIAStaggered() the vertically-averaged horizontal velocity 
vuvbar[0],vuvbar[1] from deformation is known on the regular grid.  At the end of basalSIA()
the basal sliding from an SIA-type sliding rule is in vub, vvb.  This procedure 
averages the former onto the regular grid and adds the sliding velocity.

That is, this procedure computes the SIA ``first guess'' at the vertically-averaged horizontal
velocity.  Therefore the values in \c Vec\ s \c vubar, \c vvbar are merely tentative.  The 
values in \c vuvbar are authoritative; these are PISM's estimate of \e deformation by shear
in vertical planes.

Note that communication of ghosted values must occur between velocitySIAStaggered() and this 
procedure for the averaging to work.  Only two-dimensional regular grid velocities 
are updated here.  The full three-dimensional velocity field is not updated here
but instead in horizontalVelocitySIARegular() and in vertVelocityFromIncompressibility().
 */
PetscErrorCode IceModel::velocities2DSIAToRegular() {  
  PetscErrorCode ierr;
  PetscScalar **ubar, **vbar, **uvbar[2], **ub, **vb;

  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // compute ubar,vbar on regular grid by averaging deformational on staggered grid
      // and adding basal on regular grid
      ubar[i][j] = 0.5*(uvbar[0][i-1][j] + uvbar[0][i][j]) + ub[i][j];
      vbar[i][j] = 0.5*(uvbar[1][i][j-1] + uvbar[1][i][j]) + vb[i][j];
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  return 0;
}


//! Put the (ice volume, not basal) strain-heating onto the regular grid.
/*!
At the end of velocitySIAStaggered() the volume strain-heating \f$\Sigma\f$ is available
on the staggered grid.  This procedure averages it onto the regular grid.

Note that communication of ghosted values of \c Vec \c vSigma must occur between 
velocitySIAStaggered() and this procedure for the averaging to work.
 */
PetscErrorCode IceModel::SigmaSIAToRegular() {
  // average Sigma onto regular grid for use in the temperature equation
  PetscErrorCode  ierr;
  PetscScalar ***Sigma[2], **H, ***Sigmareg;

  ierr = DAVecGetArray(grid.da3, vWork3d[2], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[3], &Sigma[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigma, &Sigmareg); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // horizontally average Sigma onto regular grid
      const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
      for (PetscInt k=0; k<ks; ++k) {
        Sigmareg[i][j][k] = 0.25 * (Sigma[0][i][j][k] + Sigma[0][i-1][j][k] +
                                    Sigma[1][i][j][k] + Sigma[1][i][j-1][k]);
      }
      for (PetscInt k=ks+1; k<grid.p->Mz; ++k) {
        Sigmareg[i][j][k] = 0.0;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da3, vWork3d[2], &Sigma[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[3], &Sigma[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigmareg); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  return 0;
}


//! Update regular grid horizontal velocities u,v at depth for SIA regions.
/*! 
The procedure velocitySIAStaggered() computes several scalar
quantities at depth (the details of which are too complicated to explain).  
These quantities correspond to three-dimensional arrays.  This procedure takes those 
quantities and computes the three-dimensional arrays for the horizontal components \f$u\f$ and 
\f$v\f$ of the velocity field.

The vertical component \f$w\f$ of the velocity field 
is computed later by vertVelocityFromIncompressibility().
 */
PetscErrorCode IceModel::horizontalVelocitySIARegular() {
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **ub, **vb;
  PetscScalar ***I[2], ***u, ***v;

  ierr = DAVecGetArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
    
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      for (PetscInt k=0; k<grid.p->Mz; ++k) {
        u[i][j][k] =  ub[i][j] - 0.25 *
          (I[0][i][j][k]*h_x[0][i][j] + I[0][i-1][j][k]*h_x[0][i-1][j] +
           I[1][i][j][k]*h_x[1][i][j] + I[1][i][j-1][k]*h_x[1][i][j-1]);
        v[i][j][k] =  vb[i][j] - 0.25 *
          (I[0][i][j][k]*h_y[0][i][j] + I[0][i-1][j][k]*h_y[0][i-1][j] +
           I[1][i][j][k]*h_y[1][i][j] + I[1][i][j-1][k]*h_y[1][i][j-1]);
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &h_x[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &h_x[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &h_y[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[3], &h_y[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[0], &I[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vWork3d[1], &I[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  return 0;
}

