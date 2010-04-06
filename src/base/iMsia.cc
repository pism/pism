// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
There are two methods for computing the surface gradient.  

The default method is to directly differentiate 
the surface elevation \f$h\f$ by the Mahaffy method \ref Mahaffy .

The alternative method, using option <c>-eta</c> which sets 
<c>transformForSurfaceGradient = PETSC_TRUE</c>, is to transform the thickness 
to something more regular and differentiate that.  We get back to the gradient 
of the surface by applying the chain rule.  In particular, as shown 
in \ref CDDSV for the flat bed and \f$n=3\f$ case, if we define
	\f[\eta = H^{(2n+2)/n}\f]
then \f$\eta\f$ is more regular near the margin than \f$H\f$.  So we compute
the surface gradient by
   \f[\nabla h = \frac{n}{(2n+2)} \eta^{(-n-2)/(2n+2)} \nabla \eta + \nabla b,\f]
recalling that \f$h = H + b\f$.  This method is only applied when \f$\eta > 0\f$
at a given point; otherwise \f$\nabla h = \nabla b\f$.

In both cases we are computing the gradient by finite differences onto 
a staggered grid.  In the method with \f$\eta\f$ we apply centered differences
using (roughly) the same method for \f$\eta\f$ and \f$b\f$ that applies 
directly to the surface elevation \f$h\f$ in the alternative method.

The resulting surface gradient on the staggered grid is in four \c Vecs,
<c>vWork2d[k]</c> for \c k=0,1,2,3; recall there are two staggered grid points
per regular grid point.  These \c Vecs are used in velocitySIAStaggered(),
basalSlidingHeatingSIA(), and horizontalVelocitySIARegular().
 */
PetscErrorCode IceModel::surfaceGradientSIA() {
  PetscErrorCode  ierr;

  const PetscScalar   dx=grid.dx, dy=grid.dy;
  PetscScalar **h_x[2], **h_y[2];

  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);

  if (transformForSurfaceGradient == PETSC_TRUE) {
    PetscScalar **eta, **bed, **H;
    const PetscScalar n = ice->exponent(), // presumably 3.0
                      etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
                      invpow  = 1.0 / etapow,
                      dinvpow = (- n - 2.0) / (2.0 * n + 2.0);
    // compute eta = H^{8/3}, which is more regular, on reg grid
    ierr = vWork2d[4].get_array(eta); CHKERRQ(ierr);
    ierr = vH.get_array(H); CHKERRQ(ierr);
#ifdef LOCAL_GHOST_UPDATE
    PetscInt GHOSTS = 2;
#else
    PetscInt GHOSTS = 0;
#endif
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        eta[i][j] = pow(H[i][j], etapow);
      }
    }
    ierr = vWork2d[4].end_access(); CHKERRQ(ierr);
    ierr = vH.end_access(); CHKERRQ(ierr);
    // communicate eta: other processors will need ghosted for d/dx and d/dy
#ifndef LOCAL_GHOST_UPDATE
    ierr = vWork2d[4].beginGhostComm(); CHKERRQ(ierr);
    ierr = vWork2d[4].endGhostComm(); CHKERRQ(ierr);
#endif

    // now use Mahaffy on eta to get grad h on staggered;
    // note   grad h = (3/8) eta^{-5/8} grad eta + grad b  because  h = H + b
    ierr = vbed.get_array(bed); CHKERRQ(ierr);
    ierr = vWork2d[4].get_array(eta); CHKERRQ(ierr);
    for (PetscInt o=0; o<2; o++) {
#ifdef LOCAL_GHOST_UPDATE
      GHOSTS = 1;
#else
      GHOSTS = 0;
#endif
      for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
	for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
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
            h_x[o][i][j] += (bed[i+1][j] - bed[i][j]) / dx;
            h_y[o][i][j] += (+ bed[i+1][j+1] + bed[i][j+1]
                             - bed[i+1][j-1] - bed[i][j-1]) / (4.0*dy);
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
            h_y[o][i][j] += (bed[i][j+1] - bed[i][j]) / dy;
            h_x[o][i][j] += (+ bed[i+1][j+1] + bed[i+1][j]
                             - bed[i-1][j+1] - bed[i-1][j]) / (4.0*dx);
          }
        }
      }
    }
    ierr = vWork2d[4].end_access(); CHKERRQ(ierr);
    ierr = vbed.end_access(); CHKERRQ(ierr);
  } else {  // if transformForSurfaceGradient == FALSE; the Mahaffy scheme
    PetscScalar **h;
    ierr = vh.get_array(h); CHKERRQ(ierr);
    for (PetscInt o=0; o<2; o++) {
#ifdef LOCAL_GHOST_UPDATE
      PetscInt GHOSTS = 1;
#else
      PetscInt GHOSTS = 0;
#endif
      for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
	for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
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
    ierr = vh.end_access(); CHKERRQ(ierr);
  } // end if (transformForSurfaceGradient)
  
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);

#ifndef LOCAL_GHOST_UPDATE
  for (int j = 0; j < 4; ++j) {
    ierr = vWork2d[j].beginGhostComm(); CHKERRQ(ierr);
    ierr = vWork2d[j].endGhostComm(); CHKERRQ(ierr);
  }
#endif

  return 0;
}


//!  Compute the vertically-averaged horizontal velocity according to the non-sliding SIA.
/*!
See the comment for massContExplicitStep() before reading the rest of this comment.

Note that one may write 
  \f[ \mathbf{q} = \bar{\mathbf{U}} H = D \nabla h + \mathbf{U}_b \cdot H\f]
in shallow ice approximation (SIA) areas.  Here \f$h\f$ is the surface elevation of the ice
\f$\mathbf{U}_b\f$ is the basal sliding velocity, and \f$D\f$ is the diffusivity (which 
is computed in this method).

At the end of this routine the value of \f$D\f$ and of the <em>deformational part of</em> 
the vertically-averaged horizontal velocity, namely \f$D \nabla h\f$, is known at all staggered 
grid points.  It is stored in an IceModelVec2Stag called \c uvbar.

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

  PetscScalar *delta, *I, *J, *K, *Sigma;
  delta = new PetscScalar[grid.Mz];
  I = new PetscScalar[grid.Mz];
  J = new PetscScalar[grid.Mz];
  K = new PetscScalar[grid.Mz];
  Sigma = new PetscScalar[grid.Mz];

  const double enhancement_factor = config.get("enhancement_factor"),
               constant_grain_size = config.get("constant_grain_size");
  const bool   do_cold_ice = config.get_flag("do_cold_ice_methods");

  PetscScalar **h_x[2], **h_y[2], **H;

  PetscScalar *ageij, *ageoffset;

  if ((realAgeForGrainSize==PETSC_TRUE) && (!config.get_flag("do_age"))) {
    PetscPrintf(grid.com,
       "PISM ERROR in IceModel::velocitySIAStaggered(): do_age not set but\n"
       "age is needed for grain-size-based flow law ...  ENDING! ...\n\n");
    PetscEnd();
  }

  const bool usetau3 = (IceFlowLawUsesGrainSize(ice) && (realAgeForGrainSize==PETSC_TRUE) && (config.get_flag("do_age")));

  const PetscTruth usesGrainSize = IceFlowLawUsesGrainSize(ice);
  
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);
  ierr = uvbar.begin_access(); CHKERRQ(ierr);

  if (usetau3) {
    ierr = tau3.begin_access(); CHKERRQ(ierr);
  }
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = Istag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Istag3[1].begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].begin_access(); CHKERRQ(ierr);
  
  PetscScalar *Enthij, *Enthoffset;
  PolyThermalGPBLDIce *gpbldi = NULL;
  if (config.get_flag("do_cold_ice_methods") == false) {
    gpbldi = dynamic_cast<PolyThermalGPBLDIce*>(ice);
    if (!gpbldi) {
      PetscPrintf(grid.com,
        "do_cold_ice_methods == false in IceMethod::velocitySIAStaggered()\n"
        "   but not using PolyThermalGPBLDIce ... ending ....\n");
      PetscEnd();
    }
  }
  ierr = Enth3.begin_access(); CHKERRQ(ierr);

  // staggered grid computation of: I, J, Sigma
  for (PetscInt o=0; o<2; o++) {
#ifdef LOCAL_GHOST_UPDATE
    PetscInt GHOSTS = 1;
#else
    PetscInt GHOSTS = 0;
#endif
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        // staggered point: o=0 is i+1/2, o=1 is j+1/2,
        //   (i,j) and (i+oi,j+oj) are reg grid neighbors of staggered pt:
        const PetscInt     oi = 1-o, oj=o;  
        const PetscScalar  slope = (o==0) ? h_x[o][i][j] : h_y[o][i][j];
        const PetscScalar  thickness = 0.5 * (H[i][j] + H[i+oi][j+oj]);

        if (thickness > 0) { 
          if (usetau3) {
            ierr = tau3.getInternalColumn(i,j,&ageij); CHKERRQ(ierr);
            ierr = tau3.getInternalColumn(i+oi,j+oj,&ageoffset); CHKERRQ(ierr);
          }
	  
	  ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
	  ierr = Enth3.getInternalColumn(i+oi,j+oj,&Enthoffset); CHKERRQ(ierr);

          // does validity check for thickness:
          const PetscInt      ks = grid.kBelowHeight(thickness);  
          const PetscScalar   alpha =
                  sqrt(PetscSqr(h_x[o][i][j]) + PetscSqr(h_y[o][i][j]));

          I[0] = 0;   J[0] = 0;   K[0] = 0;
          for (PetscInt k=0; k<=ks; ++k) {
	    // pressure added by the ice (i.e. pressure difference between the
	    // current level and the top of the column)
            const PetscScalar   pressure = ice->rho * standard_gravity * (thickness - grid.zlevels[k]);
            PetscScalar flow, grainsize = constant_grain_size;
            if (usetau3 && usesGrainSize && realAgeForGrainSize) {
              grainsize = grainSizeVostok(0.5 * (ageij[k] + ageoffset[k]));
            }
            // If the flow law does not use grain size, it will just ignore it, no harm there
	    PetscScalar E = 0.5 * (Enthij[k] + Enthoffset[k]);
            if (do_cold_ice) {
	      PetscScalar T;
	      ierr = EC->getAbsTemp(E, pressure, T); CHKERRQ(ierr);
              flow = ice->flow(alpha * pressure, T, pressure, grainsize);
            } else {
              flow = gpbldi->flowFromEnth(alpha * pressure, E, pressure, grainsize);
            }

            delta[k] = 2.0 * pressure * enhancement_factor * flow;

            // for Sigma, ignore mask value and assume SHEET; will be overwritten
            // by correctSigma() in iMssa.cc
            Sigma[k] = delta[k] * PetscSqr(alpha) * pressure;

            if (k>0) { // trapezoid rule for I[k] and K[k]
              const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
              I[k] = I[k-1] + 0.5 * dz * (delta[k-1] + delta[k]);
              K[k] = K[k-1] + 0.5 * dz * (grid.zlevels[k-1] * delta[k-1]
                                          + grid.zlevels[k] * delta[k]);
              J[k] = grid.zlevels[k] * I[k] - K[k];
            }
          }
          for (PetscInt k=ks+1; k<grid.Mz; ++k) { // above the ice
            Sigma[k] = 0.0;
            I[k] = I[ks];
            J[k] = grid.zlevels[k] * I[ks];
          }  

          // diffusivity for deformational flow (vs basal diffusivity, incorporated in ub,vb)
          const PetscScalar  Dfoffset = J[ks] + (thickness - grid.zlevels[ks]) * I[ks];

          // vertically-averaged SIA-only velocity, sans sliding;
          //   note uvbar(i,j,0) is  u  at right staggered point (i+1/2,j)
          //   but uvbar(i,j,1) is  v  at up staggered point (i,j+1/2)
          uvbar(i,j,o) = - Dfoffset * slope / thickness;
         
          ierr = Istag3[o].setInternalColumn(i, j, I); CHKERRQ(ierr);
          ierr = Sigmastag3[o].setInternalColumn(i, j, Sigma); CHKERRQ(ierr);
        } else {  // zero thickness case
          uvbar(i,j,o) = 0;
          ierr = Istag3[o].setColumn(i,j,0.0); CHKERRQ(ierr);
          ierr = Sigmastag3[o].setColumn(i,j,0.0); CHKERRQ(ierr);
        } 
      } // o
    } // j
  } // i

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = uvbar.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);

  if (usetau3) {
    ierr = tau3.end_access(); CHKERRQ(ierr);
  }
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].end_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].end_access(); CHKERRQ(ierr);
  ierr = Istag3[0].end_access(); CHKERRQ(ierr);
  ierr = Istag3[1].end_access(); CHKERRQ(ierr);

  ierr = Enth3.end_access(); CHKERRQ(ierr);

  delete [] delta;   delete [] I;   delete [] J;   delete [] K;   delete [] Sigma;

#ifndef LOCAL_GHOST_UPDATE
  ierr = uvbar.beginGhostComm(); CHKERRQ(ierr);
  ierr = uvbar.endGhostComm(); CHKERRQ(ierr);

    ierr = Sigmastag3[0].beginGhostComm(); CHKERRQ(ierr);
    ierr = Sigmastag3[1].beginGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[0].beginGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[1].beginGhostComm(); CHKERRQ(ierr);
    ierr = Sigmastag3[0].endGhostComm(); CHKERRQ(ierr);
    ierr = Sigmastag3[1].endGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[0].endGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[1].endGhostComm(); CHKERRQ(ierr);
#endif

  
  return 0;
}


//! Compute the basal sliding and frictional heating if (where) SIA sliding rule is used.
/*!
This routine is only called, by velocity(), if \f$\mu\f$=\c mu_sliding is
non-zero.

THIS KIND OF SIA SLIDING LAW IS A BAD IDEA.  THAT'S WHY \f$\mu\f$ IS SET TO 
ZERO BY DEFAULT.  See Appendix B of \ref BBssasliding for the dangers in this
mechanism.

This routine calls the SIA-type sliding law, which may return zero in the
frozen base case; see basalVelocitySIA().  The basal sliding velocity is
computed for all SIA points.  This routine also computes the basal frictional
heating.  The basal velocity \c Vecs \c vub and \c vvb and the frictional
heating \c Vec are fully over-written.  Where the ice is floating, they all
have value zero.  

See correctBasalFrictionalHeating() for the SSA contribution.
 */
PetscErrorCode IceModel::basalSlidingHeatingSIA() {
  PetscErrorCode  ierr;
  PetscScalar **h_x[2], **h_y[2], **Rb, **H;

  PISMVector2 **bvel;

  double mu_sliding = config.get("mu_sliding");
  double minimum_temperature_for_sliding = config.get("minimum_temperature_for_sliding");

  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);

  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  ierr = vel_basal.get_array(bvel); CHKERRQ(ierr);
  ierr = vRb.get_array(Rb); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  
  for (PetscInt o=0; o<2; o++) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
        if (vMask.is_floating(i,j)) {
          bvel[i][j].u = 0.0;
          bvel[i][j].v = 0.0;
          Rb[i][j] = 0.0;
        } else {
          // basal velocity from SIA-type sliding law: not recommended!
          const PetscScalar
                  myx = -grid.Lx + grid.dx * i, 
                  myy = -grid.Ly + grid.dy * j,
                  myhx = 0.25 * (  h_x[0][i][j] + h_x[0][i-1][j]
                                 + h_x[1][i][j] + h_x[1][i][j-1]),
                  myhy = 0.25 * (  h_y[0][i][j] + h_y[0][i-1][j]
                                 + h_y[1][i][j] + h_y[1][i][j-1]),
	    alpha = sqrt(PetscSqr(myhx) + PetscSqr(myhy));
	  PetscScalar T, basalC;

	  ierr = EC->getAbsTemp(Enth3.getValZ(i,j,0.0),
				EC->getPressureFromDepth(H[i][j]), T); CHKERRQ(ierr);

	  basalC = basalVelocitySIA(myx, myy, H[i][j], T, 
				    alpha, mu_sliding,
				    minimum_temperature_for_sliding);
          bvel[i][j].u = - basalC * myhx;
          bvel[i][j].v = - basalC * myhy;
          // basal frictional heating; note P * dh/dx is x comp. of basal shear stress
          // in ice streams this result will be *overwritten* by
          //   correctBasalFrictionalHeating() if useSSAVelocities==TRUE
          const PetscScalar P = ice->rho * standard_gravity * H[i][j];
          Rb[i][j] = - (P * myhx) * bvel[i][j].u - (P * myhy) * bvel[i][j].v;
        }
      }
    }
  }
  
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);

  ierr = vel_basal.end_access(); CHKERRQ(ierr);
  ierr = vRb.end_access(); CHKERRQ(ierr);
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Average staggered-grid vertically-averaged horizontal velocity onto regular grid.
/*! 
At the end of velocitySIAStaggered() the vertically-averaged horizontal velocity 
components from deformation (stored in uvbar) are known on the regular grid.
At the end of basalSIA() the basal sliding from an SIA-type sliding rule is in
vel_basal.  This procedure averages the former onto the regular grid and adds
the sliding velocity.

That is, this procedure computes the SIA "first guess" at the
vertically-averaged horizontal velocity.  Therefore the values in \c Vec\ s
\c vel_bar are merely tentative.  The values in \c uvbar are, however,
PISM's estimate of \e deformation by shear in vertical planes.

Note that communication of ghosted values must occur between calling
velocitySIAStaggered() and this procedure for the averaging to work.  Only
two-dimensional regular grid velocities are updated here.  The full
three-dimensional velocity field is not updated here but instead in
horizontalVelocitySIARegular() and in vertVelocityFromIncompressibility().
 */
PetscErrorCode IceModel::velocities2DSIAToRegular() {  
  PetscErrorCode ierr;

  double mu_sliding = config.get("mu_sliding");

  ierr = vel_bar.begin_access(); CHKERRQ(ierr);
  ierr = uvbar.begin_access(); CHKERRQ(ierr);
  if (mu_sliding == 0.0) {
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        // compute ubar,vbar on regular grid by averaging deformational onto
        //   staggered grid
        vel_bar(i,j).u = 0.5*(uvbar(i-1,j,0) + uvbar(i,j,0));
        vel_bar(i,j).v = 0.5*(uvbar(i,j-1,1) + uvbar(i,j,1));
      }
    }
  } else {
    // this case is not recommended!
    PISMVector2 **bvel;
    ierr = vel_basal.get_array(bvel); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        // as above, adding basal on regular grid
        vel_bar(i,j).u = 0.5*(uvbar(i-1,j,0) + uvbar(i,j,0)) + bvel[i][j].u;
        vel_bar(i,j).v = 0.5*(uvbar(i,j-1,1) + uvbar(i,j,1)) + bvel[i][j].v;
      }
    }
    ierr = vel_basal.end_access(); CHKERRQ(ierr);
  }
  ierr = vel_bar.end_access(); CHKERRQ(ierr);
  ierr =   uvbar.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Put the volume strain heating (dissipation heating) onto the regular grid.
/*!
At the end of velocitySIAStaggered() the volume strain-heating \f$\Sigma\f$ is
available on the staggered grid.  This procedure averages it onto the regular
grid.  \f$\Sigma\f$ is used in the temperature equation.

Communication of ghosted values of \c Vec \c vSigma must occur between 
velocitySIAStaggered() and this procedure for the averaging to work.
 */
PetscErrorCode IceModel::SigmaSIAToRegular() {
  PetscErrorCode  ierr;
  PetscScalar **H;

  PetscScalar *Sigmareg, *SigmaEAST, *SigmaWEST, *SigmaNORTH, *SigmaSOUTH;

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = Sigma3.begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
        // horizontally average Sigma onto regular grid
        const PetscInt ks = grid.kBelowHeight(H[i][j]);
        ierr = Sigma3.getInternalColumn(i,j,&Sigmareg); CHKERRQ(ierr);
        ierr = Sigmastag3[0].getInternalColumn(i,j,&SigmaEAST); CHKERRQ(ierr);
        ierr = Sigmastag3[0].getInternalColumn(i-1,j,&SigmaWEST); CHKERRQ(ierr);
        ierr = Sigmastag3[1].getInternalColumn(i,j,&SigmaNORTH); CHKERRQ(ierr);
        ierr = Sigmastag3[1].getInternalColumn(i,j-1,&SigmaSOUTH); CHKERRQ(ierr);
        for (PetscInt k = 0; k <= ks; ++k) {
          Sigmareg[k] = 0.25 * (SigmaEAST[k] + SigmaWEST[k] + SigmaNORTH[k] + SigmaSOUTH[k]);
        }
        for (PetscInt k = ks+1; k < grid.Mz; ++k) {
          Sigmareg[k] = 0.0;
        }
        // no need to call Sigma3.setInternalColumn(); already set!
      } else { // zero thickness case
        ierr = Sigma3.setColumn(i,j,0.0); CHKERRQ(ierr);
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = Sigma3.end_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].end_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].end_access(); CHKERRQ(ierr);
  
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
  PetscScalar **h_x[2], **h_y[2];
  PetscScalar *u, *v, *IEAST, *IWEST, *INORTH, *ISOUTH;

  u = new PetscScalar[grid.Mz];
  v = new PetscScalar[grid.Mz];  

  double mu_sliding = config.get("mu_sliding");

  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = Istag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Istag3[1].begin_access(); CHKERRQ(ierr);
  ierr = vel_basal.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = Istag3[0].getInternalColumn(i,j,&IEAST); CHKERRQ(ierr);
      ierr = Istag3[0].getInternalColumn(i-1,j,&IWEST); CHKERRQ(ierr);
      ierr = Istag3[1].getInternalColumn(i,j,&INORTH); CHKERRQ(ierr);
      ierr = Istag3[1].getInternalColumn(i,j-1,&ISOUTH); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        u[k] = - 0.25 * ( IEAST[k] * h_x[0][i][j] + IWEST[k] * h_x[0][i-1][j] +
                          INORTH[k] * h_x[1][i][j] + ISOUTH[k] * h_x[1][i][j-1] );
        v[k] = - 0.25 * ( IEAST[k] * h_y[0][i][j] + IWEST[k] * h_y[0][i-1][j] +
                          INORTH[k] * h_y[1][i][j] + ISOUTH[k] * h_y[1][i][j-1] );
      }

      if (mu_sliding > 0.0) {	// unusual case
        for (PetscInt k=0; k<grid.Mz; ++k) {
          u[k] += vel_basal(i,j).u;
	  v[k] += vel_basal(i,j).v;
        }
      }

      ierr = u3.setInternalColumn(i, j, u); CHKERRQ(ierr);
      ierr = v3.setInternalColumn(i, j, v); CHKERRQ(ierr);
    }
  }

  delete[] u;
  delete[] v;

  ierr = vel_basal.end_access(); CHKERRQ(ierr);
  
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  ierr = Istag3[0].end_access(); CHKERRQ(ierr);
  ierr = Istag3[1].end_access(); CHKERRQ(ierr);

  return 0;
}


//! Compute the coefficient of surface gradient, for basal sliding velocity as a function of driving stress in SIA regions.
/*!
THIS KIND OF SIA SLIDING LAW IS A BAD IDEA IN A THERMOMECHANICALLY-COUPLED
MODEL.  THAT'S WHY \f$\mu\f$ IS SET TO ZERO BY DEFAULT.                

In SIA regions (= MASK_SHEET) a basal sliding law of the form
  \f[ \mathbf{U}_b = (u_b,v_b) = - C \nabla h \f] 
is allowed.  Here \f$\mathbf{U}_b\f$ is the horizontal velocity of the base of
the ice (the "sliding velocity") and \f$h\f$ is the elevation of the ice
surface.  This procedure returns the \em positive \em coefficient \f$C\f$ in
this relationship.  This coefficient can depend of the thickness, the basal
temperature, and the horizontal location.

The default version for IceModel here is location-independent 
pressure-melting-temperature-activated linear sliding.  See Appendix B of
\ref BBssasliding for the dangers in this mechanism.

Parameter \f$\mu\f$ can be set by option \c -mu_sliding.

The returned coefficient is used in basalSlidingHeatingSIA().

This procedure is virtual and can be replaced by any derived class.
 */
PetscScalar IceModel::basalVelocitySIA(
               PetscScalar /*x*/, PetscScalar /*y*/, PetscScalar H, PetscScalar T,
               PetscScalar /*alpha*/, PetscScalar mu, PetscScalar min_T) const {
  if (T + ice->beta_CC_grad * H > min_T) {
    const PetscScalar p_over = ice->rho * standard_gravity * H;
    return mu * p_over;
  } else {
    return 0;
  }
}

