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
There are three methods for computing the surface gradient.  Which method is
controlled by configuration parameter \c surface_gradient_method which can have
values \c haseloff, \c mahaffy, or \c eta.

The most traditional method is to directly differentiate 
the surface elevation \f$h\f$ by the Mahaffy method \ref Mahaffy.  The \c haseloff
method, suggested by Marianne Haseloff, modifies the Mahaffy method only where
ice-free adjacent bedrock points are above the ice surface, and in those cases
the returned gradient component is zero.

The alternative method, when \c surface_gradient_method = \c eta, transforms the thickness 
to something more regular and differentiates that.  We get back to the gradient 
of the surface by applying the chain rule.  In particular, as shown 
in \ref CDDSV for the flat bed and \f$n=3\f$ case, if we define
	\f[\eta = H^{(2n+2)/n}\f]
then \f$\eta\f$ is more regular near the margin than \f$H\f$.  So we compute
the surface gradient by
   \f[\nabla h = \frac{n}{(2n+2)} \eta^{(-n-2)/(2n+2)} \nabla \eta + \nabla b,\f]
recalling that \f$h = H + b\f$.  This method is only applied when \f$\eta > 0\f$
at a given point; otherwise \f$\nabla h = \nabla b\f$.

In all cases we are computing the gradient by finite differences onto 
a staggered grid.  In the method with \f$\eta\f$ we apply centered differences
using (roughly) the same method for \f$\eta\f$ and \f$b\f$ that applies 
directly to the surface elevation \f$h\f$ in the \c mahaffy and \c haseloff methods.

The resulting surface gradient on the staggered grid is put in four \c Vecs,
<c>vWork2d[k]</c> for \c k=0,1,2,3; recall there are two staggered grid points
per regular grid point and two scalar components to the vector gradient.
The surface gradient values stored in these \c Vecs are used in velocitySIAStaggered(),
basalSlidingHeatingSIA(), and horizontalVelocitySIARegular().
 */
PetscErrorCode IceModel::surfaceGradientSIA() {
  PetscErrorCode  ierr;

  const string method = config.get_string("surface_gradient_method");
  // here is where the three versions are implemented, so check that it is one
  //   of the three
  if ((method!="eta") && (method!="mahaffy") && (method!="haseloff")) {
    verbPrintf(1, grid.com,
      "PISM ERROR: value of surface_gradient_method, option -gradient, not valid ... ending\n");
    PetscEnd();
  }

  if (method == "eta") {

    ierr = surface_gradient_eta(); CHKERRQ(ierr); 

  } else if (method == "haseloff") {

    ierr = surface_gradient_haseloff(); CHKERRQ(ierr);

  } else if (method == "mahaffy") {

    ierr = surface_gradient_mahaffy(); CHKERRQ(ierr);

  } else {
    SETERRQ(1, "can't happen");
  }

  return 0;
}
/*!
 * Mary Anne Mahaffy method; see [\ref Mahaffy].
 */
PetscErrorCode IceModel::surface_gradient_mahaffy() {
  PetscErrorCode ierr;

  const PetscScalar dx = grid.dx, dy = grid.dy;  // convenience

  PetscScalar **h_x[2], **h_y[2], **h;
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);
  ierr = vh.get_array(h); CHKERRQ(ierr);

  for (PetscInt o=0; o<2; o++) {
    PetscInt GHOSTS = 1;
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

  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);

  return 0;
}

/*!
 * Marianne Haseloff method: deals correctly with adjacent ice-free points with
 * bed elevations which are above the surface of the ice
 */
PetscErrorCode IceModel::surface_gradient_haseloff() {
  PetscErrorCode ierr;

  const PetscScalar Hicefree = 0.0;  // standard for ice-free, in Haseloff
  const PetscScalar dx = grid.dx, dy = grid.dy;  // convenience

  PetscScalar **h_x[2], **h_y[2];
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);

  PetscScalar **h, **b, **H;
  ierr = vbed.get_array(b); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = vh.get_array(h); CHKERRQ(ierr);
  for (PetscInt o=0; o<2; o++) {

    PetscInt GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        if (o==0) {     // If I-offset
          const bool icefreeP  = (H[i][j]     <= Hicefree),
            icefreeE  = (H[i+1][j]   <= Hicefree),
            icefreeN  = (H[i][j+1]   <= Hicefree),
            icefreeS  = (H[i][j-1]   <= Hicefree),
            icefreeNE = (H[i+1][j+1] <= Hicefree),
            icefreeSE = (H[i+1][j-1] <= Hicefree);

          PetscScalar hhE = h[i+1][j];  // east pseudo-surface elevation
          if (icefreeE  && (b[i+1][j]   > h[i][j]    ))  hhE  = h[i][j];
          if (icefreeP  && (b[i][j]     > h[i+1][j]  ))  hhE  = h[i][j];
          h_x[o][i][j] = (hhE - h[i][j]) / dx;

          PetscScalar hhN  = h[i][j+1];  // north pseudo-surface elevation
          if (icefreeN  && (b[i][j+1]   > h[i][j]    ))  hhN  = h[i][j];
          if (icefreeP  && (b[i][j]     > h[i][j+1]  ))  hhN  = h[i][j];
          PetscScalar hhS  = h[i][j-1];  // south pseudo-surface elevation
          if (icefreeS  && (b[i][j-1]   > h[i][j]    ))  hhS  = h[i][j];
          if (icefreeP  && (b[i][j]     > h[i][j-1]  ))  hhS  = h[i][j];
          PetscScalar hhNE = h[i+1][j+1];// northeast pseudo-surface elevation
          if (icefreeNE && (b[i+1][j+1] > h[i+1][j]  ))  hhNE = h[i+1][j];
          if (icefreeE  && (b[i+1][j]   > h[i+1][j+1]))  hhNE = h[i+1][j];
          PetscScalar hhSE = h[i+1][j-1];// southeast pseudo-surface elevation
          if (icefreeSE && (b[i+1][j-1] > h[i+1][j]  ))  hhSE = h[i+1][j];
          if (icefreeE  && (b[i+1][j]   > h[i+1][j-1]))  hhSE = h[i+1][j];
          h_y[o][i][j] = (hhNE + hhN - hhSE - hhS) / (4.0 * dy);
        } else {        // J-offset
          const bool icefreeP  = (H[i][j]     <= Hicefree),
            icefreeN  = (H[i][j+1]   <= Hicefree),
            icefreeE  = (H[i+1][j]   <= Hicefree),
            icefreeW  = (H[i-1][j]   <= Hicefree),
            icefreeNE = (H[i+1][j+1] <= Hicefree),
            icefreeNW = (H[i-1][j+1] <= Hicefree);

          PetscScalar hhN  = h[i][j+1];  // north pseudo-surface elevation
          if (icefreeN  && (b[i][j+1]   > h[i][j]    ))  hhN  = h[i][j];
          if (icefreeP  && (b[i][j]     > h[i][j+1]  ))  hhN  = h[i][j];
          h_y[o][i][j] = (hhN - h[i][j]) / dy;

          PetscScalar hhE  = h[i+1][j];  // east pseudo-surface elevation
          if (icefreeE  && (b[i+1][j]   > h[i][j]    ))  hhE  = h[i][j];
          if (icefreeP  && (b[i][j]     > h[i+1][j]  ))  hhE  = h[i][j];
          PetscScalar hhW  = h[i-1][j];  // west pseudo-surface elevation
          if (icefreeW  && (b[i-1][j]   > h[i][j]    ))  hhW  = h[i][j];
          if (icefreeP  && (b[i][j]     > h[i-1][j]  ))  hhW  = h[i][j];
          PetscScalar hhNE = h[i+1][j+1];// northeast pseudo-surface elevation
          if (icefreeNE && (b[i+1][j+1] > h[i][j+1]  ))  hhNE = h[i][j+1];
          if (icefreeN  && (b[i][j+1]   > h[i+1][j+1]))  hhNE = h[i][j+1];
          PetscScalar hhNW = h[i-1][j+1];// northwest pseudo-surface elevation
          if (icefreeNW && (b[i-1][j+1] > h[i][j+1]  ))  hhNW = h[i][j+1];
          if (icefreeN  && (b[i][j+1]   > h[i-1][j+1]))  hhNW = h[i][j+1];
          h_x[o][i][j] = (hhNE + hhE - hhNW - hhW) / (4.0 * dx);
        }

      } // j
    }   // i
  }     // o
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vh.end_access(); CHKERRQ(ierr);

  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the surface gradient on the staggered grid using the eta
//! transformation.
PetscErrorCode IceModel::surface_gradient_eta() {
  PetscErrorCode ierr;

  PetscScalar **h_x[2], **h_y[2];
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);

  PetscScalar **eta;
  const PetscScalar n = ice->exponent(), // presumably 3.0
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0);
  const PetscScalar dx = grid.dx, dy = grid.dy;  // convenience

  // compute eta = H^{8/3}, which is more regular, on reg grid
  ierr = vWork2d[4].get_array(eta); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);

  PetscInt GHOSTS = 2;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      eta[i][j] = pow(vH(i,j), etapow);
    }
  }
  ierr = vWork2d[4].end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  // now use Mahaffy on eta to get grad h on staggered;
  // note   grad h = (3/8) eta^{-5/8} grad eta + grad b  because  h = H + b
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  ierr = vWork2d[4].get_array(eta); CHKERRQ(ierr);
  for (PetscInt o=0; o<2; o++) {

    GHOSTS = 1;
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
          h_x[o][i][j] += vbed.diff_x_stagE(i,j);
          h_y[o][i][j] += vbed.diff_y_stagE(i,j);
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
          h_y[o][i][j] += vbed.diff_y_stagN(i,j);
          h_x[o][i][j] += vbed.diff_x_stagN(i,j);
        }
      }
    }
  }
  ierr = vWork2d[4].end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);

  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief This method pre-computes delta, the quantity used later in SIA computations.
/*!
 * The computation result is on the staggered grid.
 */
PetscErrorCode IceModel::compute_delta() {
  PetscErrorCode ierr;

  PetscScalar *delta;

  double enhancement_factor = config.get("enhancement_factor"),
               constant_grain_size = config.get("constant_grain_size");

  // put "theta" from Schoof (2003) bed smoothness calculation in vWork2d[4]
  //   and the thickness relative to the smoothed bed in vWork2d[5]; each of
  //   vh,vH,vWork2d[4,5],sia_bed_smoother->topgsmooth must have stencil width
  //   WIDE_GHOSTS for this too work
  const PetscInt WIDE_GHOSTS = 2;
  ierr = sia_bed_smoother->get_theta(
    vh, config.get("Glen_exponent"), WIDE_GHOSTS, &vWork2d[4]); CHKERRQ(ierr);
  ierr = sia_bed_smoother->get_smoothed_thk(
    vh, vH, WIDE_GHOSTS, &vWork2d[5]); CHKERRQ(ierr);

  // surface gradient
  PetscScalar **h_x[2], **h_y[2];
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);

  PetscScalar **thk_smoothed, **theta_array;  
  ierr = vWork2d[4].get_array(theta_array); CHKERRQ(ierr);
  ierr = vWork2d[5].get_array(thk_smoothed); CHKERRQ(ierr);

  // delta:
  ierr = DeltaStag3[0].begin_access(); CHKERRQ(ierr);
  ierr = DeltaStag3[1].begin_access(); CHKERRQ(ierr);

  // some flow laws use grainsize, and even need age to update grainsize
  if ((realAgeForGrainSize==PETSC_TRUE) && (!config.get_flag("do_age"))) {
    PetscPrintf(grid.com,
       "PISM ERROR in IceModel::velocitySIAStaggered(): do_age not set but\n"
       "age is needed for grain-size-based flow law ...  ENDING! ...\n\n");
    PetscEnd();
  }
  const bool usetau3 =    (IceFlowLawUsesGrainSize(ice)
                       && (realAgeForGrainSize==PETSC_TRUE)
                       && (config.get_flag("do_age"))),
             usesGrainSize = IceFlowLawUsesGrainSize(ice);
  // age:
  PetscScalar *ageij, *ageoffset;
  if (usetau3) {
    ierr = tau3.begin_access(); CHKERRQ(ierr);
  }

  // enthalpy:
  bool e_age_coupling = config.get_flag("do_e_age_coupling");
  PetscScalar *Enthij, *Enthoffset;
  ierr = Enth3.begin_access(); CHKERRQ(ierr);

  // staggered grid computation of delta:
  for (PetscInt o=0; o<2; o++) {
    PetscInt GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        // staggered point: o=0 is i+1/2, o=1 is j+1/2,
        //   (i,j) and (i+oi,j+oj) are regular grid neighbors of staggered pt:
        const PetscInt     oi = 1-o, oj=o;  

        const PetscScalar
          thickness = 0.5 * ( thk_smoothed[i][j] + thk_smoothed[i+oi][j+oj] );

        if (thickness > 0) { 
          if (usetau3) {
            ierr = tau3.getInternalColumn(i,j,&ageij); CHKERRQ(ierr);
            ierr = tau3.getInternalColumn(i+oi,j+oj,&ageoffset); CHKERRQ(ierr);
          }
	  
          ierr = Enth3.getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
          ierr = Enth3.getInternalColumn(i+oi,j+oj,&Enthoffset); CHKERRQ(ierr);

          ierr = DeltaStag3[o].getInternalColumn(i,j,&delta); CHKERRQ(ierr);

          // does validity check for thickness:
          const PetscInt      ks = grid.kBelowHeight(thickness);  
          const PetscScalar   alpha =
                  sqrt(PetscSqr(h_x[o][i][j]) + PetscSqr(h_y[o][i][j]));
          const PetscReal     theta = 0.5 * ( theta_array[i][j] + theta_array[i+oi][j+oj] );

          for (PetscInt k=0; k<=ks; ++k) {
            PetscReal depth = thickness - grid.zlevels[k];
            // pressure added by the ice (i.e. pressure difference between the
            // current level and the top of the column)
            const PetscScalar   pressure = EC->getPressureFromDepth(depth);

            PetscScalar flow, grainsize = constant_grain_size;
            if (usetau3 && usesGrainSize && realAgeForGrainSize) {
              grainsize = grainSizeVostok(0.5 * (ageij[k] + ageoffset[k]));
            }
          
            if (e_age_coupling) {
              PetscReal age = 0.5 * (ageij[k] + ageoffset[k]) / secpera; // age in years
              if (age < 11000) {
                enhancement_factor = 1;
              } else {
                enhancement_factor = 3; 
              }
            } 
            // If the flow law does not use grain size, it will just ignore it, no harm there

            PetscScalar E = 0.5 * (Enthij[k] + Enthoffset[k]);
            flow = ice->flow_from_enth(alpha * pressure, E, pressure, grainsize);

            delta[k] = enhancement_factor * theta * 2.0 * pressure * flow;
          }

          for (PetscInt k=ks+1; k<grid.Mz; ++k) { // above the ice
            delta[k] = 0.0;
          }  

        } else {  // zero thickness case
          ierr = DeltaStag3[o].setColumn(i,j,0.0); CHKERRQ(ierr); 
        } 
      } // o
    } // j
  } // i

  // delta:
  ierr = DeltaStag3[1].end_access(); CHKERRQ(ierr);
  ierr = DeltaStag3[0].end_access(); CHKERRQ(ierr);

  // surface gradient:
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);

  // theta_array and thk_smoothed
  ierr = vWork2d[4].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[5].end_access(); CHKERRQ(ierr);

  // age
  if (usetau3) {
    ierr = tau3.end_access(); CHKERRQ(ierr);
  }

  // enthalpy
  ierr = Enth3.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes SIA deformational flow diffusivity on the staggered grid.
//! Uses DeltaStag3.
PetscErrorCode IceModel::compute_diffusivity_staggered() {
  PetscErrorCode ierr;

  PetscScalar *delta;

  // put the ice thickness relative to the smoothed bed in vWork2d[5]. See
  // [\ref Schoofbasaltopg2003].
  const PetscInt WIDE_GHOSTS = 2;
  ierr = sia_bed_smoother->get_smoothed_thk(vh, vH, WIDE_GHOSTS,
                                            &vWork2d[5]); CHKERRQ(ierr);

  ierr = vWork2d[5].begin_access(); CHKERRQ(ierr);
  ierr = vWork2dStag.begin_access(); CHKERRQ(ierr);

  ierr = DeltaStag3[0].begin_access(); CHKERRQ(ierr);
  ierr = DeltaStag3[1].begin_access(); CHKERRQ(ierr);

  PetscScalar Dmax = 0.0;
  // staggered grid computation of diffusivity
  for (PetscInt o=0; o<2; o++) {
    PetscInt GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        // staggered point: o=0 is i+1/2, o=1 is j+1/2,
        //   (i,j) and (i+oi,j+oj) are reg grid neighbors of staggered pt:
        const PetscInt     oi = 1-o, oj=o;  

        ierr = DeltaStag3[o].getInternalColumn(i,j,&delta); CHKERRQ(ierr);

        const PetscScalar
          thickness = 0.5 * ( vWork2d[5](i,j) + vWork2d[5](i+oi,j+oj) );

        if (thickness > 0) { 
          // does validity check for thickness:
          const PetscInt ks = grid.kBelowHeight(thickness);  

          PetscScalar Doffset = 0.0;  // diffusivity for deformational SIA flow
          for (PetscInt k = 1; k <= ks; ++k) {
            PetscReal depth = thickness - grid.zlevels[k];
              const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
              Doffset += 0.5 * dz * ((depth + dz) * delta[k-1] + depth * delta[k]);
          }

          // finish off D with (1/2) dz (0 + (H-z[ks])*delta[ks]), but dz=H-z[ks]:
          const PetscScalar dz = thickness - grid.zlevels[ks];
          Doffset += 0.5 * dz * dz * delta[ks];

          vWork2dStag(i,j,o) = Doffset;
          Dmax = PetscMax(Dmax, Doffset);
        } else {  // zero thickness case

          vWork2dStag(i,j,o) = 0.0;  // zero diffusivity

        } 
      } // o
    } // j
  } // i

  ierr = vWork2dStag.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[5].end_access(); CHKERRQ(ierr);

  ierr = DeltaStag3[0].end_access(); CHKERRQ(ierr);
  ierr = DeltaStag3[1].end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&Dmax, &gDmax, grid.com); CHKERRQ(ierr);
  
  return 0;
}

//! \brief Computes SIA deformational strain heating.
PetscErrorCode IceModel::compute_sigma() {
  PetscErrorCode ierr;
  PetscScalar *Sigma, *delta;

  // put the ice thickness relative to the smoothed bed in vWork2d[5]. See
  // [\ref Schoofbasaltopg2003].
  const PetscInt WIDE_GHOSTS = 2;
  ierr = sia_bed_smoother->get_smoothed_thk(vh, vH, WIDE_GHOSTS,
                                            &vWork2d[5]); CHKERRQ(ierr);

  ierr = DeltaStag3[0].begin_access(); CHKERRQ(ierr);
  ierr = DeltaStag3[1].begin_access(); CHKERRQ(ierr);

  ierr = Sigmastag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].begin_access(); CHKERRQ(ierr);

  // surface gradient
  PetscScalar **h_x[2], **h_y[2];
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);

  ierr = vWork2d[5].begin_access(); CHKERRQ(ierr);

  for (PetscInt o = 0; o < 2; ++o) {
    PetscInt GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        const PetscInt oi = 1-o, oj=o;  

        ierr = DeltaStag3[o].getInternalColumn(i,j,&delta); CHKERRQ(ierr);
        ierr = Sigmastag3[o].getInternalColumn(i,j,&Sigma); CHKERRQ(ierr);

        const PetscScalar
          thickness = 0.5 * ( vWork2d[5](i,j) + vWork2d[5](i+oi,j+oj) );

        const PetscInt ks = grid.kBelowHeight(thickness);  

        const PetscScalar alpha_squared =
          PetscSqr(h_x[o][i][j]) + PetscSqr(h_y[o][i][j]);

        for (PetscInt k=0; k<=ks; ++k) {
          PetscReal depth = thickness - grid.zlevels[k];
          PetscReal pressure = EC->getPressureFromDepth(depth);

          Sigma[k] = delta[k] * alpha_squared * pressure;
        }

        for (PetscInt k=ks+1; k<grid.Mz; ++k) { // above the ice
          Sigma[k] = 0.0;
        }  
        
      }
    }
  }

  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  
  ierr = DeltaStag3[1].end_access(); CHKERRQ(ierr);
  ierr = DeltaStag3[0].end_access(); CHKERRQ(ierr);

  // Now transfer Sigma from the staggered onto the regular grid.
  PetscScalar *Sigmareg, *SigmaEAST, *SigmaWEST, *SigmaNORTH, *SigmaSOUTH;
  ierr = Sigma3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscReal thickness = vWork2d[5](i,j);
      if (thickness > 0.0) {
        // horizontally average Sigma onto regular grid
        const PetscInt ks = grid.kBelowHeight(thickness);
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
      } else { // zero thickness case
        ierr = Sigma3.setColumn(i,j,0.0); CHKERRQ(ierr);
      }
    }
  }
  ierr = Sigma3.end_access(); CHKERRQ(ierr);

  ierr = Sigmastag3[1].end_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].end_access(); CHKERRQ(ierr);

  ierr = vWork2d[5].end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::compute_I_staggered() {
  PetscErrorCode ierr;
  PetscScalar *I, *delta;

  // put the ice thickness relative to the smoothed bed in vWork2d[5]. See
  // [\ref Schoofbasaltopg2003].
  const PetscInt WIDE_GHOSTS = 2;
  ierr = sia_bed_smoother->get_smoothed_thk(vh, vH, WIDE_GHOSTS,
                                            &vWork2d[5]); CHKERRQ(ierr);

  // delta
  ierr = DeltaStag3[0].begin_access(); CHKERRQ(ierr);
  ierr = DeltaStag3[1].begin_access(); CHKERRQ(ierr);
  // I (result)
  ierr = Istag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Istag3[1].begin_access(); CHKERRQ(ierr);
  // smoothed thk
  ierr = vWork2d[5].begin_access(); CHKERRQ(ierr);

  for (PetscInt o = 0; o < 2; ++o) {
    PetscInt GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        const PetscInt oi = 1-o, oj=o;  
        const PetscScalar
          thickness = 0.5 * ( vWork2d[5](i,j) + vWork2d[5](i+oi,j+oj) );

        ierr = DeltaStag3[o].getInternalColumn(i,j,&delta); CHKERRQ(ierr);
        ierr = Istag3[o].getInternalColumn(i,j,&I); CHKERRQ(ierr);

        const PetscInt ks = grid.kBelowHeight(thickness);  

        I[0] = 0.0;
        for (int k = 1; k < ks; ++k) {
          const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
          I[k] = I[k-1] + 0.5 * dz * (delta[k-1] + delta[k]);
        }

        for (PetscInt k=ks + 1; k < grid.Mz; ++k) { // above the ice
          I[k] = I[ks];
        }  
        
      }
    }
  }
  
  ierr = vWork2d[5].end_access(); CHKERRQ(ierr);

  ierr = Istag3[1].end_access(); CHKERRQ(ierr);
  ierr = Istag3[0].end_access(); CHKERRQ(ierr);

  ierr = DeltaStag3[1].end_access(); CHKERRQ(ierr);
  ierr = DeltaStag3[0].end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes vertically-integrated horizontal SIA velocities.
PetscErrorCode IceModel::velocity_SIA_2D() {
  PetscErrorCode ierr;

  const PetscInt WIDE_GHOSTS = 2;
  ierr = sia_bed_smoother->get_smoothed_thk(vh, vH, WIDE_GHOSTS,
                                            &vWork2d[5]); CHKERRQ(ierr);

  PetscScalar **h_x[2], **h_y[2];  
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);
  // thk_smoothed
  ierr = vWork2d[5].begin_access(); CHKERRQ(ierr);
  // diffusivity (staggered)
  ierr = vWork2dStag.begin_access(); CHKERRQ(ierr);

  ierr = uvbar.begin_access(); CHKERRQ(ierr);
  for (PetscInt o=0; o<2; o++) {
    PetscInt GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        // staggered point: o=0 is i+1/2, o=1 is j+1/2,
        //   (i,j) and (i+oi,j+oj) are reg grid neighbors of staggered pt:
        const PetscInt     oi = 1-o, oj=o;  
        const PetscScalar  slope = (o==0) ? h_x[o][i][j] : h_y[o][i][j];

        const PetscScalar
          thickness = 0.5 * ( vWork2d[5](i,j) + vWork2d[5](i+oi,j+oj) ),
          diffusivity = vWork2dStag(i,j,o);

        if (thickness > 0) { 
          // vertically-averaged SIA-only velocity, sans sliding; note
          //   uvbar(i,j,0) is  u  at E (east)  staggered point (i+1/2,j)
          //   uvbar(i,j,1) is  v  at N (north) staggered point (i,j+1/2)
          uvbar(i,j,o) = - diffusivity * slope / thickness;
        } else {
          uvbar(i,j,o) = 0.0;
        }
      } // j
    }   // i
  }     // o
  
  ierr = uvbar.end_access(); CHKERRQ(ierr);  

  return 0;
}

//! \brief Computes horizontal components of the SIA velocity within the ice.
PetscErrorCode IceModel::velocity_SIA_3D() {
  
  return 0;
}

//!  Compute the vertically-averaged horizontal velocity according to the non-sliding SIA.
/*!
See the comment for massContExplicitStep() before reading the rest of this comment.

In shallow ice approximation (SIA) areas, one may write either of two forms for
the vertically-integrated mass flux:
  \f[ \mathbf{q} = \bar{\mathbf{U}} H = - D \nabla h + \mathbf{U}_b \cdot H.\f]
Here \f$h\f$ is the surface elevation of the ice
\f$\mathbf{U}_b\f$ is the basal sliding velocity, and \f$D\f$ is the diffusivity,
which is computed in this method.

At the end of this routine the value of \f$D\f$ and of the <em>deformational
part of</em>  the vertically-averaged horizontal velocity \f$\bar{\mathbf{U}}\f$,
namely \f$- D \nabla h / H\f$, are known at all staggered grid points.  The
latter is stored in an IceModelVec2Stag called \c uvbar.

Two vertical integrals are computed.  Both are in terms of this internal quantity,
	\f[\delta(z) = e\,2\rho g (H-z) \,F,\f]
where \f$F\f$ is the result of the flow law, which depends on pressure and stress;
see IceFlowLaw::flow_from_enth().

One integral is evaluated at every level in the ice,
	\f[I(z) = \int_0^z \delta(z')\,dz'.\f]
These values are then stored in \c IceModelVec3 \c Istag3[2] and used later by
IceModel::horizontalVelocitySIARegular() if the horizontal velocity is needed
"at depth", that is, within the three dimensional ice fluid.

The other integral is the diffusivity of the SIA, used to compute
\c IceModelVec2Stag \c uvbar:
	\f[D = \int_0^H (H-z)\delta(z)\,dz.\f]
See [\ref BBL] on the meaning of \f$D\f$.  This method puts the maximum, over
all staggered points, of the diffusivities into the global variable gDmax.
Because the mass continuity PDE is actually nonlinear, this diffusivity \f$D\f$
changes at every time step.

Both integrals are approximated by the trapezoid rule.

The method surfaceGradientSIA() should be called before this one.  The
staggered-grid values of the surface gradient  \f$\nabla h\f$ are assumed to
reside in vWork2d[0],..,vWork2d[3].  In any case, all schemes here are essentially
"type I" in the classification described in (Hindmarsh and Payne 1996).
Note that \f$\nabla h\f$ is needed on the staggered grid although the surface 
elevation \f$h\f$ and thickness \f$H\f$ are known on the regular grid.  

This routine also computes the strain-heating \f$\Sigma\f$ in the SIA model.
This means the deformational-heating, not the basal friction.  The staggered
grid value of \f$\Sigma\f$ is computed using the formula appropriate
to the SIA, and is put in \c IceModelVec3 \c Sigmastag3[2].  See correctSigma()
for how the SIA+SSA hybrid \f$\Sigma\f$ is computed.
 */
PetscErrorCode IceModel::velocitySIAStaggered() {
  PetscErrorCode  ierr;

  PetscScalar *delta, *I, *Sigma;
  delta = new PetscScalar[grid.Mz];
  I = new PetscScalar[grid.Mz];
  Sigma = new PetscScalar[grid.Mz];

  const double enhancement_factor = config.get("enhancement_factor"),
               constant_grain_size = config.get("constant_grain_size");

  // put "theta" from Schoof (2003) bed smoothness calculation in vWork2d[4]
  //   and the thickness relative to the smoothed bed in vWork2d[5]; each of
  //   vh,vH,vWork2d[4,5],sia_bed_smoother->topgsmooth must have stencil width
  //   WIDE_GHOSTS for this too work
  const PetscInt WIDE_GHOSTS = 2;
  ierr = sia_bed_smoother->get_theta(
    vh, config.get("Glen_exponent"), WIDE_GHOSTS, &vWork2d[4]); CHKERRQ(ierr);
  ierr = sia_bed_smoother->get_smoothed_thk(
    vh, vH, WIDE_GHOSTS, &vWork2d[5]); CHKERRQ(ierr);

  PetscScalar **h_x[2], **h_y[2];  
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);
  ierr = vWork2d[4].begin_access(); CHKERRQ(ierr);
  ierr = vWork2d[5].begin_access(); CHKERRQ(ierr);
  ierr = uvbar.begin_access(); CHKERRQ(ierr);
  ierr = vWork2dStag.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = Istag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Istag3[1].begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].begin_access(); CHKERRQ(ierr);

  // some flow laws use grainsize, and even need age to update grainsize
  if ((realAgeForGrainSize==PETSC_TRUE) && (!config.get_flag("do_age"))) {
    PetscPrintf(grid.com,
       "PISM ERROR in IceModel::velocitySIAStaggered(): do_age not set but\n"
       "age is needed for grain-size-based flow law ...  ENDING! ...\n\n");
    PetscEnd();
  }
  const bool usetau3 =    (IceFlowLawUsesGrainSize(ice)
                       && (realAgeForGrainSize==PETSC_TRUE)
                       && (config.get_flag("do_age"))),
             usesGrainSize = IceFlowLawUsesGrainSize(ice);
  PetscScalar *ageij, *ageoffset;
  if (usetau3) {
    ierr = tau3.begin_access(); CHKERRQ(ierr);
  }

  // some flow laws use enthalpy while some ("cold ice methods") use temperature
  PetscScalar *Enthij, *Enthoffset;

  ierr = Enth3.begin_access(); CHKERRQ(ierr);

  PetscScalar Dmax = 0.0;
  // staggered grid computation of: uvbar, I, Sigma
  for (PetscInt o=0; o<2; o++) {
    PetscInt GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        // staggered point: o=0 is i+1/2, o=1 is j+1/2,
        //   (i,j) and (i+oi,j+oj) are reg grid neighbors of staggered pt:
        const PetscInt     oi = 1-o, oj=o;  
        const PetscScalar  slope = (o==0) ? h_x[o][i][j] : h_y[o][i][j];

        const PetscScalar
          thickness = 0.5 * ( vWork2d[5](i,j) + vWork2d[5](i+oi,j+oj) );

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
          const PetscReal     theta = 0.5 * ( vWork2d[4](i,j) + vWork2d[4](i+oi,j+oj) );

          I[0] = 0;
          PetscScalar  Dfoffset = 0.0;  // diffusivity for deformational SIA flow
          for (PetscInt k=0; k<=ks; ++k) {
            PetscReal depth = thickness - grid.zlevels[k];
            // pressure added by the ice (i.e. pressure difference between the
            // current level and the top of the column)
            const PetscScalar   pressure = ice->rho * standard_gravity * depth;
            PetscScalar flow, grainsize = constant_grain_size;
            if (usetau3 && usesGrainSize && realAgeForGrainSize) {
              grainsize = grainSizeVostok(0.5 * (ageij[k] + ageoffset[k]));
            }
            // If the flow law does not use grain size, it will just ignore it, no harm there
            PetscScalar E = 0.5 * (Enthij[k] + Enthoffset[k]);
            flow = ice->flow_from_enth(alpha * pressure, E, pressure, grainsize);

            delta[k] = enhancement_factor * theta * 2.0 * pressure * flow;

            // for Sigma, ignore mask value and assume SHEET; will be overwritten
            // by correctSigma() in iMssa.cc
            Sigma[k] = delta[k] * PetscSqr(alpha) * pressure;

            if (k>0) { // trapezoid rule for I[k] and K[k]
              const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
              I[k] = I[k-1] + 0.5 * dz * (delta[k-1] + delta[k]);
              Dfoffset += 0.5 * dz * ((depth + dz) * delta[k-1] + depth * delta[k]);
            }
          }
          // finish off D with (1/2) dz (0 + (H-z[ks])*delta[ks]), but dz=H-z[ks]:
          const PetscScalar dz = thickness - grid.zlevels[ks];
          Dfoffset += 0.5 * dz * dz * delta[ks];

          vWork2dStag(i,j,o) = Dfoffset;
          Dmax = PetscMax(Dmax, Dfoffset);

          for (PetscInt k=ks+1; k<grid.Mz; ++k) { // above the ice
            Sigma[k] = 0.0;
            I[k] = I[ks];
          }  
          ierr = Istag3[o].setInternalColumn(i, j, I); CHKERRQ(ierr);
          ierr = Sigmastag3[o].setInternalColumn(i, j, Sigma); CHKERRQ(ierr);

          // vertically-averaged SIA-only velocity, sans sliding; note
          //   uvbar(i,j,0) is  u  at E (east)  staggered point (i+1/2,j)
          //   uvbar(i,j,1) is  v  at N (north) staggered point (i,j+1/2)
          uvbar(i,j,o) = - Dfoffset * slope / thickness;
         
        } else {  // zero thickness case
          uvbar(i,j,o) = 0.0;
          vWork2dStag(i,j,o) = 0.0;  // zero diffusivity
          ierr = Istag3[o].setColumn(i,j,0.0); CHKERRQ(ierr);
          ierr = Sigmastag3[o].setColumn(i,j,0.0); CHKERRQ(ierr);
        } 
      } // o
    } // j
  } // i

  ierr = uvbar.end_access(); CHKERRQ(ierr);
  ierr = vWork2dStag.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[4].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[5].end_access(); CHKERRQ(ierr);
  if (usetau3) {
    ierr = tau3.end_access(); CHKERRQ(ierr);
  }
  ierr = w3.end_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].end_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].end_access(); CHKERRQ(ierr);
  ierr = Istag3[0].end_access(); CHKERRQ(ierr);
  ierr = Istag3[1].end_access(); CHKERRQ(ierr);

  ierr = Enth3.end_access(); CHKERRQ(ierr);

  delete [] delta;   delete [] I;   delete [] Sigma;

  ierr = PetscGlobalMax(&Dmax, &gDmax, grid.com); CHKERRQ(ierr);
  
  return 0;
}


//! Compute the basal sliding and frictional heating if (where) SIA sliding rule is used.
/*!
This routine is only called, by velocity(), if \f$\mu\f$=\c mu_sliding is
non-zero.

THIS KIND OF SIA SLIDING LAW IS A BAD IDEA.  THAT'S WHY \f$\mu\f$ IS SET TO 
ZERO BY DEFAULT.  See Appendix B of [\ref BBssasliding] for the dangers in this
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
  PetscScalar **h_x[2], **h_y[2], **Rb;

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
  ierr = vh.begin_access(); CHKERRQ(ierr);
  ierr = sia_bed_smoother->topgsmooth.begin_access(); CHKERRQ(ierr);
  
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

          // change r1200: new meaning of H
          const PetscScalar H = vh(i,j) - sia_bed_smoother->topgsmooth(i,j);

          ierr = EC->getAbsTemp(Enth3.getValZ(i,j,0.0),
                                EC->getPressureFromDepth(H), T); CHKERRQ(ierr);

          basalC = basalVelocitySIA(myx, myy, H, T, 
				    alpha, mu_sliding,
				    minimum_temperature_for_sliding);
          bvel[i][j].u = - basalC * myhx;
          bvel[i][j].v = - basalC * myhy;
          // basal frictional heating; note P * dh/dx is x comp. of basal shear stress
          // in ice streams this result will be *overwritten* by
          //   correctBasalFrictionalHeating() if useSSAVelocities==TRUE
          const PetscScalar P = ice->rho * standard_gravity * H;
          Rb[i][j] = - (P * myhx) * bvel[i][j].u - (P * myhy) * bvel[i][j].v;
        }
      }
    }
  }
  
  ierr = vh.end_access(); CHKERRQ(ierr);
  ierr = sia_bed_smoother->topgsmooth.end_access(); CHKERRQ(ierr);
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

Only two-dimensional regular grid velocities are updated here. The full
three-dimensional velocity field is not updated here but instead in
horizontalVelocitySIARegular() and in vertVelocityFromIncompressibility().
 */
PetscErrorCode IceModel::velocities2DSIAToRegular() {  
  PetscErrorCode ierr;
  double mu_sliding = config.get("mu_sliding");

  ierr = uvbar.staggered_to_regular(vel_bar); CHKERRQ(ierr); 

  if (mu_sliding > 0.0) {
    ierr = vel_bar.begin_access(); CHKERRQ(ierr);
    ierr = vel_basal.begin_access(); CHKERRQ(ierr);
    for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
        vel_bar(i,j).u += vel_basal(i,j).u;
        vel_bar(i,j).v += vel_basal(i,j).v;
      }
    }
    ierr = vel_basal.end_access(); CHKERRQ(ierr);
    ierr = vel_bar.end_access(); CHKERRQ(ierr);
  }

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
  PetscScalar *Sigmareg, *SigmaEAST, *SigmaWEST, *SigmaNORTH, *SigmaSOUTH;

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = Sigma3.begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[0].begin_access(); CHKERRQ(ierr);
  ierr = Sigmastag3[1].begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vH(i,j) > 0.0) {
        // horizontally average Sigma onto regular grid
        const PetscInt ks = grid.kBelowHeight(vH(i,j));
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
quantities at depth, including \f$I(z)\f$.  This procedure takes \f$I(z)\f$ and
the surface slope, both known on the staggered grid, and computes
the three-dimensional arrays for the horizontal components \f$u\f$ and 
\f$v\f$ of the velocity field:
	\f[(u(z),v(z)) = - I(z) (h_x, h_y).\f]

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
[\ref BBssasliding] for the dangers in this mechanism.

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

