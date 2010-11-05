// Copyright (C) 2004--2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "SIAFD.hh"

//! \brief Compute the ice surface gradient for the SIA.
/*! 
  There are three methods for computing the surface gradient. Which method is
  controlled by configuration parameter \c surface_gradient_method which can
  have values \c haseloff, \c mahaffy, or \c eta.

  The most traditional method is to directly differentiate the surface
  elevation \f$h\f$ by the Mahaffy method \ref Mahaffy. The \c haseloff method,
  suggested by Marianne Haseloff, modifies the Mahaffy method only where
  ice-free adjacent bedrock points are above the ice surface, and in those
  cases the returned gradient component is zero.

  The alternative method, when \c surface_gradient_method = \c eta, transforms
  the thickness to something more regular and differentiates that. We get back
  to the gradient of the surface by applying the chain rule. In particular, as
  shown in \ref CDDSV for the flat bed and \f$n=3\f$ case, if we define

  \f[\eta = H^{(2n+2)/n}\f]

  then \f$\eta\f$ is more regular near the margin than \f$H\f$. So we compute
  the surface gradient by

  \f[\nabla h = \frac{n}{(2n+2)} \eta^{(-n-2)/(2n+2)} \nabla \eta + \nabla b,\f]

  recalling that \f$h = H + b\f$. This method is only applied when \f$\eta >
  0\f$ at a given point; otherwise \f$\nabla h = \nabla b\f$.

  In all cases we are computing the gradient by finite differences onto a
  staggered grid. In the method with \f$\eta\f$ we apply centered differences
  using (roughly) the same method for \f$\eta\f$ and \f$b\f$ that applies
  directly to the surface elevation \f$h\f$ in the \c mahaffy and \c haseloff
  methods.
*/
PetscErrorCode SIAFD::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode  ierr;

  const string method = config.get_string("surface_gradient_method");

  if ((method != "eta") && (method != "mahaffy") && (method != "haseloff")) {
    verbPrintf(1, grid.com,
      "PISM ERROR: value of surface_gradient_method, option -gradient, not valid ... ending\n");
    PetscEnd();
  }

  if (method == "eta") {

    ierr = surface_gradient_eta(h_x, h_y); CHKERRQ(ierr); 

  } else if (method == "haseloff") {

    ierr = surface_gradient_haseloff(h_x, h_y); CHKERRQ(ierr);

  } else if (method == "mahaffy") {

    ierr = surface_gradient_mahaffy(h_x, h_y); CHKERRQ(ierr);

  } else {
    SETERRQ(1, "can't happen");
  }

  return 0;
}


//! \brief Compute the ice surface gradient using the eta-transformation.
PetscErrorCode SIAFD::surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;

  const PetscScalar n = ice.exponent(), // presumably 3.0
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0);
  const PetscScalar dx = grid.dx, dy = grid.dy;  // convenience
  IceModelVec2S eta = tmp1;

  // compute eta = H^{8/3}, which is more regular, on reg grid
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = eta.begin_access(); CHKERRQ(ierr);
  PetscInt GHOSTS = 2;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      eta(i,j) = pow((*thickness)(i,j), etapow);
    }
  }
  ierr = eta.end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);
  
  // now use Mahaffy on eta to get grad h on staggered;
  // note   grad h = (3/8) eta^{-5/8} grad eta + grad b  because  h = H + b
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = eta.begin_access(); CHKERRQ(ierr);
  for (PetscInt o=0; o<2; o++) {

    GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        if (o==0) {     // If I-offset
          const PetscScalar mean_eta = 0.5 * (eta(i+1,j) + eta(i,j));
          if (mean_eta > 0.0) {
            const PetscScalar factor = invpow * pow(mean_eta, dinvpow);
            h_x(i,j,o) = factor * (eta(i+1,j) - eta(i,j)) / dx;
            h_y(i,j,o) = factor * (+ eta(i+1,j+1) + eta(i,j+1)
                                   - eta(i+1,j-1) - eta(i,j-1)) / (4.0*dy);
          } else {
            h_x(i,j,o) = 0.0;
            h_y(i,j,o) = 0.0;
          }
          // now add bed slope to get actual h_x,h_y
          h_x(i,j,o) += bed->diff_x_stagE(i,j);
          h_y(i,j,o) += bed->diff_y_stagE(i,j);
        } else {        // J-offset
          const PetscScalar mean_eta = 0.5 * (eta(i,j+1) + eta(i,j));
          if (mean_eta > 0.0) {
            const PetscScalar factor = invpow * pow(mean_eta, dinvpow);
            h_y(i,j,o) = factor * (eta(i,j+1) - eta(i,j)) / dy;
            h_x(i,j,o) = factor * (+ eta(i+1,j+1) + eta(i+1,j)
                                   - eta(i-1,j+1) - eta(i-1,j)) / (4.0*dx);
          } else {
            h_y(i,j,o) = 0.0;
            h_x(i,j,o) = 0.0;
          }
          // now add bed slope to get actual h_x,h_y
          h_y(i,j,o) += bed->diff_y_stagN(i,j);
          h_x(i,j,o) += bed->diff_x_stagN(i,j);
        }
      }
    }
  }
  ierr = eta.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);

  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the ice surface gradient using Marianne Haseloff's approach.
/*!
 * Deals correctly with adjacent ice-free points with bed elevations which are
 * above the surface of the ice
 */
PetscErrorCode SIAFD::surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;

  const PetscScalar Hicefree = 0.0;  // standard for ice-free, in Haseloff
  const PetscScalar dx = grid.dx, dy = grid.dy;  // convenience

  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);

  PetscScalar **h, **b, **H;
  ierr = bed->get_array(b); CHKERRQ(ierr);
  ierr = thickness->get_array(H); CHKERRQ(ierr);
  ierr = surface->get_array(h); CHKERRQ(ierr);
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
          h_x(i,j,o) = (hhE - h[i][j]) / dx;

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
          h_y(i,j,o) = (hhNE + hhN - hhSE - hhS) / (4.0 * dy);
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
          h_y(i,j,o) = (hhN - h[i][j]) / dy;

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
          h_x(i,j,o) = (hhNE + hhE - hhNW - hhW) / (4.0 * dx);
        }

      } // j
    }   // i
  }     // o
  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr =       bed->end_access(); CHKERRQ(ierr);
  ierr =   surface->end_access(); CHKERRQ(ierr);

  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the ice surface gradient using the Mary Anne Mahaffy method;
//! see [\ref Mahaffy].
PetscErrorCode SIAFD::surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;
  const PetscScalar dx = grid.dx, dy = grid.dy;  // convenience

  PetscScalar **h;
  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);
  ierr = surface->get_array(h); CHKERRQ(ierr);

  for (PetscInt o=0; o<2; o++) {
    PetscInt GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        if (o==0) {     // If I-offset
          h_x(i,j,o) = (h[i+1][j] - h[i][j]) / dx;
          h_y(i,j,o) = (+ h[i+1][j+1] + h[i][j+1]
                        - h[i+1][j-1] - h[i][j-1]) / (4.0*dy);
        } else {        // J-offset
          h_y(i,j,o) = (h[i][j+1] - h[i][j]) / dy;
          h_x(i,j,o) = (+ h[i+1][j+1] + h[i+1][j]
                        - h[i-1][j+1] - h[i-1][j]) / (4.0*dx);
        }
      }
    }
  }

  ierr = surface->end_access(); CHKERRQ(ierr);
  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the SIA flux. If fast == false, also store delta on the staggered grid.
PetscErrorCode SIAFD::compute_diffusive_flux(IceModelVec2Stag &result, bool fast) {
  PetscErrorCode  ierr;
  IceModelVec2S thk_smooth = tmp1,
    theta = tmp2;
  IceModelVec2Stag h_x = tmp3, h_y = tmp4;

  PetscScalar *delta;

  delta = new PetscScalar[grid.Mz];

  if (!fast) {
    ierr = delta_staggered[0].begin_access(); CHKERRQ(ierr);
    ierr = delta_staggered[1].begin_access(); CHKERRQ(ierr);
  }

  const double enhancement_factor = config.get("enhancement_factor"),
    constant_grain_size = config.get("constant_grain_size"),
    standard_gravity = config.get("standard_gravity");

  bool compute_grain_size_using_age = config.get_flag("compute_grain_size_using_age");

  const bool use_age = (IceFlowLawUsesGrainSize(&ice) &&
                        compute_grain_size_using_age &&
                        config.get_flag("do_age"));

  // get "theta" from Schoof (2003) bed smoothness calculation and the
  //   thickness relative to the smoothed bed; each IceModelVec2S involved must
  //   have stencil width WIDE_GHOSTS for this too work
  const PetscInt WIDE_GHOSTS = 2;
  ierr = sia_bed_smoother->get_theta(*surface, config.get("Glen_exponent"),
                                     WIDE_GHOSTS, &theta); CHKERRQ(ierr);

  ierr = sia_bed_smoother->get_smoothed_thk(*surface, *thickness, WIDE_GHOSTS,
                                            &thk_smooth); CHKERRQ(ierr);

  ierr = theta.begin_access(); CHKERRQ(ierr);
  ierr = thk_smooth.begin_access(); CHKERRQ(ierr);
  ierr = diffusive_flux.begin_access(); CHKERRQ(ierr);

  // some flow laws use grainsize, and even need age to update grainsize
  if (compute_grain_size_using_age && (!config.get_flag("do_age"))) {
    PetscPrintf(grid.com,
                "PISM ERROR in IceModel::velocitySIAStaggered(): do_age not set but\n"
                "age is needed for grain-size-based flow law ...  ENDING! ...\n\n");
    PetscEnd();
  }

  PetscScalar *ageij, *ageoffset;
  if (use_age) {
    ierr = age->begin_access(); CHKERRQ(ierr);
  }

  // some flow laws use enthalpy while some ("cold ice methods") use temperature
  PetscScalar *Enthij, *Enthoffset;
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);

  PetscScalar my_D_max = 0.0;
  // staggered grid computation of: diffusive_flux, I, Sigma
  for (PetscInt o=0; o<2; o++) {
    PetscInt GHOSTS = 1;
    for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        // staggered point: o=0 is i+1/2, o=1 is j+1/2, (i,j) and (i+oi,j+oj)
        //   are regular grid neighbors of a staggered point:
        const PetscInt oi = 1 - o, oj = o;  

        const PetscScalar
          thickness = 0.5 * ( thk_smooth(i,j) + thk_smooth(i+oi,j+oj) );

        // zero thickness case:
        if (thickness == 0.0) {
          diffusive_flux(i,j,o) = 0.0;
          if (!fast) {
            ierr = delta_staggered[o].setColumn(i, j, 0.0); CHKERRQ(ierr);
          }
          continue;
        }

        if (use_age) {
          ierr = age->getInternalColumn(i, j, &ageij); CHKERRQ(ierr);
          ierr = age->getInternalColumn(i+oi, j+oj, &ageoffset); CHKERRQ(ierr);
        }
	  
        ierr = enthalpy->getInternalColumn(i, j, &Enthij); CHKERRQ(ierr);
        ierr = enthalpy->getInternalColumn(i+oi, j+oj, &Enthoffset); CHKERRQ(ierr);

        const PetscScalar slope = (o==0) ? h_x(i,j,o) : h_y(i,j,o);
        const PetscInt      ks = grid.kBelowHeight(thickness);  
        const PetscScalar   alpha =
          sqrt(PetscSqr(h_x(i,j,o)) + PetscSqr(h_y(i,j,o)));
        const PetscReal theta_local = 0.5 * ( theta(i,j) + theta(i+oi,j+oj) );

        PetscScalar  Dfoffset = 0.0;  // diffusivity for deformational SIA flow
        for (PetscInt k = 0; k <= ks; ++k) {
          PetscReal depth = thickness - grid.zlevels[k];
          // pressure added by the ice (i.e. pressure difference between the
          // current level and the top of the column)
          const PetscScalar   pressure = ice.rho * standard_gravity * depth;

          PetscScalar flow, grainsize = constant_grain_size;
          if (use_age) {
            grainsize = grainSizeVostok(0.5 * (ageij[k] + ageoffset[k]));
          }
          // If the flow law does not use grain size, it will just ignore it,
          // no harm there
          PetscScalar E = 0.5 * (Enthij[k] + Enthoffset[k]);
          flow = ice.flow_from_enth(alpha * pressure, E, pressure, grainsize);

          delta[k] = enhancement_factor * theta_local * 2.0 * pressure * flow;

          if (k > 0) { // trapezoid rule for I[k] and K[k]
            const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
            Dfoffset += 0.5 * dz * ((depth + dz) * delta[k-1] + depth * delta[k]);
          }
        }
        // finish off D with (1/2) dz (0 + (H-z[ks])*delta[ks]), but dz=H-z[ks]:
        const PetscScalar dz = thickness - grid.zlevels[ks];
        Dfoffset += 0.5 * dz * dz * delta[ks];

        my_D_max = PetscMax(my_D_max, Dfoffset);

        // vertically-averaged SIA-only flux, sans sliding; note
        //   diffusive_flux(i,j,0) is  u  at E (east)  staggered point (i+1/2,j)
        //   diffusive_flux(i,j,1) is  v  at N (north) staggered point (i,j+1/2)
        diffusive_flux(i,j,o) = - Dfoffset * slope;

        // if doing the full update, fill the delta column above the ice and
        // store it:
        if (!fast) {
          for (PetscInt k = ks + 1; k < grid.Mz; ++k) { // above the ice
            delta[k] = 0.0;
          }  
          ierr = delta_staggered[o].setInternalColumn(i,j,delta); CHKERRQ(ierr);
        }
      } // o
    } // j
  } // i

  ierr = diffusive_flux.end_access(); CHKERRQ(ierr);
  ierr = theta.end_access(); CHKERRQ(ierr);
  ierr = thk_smooth.end_access(); CHKERRQ(ierr);

  if (use_age) {
    ierr = age->end_access(); CHKERRQ(ierr);
  }

  ierr = enthalpy->end_access(); CHKERRQ(ierr);

  delete [] delta;

  if (!fast) {
    ierr = delta_staggered[1].end_access(); CHKERRQ(ierr);
    ierr = delta_staggered[0].end_access(); CHKERRQ(ierr);
  }

  ierr = PetscGlobalMax(&my_D_max, &D_max, grid.com); CHKERRQ(ierr);
  
  return 0;
}

//! Use the Vostok core as a source of a relationship between the age of the ice and the grain size.
/*! 
A data set is interpolated here.  The intention is that the softness of the ice has
nontrivial dependence on its age, through its grainsize, because of variable dustiness
of the global climate.  The grainsize is partly determined by at which point in 
the glacial cycle the given ice fell as snow.

The data is from \ref DeLaChapelleEtAl98 and 
\ref LipenkovEtAl89 .  In particular, Figure A2 in the former reference was
hand-sampled with an attempt to include the ``wiggles'' in that figure.  Ages of
the oldest ice (>= 300 ka) were estimated in a necessarily ad hoc way.  The 
age value of 10000 ka was added simply to give interpolation for very old ice;
ages beyond that get constant extrapolation.  Linear interpolation is done between
the samples.
 */
PetscScalar SIAFD::grainSizeVostok(PetscScalar age) const {
  const PetscInt numPoints = 22;
  const PetscScalar ageAt[numPoints] = {  // ages in ka
    0.0000e+00, 5.0000e+01, 1.0000e+02, 1.2500e+02, 1.5000e+02,
    1.5800e+02, 1.6500e+02, 1.7000e+02, 1.8000e+02, 1.8800e+02,
    2.0000e+02, 2.2500e+02, 2.4500e+02, 2.6000e+02, 3.0000e+02,
    3.2000e+02, 3.5000e+02, 4.0000e+02, 5.0000e+02, 6.0000e+02,
    8.0000e+02, 1.0000e+04 };
  const PetscScalar gsAt[numPoints] = {   // grain sizes in m
    1.8000e-03, 2.2000e-03, 3.0000e-03, 4.0000e-03, 4.3000e-03,
    3.0000e-03, 3.0000e-03, 4.6000e-03, 3.4000e-03, 3.3000e-03,
    5.9000e-03, 6.2000e-03, 5.4000e-03, 6.8000e-03, 3.5000e-03,
    6.0000e-03, 8.0000e-03, 8.3000e-03, 3.6000e-03, 3.8000e-03,
    9.5000e-03, 1.0000e-02 };
  const PetscScalar a = age * 1.0e-3 / secpera; // Age in ka
  PetscInt l = 0;               // Left end of the binary search
  PetscInt r = numPoints - 1;   // Right end

  // If we are out of range
  if (a < ageAt[l]) {
    return gsAt[l];
  } else if (a > ageAt[r]) {
    return gsAt[r];
  }
  // Binary search for the interval
  while (r > l + 1) {
    const PetscInt j = (r + l) / 2;
    if (a < ageAt[j]) {
      r = j;
    } else {
      l = j;
    }
  }
  if ((r == l) || (PetscAbsReal(r - l) > 1)) {
    PetscPrintf(grid.com, "binary search in grainSizeVostok: oops.\n");
  }
  // Linear interpolation on the interval
  return gsAt[l] + (a - ageAt[l]) * (gsAt[r] - gsAt[l]) / (ageAt[r] - ageAt[l]);
}
