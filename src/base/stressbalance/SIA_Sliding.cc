// Copyright (C) 2004--2012 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "SIA_Sliding.hh"
#include "Mask.hh"
#include "enthalpyConverter.hh"
#include "PISMVars.hh"
#include "flowlaw_factory.hh"

PetscErrorCode SIA_Sliding::allocate() {
  PetscErrorCode ierr;
  PetscInt WIDE_STENCIL = grid.max_stencil_width;

  for (int i = 0; i < 2; ++i) {
    char namestr[30];

    ierr = work_2d_stag[i].create(grid, "work_vector", true); CHKERRQ(ierr);
    snprintf(namestr, sizeof(namestr), "work_vector_2d_stag_%d", i);
    ierr = work_2d_stag[i].set_name(namestr); CHKERRQ(ierr);

  }

  ierr = work_2d.create(grid, "work_vector_2d", true, WIDE_STENCIL); CHKERRQ(ierr);

  {
    IceFlowLawFactory ice_factory(grid.com, "sia_", config, &EC);

    ierr = ice_factory.setType(config.get_string("sia_flow_law").c_str()); CHKERRQ(ierr);

    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ierr = ice_factory.create(&flow_law); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode SIA_Sliding::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = ShallowStressBalance::init(vars); CHKERRQ(ierr);

  standard_gravity = config.get("standard_gravity");
  verification_mode = config.get_flag("sia_sliding_verification_mode");

  if (config.has("EISMINT_II_experiment"))
    eisII_experiment = config.get_string("EISMINT_II_experiment");

  thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (surface == NULL) SETERRQ(grid.com, 1, "surface_altitude is not available");

  bed = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (bed == NULL) SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(grid.com, 1, "enthalpy is not available");

  return 0;
}

//! Compute the basal sliding and frictional heating if (where) SIA sliding rule is used.
/*!
  THIS KIND OF SIA SLIDING LAW IS A BAD IDEA. THAT'S WHY \f$\mu\f$ IS SET TO
  ZERO BY DEFAULT. See Appendix B of [\ref BBssasliding] for the dangers in
  this mechanism.

  This routine calls the SIA-type sliding law, which may return zero in the
  frozen base case; see basalVelocitySIA(). The basal sliding velocity is
  computed for all SIA points. This routine also computes the basal frictional
  heating.

  The strain heating contribution is ignored by this code.
 */
PetscErrorCode SIA_Sliding::update(bool /*fast*/) {
  PetscErrorCode ierr;
  IceModelVec2Stag h_x = work_2d_stag[0], h_y = work_2d_stag[1];

  ierr = compute_surface_gradient(h_x, h_y); CHKERRQ(ierr);

  ierr = D2.set(0.0); CHKERRQ(ierr);

  double mu_sliding = config.get("mu_sliding"),
    minimum_temperature_for_sliding = config.get("minimum_temperature_for_sliding"),
    ice_rho = config.get("ice_density");

  MaskQuery m(*mask);

  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);

  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = surface->begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = basal_frictional_heating.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      if (m.ocean(i,j)) {
        velocity(i,j).u = 0.0;
        velocity(i,j).v = 0.0;
        basal_frictional_heating(i,j) = 0.0;
      } else {
        // basal velocity from SIA-type sliding law: not recommended!
        const PetscScalar
          myx = grid.x[i],
          myy = grid.y[j],
          myhx = 0.25 * (  h_x(i,j,0) + h_x(i-1,j,0)
                           + h_x(i,j,1) + h_x(i,j-1,1)),
          myhy = 0.25 * (  h_y(i,j,0) + h_y(i-1,j,0)
                           + h_y(i,j,1) + h_y(i,j-1,1)),
          alpha = sqrt(PetscSqr(myhx) + PetscSqr(myhy));
        PetscScalar T, basalC;

        // change r1200: new meaning of H
        const PetscScalar H = (*surface)(i,j) - (*bed)(i,j);

        ierr = EC.getAbsTemp(enthalpy->getValZ(i,j,0.0),
                             EC.getPressureFromDepth(H), T); CHKERRQ(ierr);

        basalC = basalVelocitySIA(myx, myy, H, T,
                                  alpha, mu_sliding,
                                  minimum_temperature_for_sliding);
        velocity(i,j).u = - basalC * myhx;
        velocity(i,j).v = - basalC * myhy;
        // basal frictional heating; note P * dh/dx is x comp. of basal shear stress
        // in ice streams this result will be *overwritten* by
        //   correctBasalFrictionalHeating() if useSSAVelocities==TRUE
        const PetscScalar P = ice_rho * standard_gravity * H;
        basal_frictional_heating(i,j) = - (P * myhx) * velocity(i,j).u - (P * myhy) * velocity(i,j).v;
      }
    }
  }

  ierr = velocity.end_access(); CHKERRQ(ierr);
  ierr = basal_frictional_heating.end_access(); CHKERRQ(ierr);

  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  ierr = surface->end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = enthalpy->end_access(); CHKERRQ(ierr);

  ierr = velocity.beginGhostComm(); CHKERRQ(ierr);
  ierr = velocity.endGhostComm(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the coefficient of surface gradient, for basal sliding
//! velocity as a function of driving stress in SIA regions.
/*!
  THIS KIND OF SIA SLIDING LAW IS A BAD IDEA IN A THERMOMECHANICALLY-COUPLED
  MODEL.  THAT'S WHY \f$\mu\f$ IS SET TO ZERO BY DEFAULT.

  We allow the SIA sliding law of the form
  \f[ \mathbf{U}_b = (u_b,v_b) = - C \nabla h. \f]
  Here \f$\mathbf{U}_b\f$ is the horizontal velocity of the base of
  the ice (the "sliding velocity") and \f$h\f$ is the elevation of the ice
  surface.  This procedure returns the \em positive \em coefficient \f$C\f$ in
  this relationship.  This coefficient can depend of the thickness, the basal
  temperature, and the horizontal location.

  The default version for IceModel here is location-independent
  pressure-melting-temperature-activated linear sliding.  See Appendix B of
  [\ref BBssasliding] for the dangers in this mechanism.

  Parameter \f$\mu\f$ can be set by option \c -mu_sliding.

  The returned coefficient is used in update() (above).
*/
PetscScalar SIA_Sliding::basalVelocitySIA(PetscScalar xIN, PetscScalar yIN,
                                          PetscScalar H, PetscScalar T,
                                          PetscScalar /*alpha*/, PetscScalar mu,
                                          PetscScalar min_T) const {
  PetscReal ice_rho = config.get("ice_density"),
    beta_CC_grad = config.get("beta_CC") * ice_rho * config.get("standard_gravity");

  if (verification_mode) {
    // test 'E' mode
    const PetscScalar r1 = 200e3, r2 = 700e3,   /* define region of sliding */
      theta1 = 10 * (pi/180), theta2 = 40 * (pi/180);
    const PetscScalar x = fabs(xIN), y = fabs(yIN);
    const PetscScalar r = sqrt(x * x + y * y);
    PetscScalar       theta;
    if (x < 1.0)
      theta = pi / 2.0;
    else
      theta = atan(y / x);

    if ((r > r1) && (r < r2) && (theta > theta1) && (theta < theta2)) {
      // now INSIDE sliding region
      const PetscScalar rbot = (r2 - r1) * (r2 - r1),
        thetabot = (theta2 - theta1) * (theta2 - theta1);
      const PetscScalar mu_max = 2.5e-11; /* Pa^-1 m s^-1; max sliding coeff */
      PetscScalar muE = mu_max * (4.0 * (r - r1) * (r2 - r) / rbot)
        * (4.0 * (theta - theta1) * (theta2 - theta) / thetabot);
      return muE * ice_rho * standard_gravity * H;
    } else
      return 0.0;
  }

  if ((eisII_experiment == "G") || (eisII_experiment == "H")) {
    const PetscScalar  Bfactor = 1e-3 / secpera; // m s^-1 Pa^-1
    PetscReal pressure = EC.getPressureFromDepth(H), E;
    EC.getEnthPermissive(T, 0.0, pressure, E);

    if (eisII_experiment == "G") {
      return Bfactor * ice_rho * standard_gravity * H;
    } else if (eisII_experiment == "H") {
      if (EC.isTemperate(E, pressure)) {
        return Bfactor * ice_rho * standard_gravity * H; // ditto case G
      } else {
        return 0.0;
      }
    }
    return 0.0;  // zero sliding for other tests
  }

  // the "usual" case:
  if (T + beta_CC_grad * H > min_T) {
    const PetscScalar p_over = ice_rho * standard_gravity * H;
    return mu * p_over;
  } else {
    return 0;
  }
}

PetscErrorCode SIA_Sliding::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode  ierr;

  const string method = config.get_string("surface_gradient_method");

  if (method == "eta") {

    ierr = surface_gradient_eta(h_x, h_y); CHKERRQ(ierr);

  } else if (method == "haseloff") {

    ierr = surface_gradient_haseloff(h_x, h_y); CHKERRQ(ierr);

  } else if (method == "mahaffy") {

    ierr = surface_gradient_mahaffy(h_x, h_y); CHKERRQ(ierr);

  } else {
    verbPrintf(1, grid.com,
               "PISM ERROR: value of surface_gradient_method, option -gradient %s, not valid ... ending\n",
               method.c_str());
    PISMEnd();
  }

  return 0;
}

//! \brief Compute the ice surface gradient using the eta-transformation.
PetscErrorCode SIA_Sliding::surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;

  const PetscScalar n = flow_law->exponent(), // presumably 3.0
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0);
  const PetscScalar dx = grid.dx, dy = grid.dy;  // convenience
  IceModelVec2S eta = work_2d;

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
PetscErrorCode SIA_Sliding::surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
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
PetscErrorCode SIA_Sliding::surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
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
