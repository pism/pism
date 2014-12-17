// Copyright (C) 2004--2014 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "error_handling.hh"

namespace pism {

SIA_Sliding::SIA_Sliding(IceGrid &g, EnthalpyConverter &e)
  : ShallowStressBalance(g, e)
{
  verification_mode = false;
  eisII_experiment = "";

  const unsigned int WIDE_STENCIL = m_config.get("grid_max_stencil_width");

  for (int i = 0; i < 2; ++i) {
    char namestr[30];

    work_2d_stag[i].create(m_grid, "work_vector", WITH_GHOSTS);
    snprintf(namestr, sizeof(namestr), "work_vector_2d_stag_%d", i);
    work_2d_stag[i].set_name(namestr);

  }

  work_2d.create(m_grid, "work_vector_2d", WITH_GHOSTS, WIDE_STENCIL);

  {
    IceFlowLawFactory ice_factory(m_grid.com, "sia_", m_config, &EC);

    ice_factory.setType(m_config.get_string("sia_flow_law"));

    ice_factory.setFromOptions();
    flow_law = ice_factory.create();
  }
}

SIA_Sliding::~SIA_Sliding()
{
  if (flow_law != NULL) {
    delete flow_law;
    flow_law = NULL;
  }
}

void SIA_Sliding::init(Vars &vars) {

  ShallowStressBalance::init(vars);

  standard_gravity = m_config.get("standard_gravity");
  verification_mode = m_config.get_flag("sia_sliding_verification_mode");

  if (m_config.is_set("EISMINT_II_experiment")) {
    eisII_experiment = m_config.get_string("EISMINT_II_experiment");
  }

  thickness = vars.get_2d_scalar("land_ice_thickness");
  mask      = vars.get_2d_mask("mask");
  surface   = vars.get_2d_scalar("surface_altitude");
  bed       = vars.get_2d_scalar("bedrock_altitude");
  enthalpy  = vars.get_3d_scalar("enthalpy");
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
void SIA_Sliding::update(bool fast, IceModelVec2S &melange_back_pressure) {
  IceModelVec2Stag &h_x = work_2d_stag[0], &h_y = work_2d_stag[1];

  (void) fast;
  (void) melange_back_pressure;

  compute_surface_gradient(h_x, h_y);

  double mu_sliding = m_config.get("mu_sliding"),
    minimum_temperature_for_sliding = m_config.get("minimum_temperature_for_sliding"),
    ice_rho = m_config.get("ice_density");

  MaskQuery m(*mask);

  IceModelVec::AccessList list;
  list.add(h_x);
  list.add(h_y);

  list.add(*mask);
  list.add(*surface);
  list.add(*bed);
  list.add(*enthalpy);

  list.add(m_velocity);
  list.add(basal_frictional_heating);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m.ocean(i,j)) {
      m_velocity(i,j).u = 0.0;
      m_velocity(i,j).v = 0.0;
      basal_frictional_heating(i,j) = 0.0;
    } else {
      // basal velocity from SIA-type sliding law: not recommended!
      const double
        myx = m_grid.x(i),
        myy = m_grid.y(j),
        myhx = 0.25 * (h_x(i,j,0) + h_x(i-1,j,0)
                       + h_x(i,j,1) + h_x(i,j-1,1)),
        myhy = 0.25 * (h_y(i,j,0) + h_y(i-1,j,0)
                       + h_y(i,j,1) + h_y(i,j-1,1)),
        alpha = sqrt(PetscSqr(myhx) + PetscSqr(myhy));

      // change r1200: new meaning of H
      const double H = (*surface)(i,j) - (*bed)(i,j);

      double T = EC.getAbsTemp(enthalpy->getValZ(i,j,0.0), EC.getPressureFromDepth(H));

      double basalC = basalVelocitySIA(myx, myy, H, T,
                                       alpha, mu_sliding,
                                       minimum_temperature_for_sliding);
      m_velocity(i,j).u = - basalC * myhx;
      m_velocity(i,j).v = - basalC * myhy;
      // basal frictional heating; note P * dh/dx is x comp. of basal shear stress
      // in ice streams this result will be *overwritten* by
      //   correctBasalFrictionalHeating() if useSSAVelocities==TRUE
      const double P = ice_rho * standard_gravity * H;
      basal_frictional_heating(i,j) = - (P * myhx) * m_velocity(i,j).u - (P * myhy) * m_velocity(i,j).v;
    }
  }

  m_velocity.update_ghosts();
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

  Parameter \f$\mu\f$ can be set by option `-mu_sliding`.

  The returned coefficient is used in update() (above).
*/
double SIA_Sliding::basalVelocitySIA(double xIN, double yIN,
                                          double H, double T,
                                          double /*alpha*/, double mu,
                                          double min_T) const {
  double ice_rho = m_config.get("ice_density"),
    beta_CC_grad = m_config.get("beta_CC") * ice_rho * m_config.get("standard_gravity"),
    secpera = m_grid.convert(1.0, "year", "seconds");

  if (verification_mode) {
    // test 'E' mode
    const double r1 = 200e3, r2 = 700e3,   /* define region of sliding */
      theta1 = 10 * (M_PI/180), theta2 = 40 * (M_PI/180);
    const double x = fabs(xIN), y = fabs(yIN);
    const double r = sqrt(x * x + y * y);
    double       theta;
    if (x < 1.0) {
      theta = M_PI / 2.0;
    } else {
      theta = atan(y / x);
    }

    if ((r > r1) && (r < r2) && (theta > theta1) && (theta < theta2)) {
      // now INSIDE sliding region
      const double rbot = (r2 - r1) * (r2 - r1),
        thetabot = (theta2 - theta1) * (theta2 - theta1);
      const double mu_max = 2.5e-11; /* Pa^-1 m s^-1; max sliding coeff */
      double muE = mu_max * (4.0 * (r - r1) * (r2 - r) / rbot)
        * (4.0 * (theta - theta1) * (theta2 - theta) / thetabot);
      return muE * ice_rho * standard_gravity * H;
    } else {
      return 0.0;
    }
  }

  if ((eisII_experiment == "G") || (eisII_experiment == "H")) {
    const double  Bfactor = 1e-3 / secpera; // m s^-1 Pa^-1
    double pressure = EC.getPressureFromDepth(H);
    double E = EC.getEnthPermissive(T, 0.0, pressure);

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
    const double p_over = ice_rho * standard_gravity * H;
    return mu * p_over;
  } else {
    return 0;
  }
}

void SIA_Sliding::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {

  const std::string method = m_config.get_string("surface_gradient_method");

  if (method == "eta") {

    surface_gradient_eta(h_x, h_y);

  } else if (method == "haseloff") {

    surface_gradient_haseloff(h_x, h_y);

  } else if (method == "mahaffy") {

    surface_gradient_mahaffy(h_x, h_y);

  } else {
    throw RuntimeError::formatted("value of surface_gradient_method, option -gradient %s, is not valid",
                                  method.c_str());
  }
}

//! \brief Compute the ice surface gradient using the eta-transformation.
void SIA_Sliding::surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  const double n = flow_law->exponent(), // presumably 3.0
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0);
  const double dx = m_grid.dx(), dy = m_grid.dy();  // convenience
  IceModelVec2S &eta = work_2d;

  // compute eta = H^{8/3}, which is more regular, on reg grid

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(eta);

  unsigned int GHOSTS = eta.get_stencil_width();
  assert(thickness->get_stencil_width() >= GHOSTS);

  for (PointsWithGhosts p(m_grid, GHOSTS); p; p.next()) {
    const int i = p.i(), j = p.j();

    eta(i,j) = pow((*thickness)(i,j), etapow);
  }

  list.add(h_x);
  list.add(h_y);
  list.add(*bed);

  // now use Mahaffy on eta to get grad h on the staggered grid;
  // note   grad h = (3/8) eta^{-5/8} grad eta + grad b  because  h = H + b

  assert(h_x.get_stencil_width()  >= 1);
  assert(h_y.get_stencil_width()  >= 1);
  assert(eta.get_stencil_width()  >= 2);
  assert(bed->get_stencil_width() >= 2);

  for (int o=0; o<2; o++) {

    for (PointsWithGhosts p(m_grid, GHOSTS); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (o==0) {     // If I-offset
        const double mean_eta = 0.5 * (eta(i+1,j) + eta(i,j));
        if (mean_eta > 0.0) {
          const double factor = invpow * pow(mean_eta, dinvpow);
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
        const double mean_eta = 0.5 * (eta(i,j+1) + eta(i,j));
        if (mean_eta > 0.0) {
          const double factor = invpow * pow(mean_eta, dinvpow);
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

//! \brief Compute the ice surface gradient using Marianne Haseloff's approach.
/*!
 * Deals correctly with adjacent ice-free points with bed elevations which are
 * above the surface of the ice
 */
void SIA_Sliding::surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  const double Hicefree = 0.0;  // standard for ice-free, in Haseloff
  const double dx = m_grid.dx(), dy = m_grid.dy();  // convenience

  IceModelVec2S
    &b = *bed,
    &H = *thickness,
    &h = *surface;

  IceModelVec::AccessList list;
  list.add(h_x);
  list.add(h_y);
  list.add(*bed);
  list.add(*thickness);
  list.add(*surface);

  assert(h_x.get_stencil_width() >= 1);
  assert(h_y.get_stencil_width() >= 1);
  assert(bed->get_stencil_width()       >= 2);
  assert(thickness->get_stencil_width() >= 2);
  assert(surface->get_stencil_width()   >= 2);

  for (int o=0; o<2; o++) {

    for (PointsWithGhosts p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (o==0) {     // If I-offset
        const bool icefreeP  = (H(i,j)     <= Hicefree),
          icefreeE  = (H(i+1,j)   <= Hicefree),
          icefreeN  = (H(i,j+1)   <= Hicefree),
          icefreeS  = (H(i,j-1)   <= Hicefree),
          icefreeNE = (H(i+1,j+1) <= Hicefree),
          icefreeSE = (H(i+1,j-1) <= Hicefree);

        double hhE = h(i+1,j);  // east pseudo-surface elevation
        if (icefreeE  && (b(i+1,j)   > h(i,j))) {
          hhE  = h(i,j);
        }
        if (icefreeP  && (b(i,j)     > h(i+1,j))) {
          hhE  = h(i,j);
        }
        h_x(i,j,o) = (hhE - h(i,j)) / dx;

        double hhN  = h(i,j+1);  // north pseudo-surface elevation
        if (icefreeN  && (b(i,j+1)   > h(i,j))) {
          hhN  = h(i,j);
        }
        if (icefreeP  && (b(i,j)     > h(i,j+1))) {
          hhN  = h(i,j);
        }
        double hhS  = h(i,j-1);  // south pseudo-surface elevation
        if (icefreeS  && (b(i,j-1)   > h(i,j))) {
          hhS  = h(i,j);
        }
        if (icefreeP  && (b(i,j)     > h(i,j-1))) {
          hhS  = h(i,j);
        }
        double hhNE = h(i+1,j+1);// northeast pseudo-surface elevation
        if (icefreeNE && (b(i+1,j+1) > h(i+1,j))) {
          hhNE = h(i+1,j);
        }
        if (icefreeE  && (b(i+1,j)   > h(i+1,j+1))) {
          hhNE = h(i+1,j);
        }
        double hhSE = h(i+1,j-1);// southeast pseudo-surface elevation
        if (icefreeSE && (b(i+1,j-1) > h(i+1,j))) {
          hhSE = h(i+1,j);
        }
        if (icefreeE  && (b(i+1,j)   > h(i+1,j-1))) {
          hhSE = h(i+1,j);
        }
        h_y(i,j,o) = (hhNE + hhN - hhSE - hhS) / (4.0 * dy);
      } else {        // J-offset
        const bool icefreeP  = (H(i,j)     <= Hicefree),
          icefreeN  = (H(i,j+1)   <= Hicefree),
          icefreeE  = (H(i+1,j)   <= Hicefree),
          icefreeW  = (H(i-1,j)   <= Hicefree),
          icefreeNE = (H(i+1,j+1) <= Hicefree),
          icefreeNW = (H(i-1,j+1) <= Hicefree);

        double hhN  = h(i,j+1);  // north pseudo-surface elevation
        if (icefreeN  && (b(i,j+1)   > h(i,j))) {
          hhN  = h(i,j);
        }
        if (icefreeP  && (b(i,j)     > h(i,j+1))) {
          hhN  = h(i,j);
        }
        h_y(i,j,o) = (hhN - h(i,j)) / dy;

        double hhE  = h(i+1,j);  // east pseudo-surface elevation
        if (icefreeE  && (b(i+1,j)   > h(i,j))) {
          hhE  = h(i,j);
        }
        if (icefreeP  && (b(i,j)     > h(i+1,j))) {
          hhE  = h(i,j);
        }
        double hhW  = h(i-1,j);  // west pseudo-surface elevation
        if (icefreeW  && (b(i-1,j)   > h(i,j))) {
          hhW  = h(i,j);
        }
        if (icefreeP  && (b(i,j)     > h(i-1,j))) {
          hhW  = h(i,j);
        }
        double hhNE = h(i+1,j+1);// northeast pseudo-surface elevation
        if (icefreeNE && (b(i+1,j+1) > h(i,j+1))) {
          hhNE = h(i,j+1);
        }
        if (icefreeN  && (b(i,j+1)   > h(i+1,j+1))) {
          hhNE = h(i,j+1);
        }
        double hhNW = h(i-1,j+1);// northwest pseudo-surface elevation
        if (icefreeNW && (b(i-1,j+1) > h(i,j+1))) {
          hhNW = h(i,j+1);
        }
        if (icefreeN  && (b(i,j+1)   > h(i-1,j+1))) {
          hhNW = h(i,j+1);
        }
        h_x(i,j,o) = (hhNE + hhE - hhNW - hhW) / (4.0 * dx);
      }

    }
  }     // o
}

//! \brief Compute the ice surface gradient using the Mary Anne Mahaffy method;
//! see [\ref Mahaffy].
void SIA_Sliding::surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  const double dx = m_grid.dx(), dy = m_grid.dy();  // convenience

  IceModelVec2S &h = *surface;

  IceModelVec::AccessList list;
  list.add(h_x);
  list.add(h_y);
  list.add(*surface);

  assert(h_x.get_stencil_width() >= 1);
  assert(h_y.get_stencil_width() >= 1);
  assert(surface->get_stencil_width() >= 2);

  for (int o=0; o<2; o++) {
    for (PointsWithGhosts p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (o==0) {     // If I-offset
        h_x(i,j,o) = (h(i+1,j) - h(i,j)) / dx;
        h_y(i,j,o) = (+ h(i+1,j+1) + h(i,j+1)
                      - h(i+1,j-1) - h(i,j-1)) / (4.0*dy);
      } else {        // J-offset
        h_y(i,j,o) = (h(i,j+1) - h(i,j)) / dy;
        h_x(i,j,o) = (+ h(i+1,j+1) + h(i+1,j)
                      - h(i-1,j+1) - h(i-1,j)) / (4.0*dx);
      }
    }
  }
}

} // end of namespace pism
