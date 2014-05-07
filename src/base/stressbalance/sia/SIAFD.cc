// Copyright (C) 2004--2014 Jed Brown, Craig Lingle, Ed Bueler and Constantine Khroulev
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

#include "SIAFD.hh"
#include "Mask.hh"
#include "PISMBedSmoother.hh"
#include "enthalpyConverter.hh"
#include "PISMVars.hh"
#include "flowlaw_factory.hh"

namespace pism {

SIAFD::~SIAFD() {
  delete bed_smoother;
  if (flow_law != NULL) {
    delete flow_law;
    flow_law = NULL;
  }
}

//! \brief Allocate the SIAFD module.
PetscErrorCode SIAFD::allocate() {
  PetscErrorCode ierr;

  // 2D temporary storage:
  for (int i = 0; i < 2; ++i) {
    char namestr[30];

    ierr = work_2d[i].create(grid, "work_vector", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
    ierr = work_2d_stag[i].create(grid, "work_vector", WITH_GHOSTS); CHKERRQ(ierr);

    snprintf(namestr, sizeof(namestr), "work_vector_2d_%d", i);
    ierr = work_2d[i].set_name(namestr); CHKERRQ(ierr);

    for (int j = 0; j < 2; ++j) {
      snprintf(namestr, sizeof(namestr), "work_vector_2d_stag_%d_%d", i, j);
      ierr = work_2d_stag[i].set_name(namestr, j); CHKERRQ(ierr);
    }
  }

  ierr = delta[0].create(grid, "delta_0", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = delta[1].create(grid, "delta_1", WITH_GHOSTS); CHKERRQ(ierr);

  // 3D temporary storage:
  ierr = work_3d[0].create(grid, "work_3d_0", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = work_3d[1].create(grid, "work_3d_1", WITH_GHOSTS); CHKERRQ(ierr);

  // bed smoother
  bed_smoother = new BedSmoother(grid, config, WIDE_STENCIL);

  second_to_kiloyear = grid.convert(1, "second", "1000 years");

  {
    IceFlowLawFactory ice_factory(grid.com, "sia_", config, &EC);

    ierr = ice_factory.setType(config.get_string("sia_flow_law")); CHKERRQ(ierr);

    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ierr = ice_factory.create(&flow_law); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Initialize the SIA module.
PetscErrorCode SIAFD::init(Vars &vars) {
  PetscErrorCode ierr;

  ierr = SSB_Modifier::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the SIA stress balance modifier...\n"); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
                    "  [using the %s flow law]\n", flow_law->name().c_str()); CHKERRQ(ierr);

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (surface == NULL) SETERRQ(grid.com, 1, "surface_altitude is not available");

  bed = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (bed == NULL) SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(grid.com, 1, "enthalpy is not available");

  if (config.get_flag("do_age")) {
    age = dynamic_cast<IceModelVec3*>(vars.get("age"));
    if (age == NULL) SETERRQ(grid.com, 1, "age is not available");
  } else {
    age = NULL;
  }

  // set bed_state_counter to -1 so that the smoothed bed is computed the first
  // time update() is called.
  bed_state_counter = -1;
  return 0;
}

//! \brief Do the update; if fast == true, skip the update of 3D velocities and
//! strain heating.
PetscErrorCode SIAFD::update(IceModelVec2V *vel_input, bool fast) {
  PetscErrorCode ierr;
  IceModelVec2Stag &h_x = work_2d_stag[0], &h_y = work_2d_stag[1];

  // Check if the smoothed bed computed by BedSmoother is out of date and
  // recompute if necessary.
  if (bed->get_state_counter() > bed_state_counter) {
    ierr = bed_smoother->preprocess_bed(*bed); CHKERRQ(ierr);
    bed_state_counter = bed->get_state_counter();
  }

  ierr = compute_surface_gradient(h_x, h_y); CHKERRQ(ierr);

  ierr = compute_diffusive_flux(h_x, h_y, diffusive_flux, fast); CHKERRQ(ierr);

  if (!fast) {
    ierr = compute_3d_horizontal_velocity(h_x, h_y, vel_input, u, v); CHKERRQ(ierr);
  }

  return 0;
}


//! \brief Compute the ice surface gradient for the SIA.
/*!
  There are three methods for computing the surface gradient. Which method is
  controlled by configuration parameter `surface_gradient_method` which can
  have values `haseloff`, `mahaffy`, or `eta`.

  The most traditional method is to directly differentiate the surface
  elevation \f$h\f$ by the Mahaffy method [\ref Mahaffy]. The `haseloff` method,
  suggested by Marianne Haseloff, modifies the Mahaffy method only where
  ice-free adjacent bedrock points are above the ice surface, and in those
  cases the returned gradient component is zero.

  The alternative method, when `surface_gradient_method` = `eta`, transforms
  the thickness to something more regular and differentiates that. We get back
  to the gradient of the surface by applying the chain rule. In particular, as
  shown in [\ref CDDSV] for the flat bed and \f$n=3\f$ case, if we define

  \f[\eta = H^{(2n+2)/n}\f]

  then \f$\eta\f$ is more regular near the margin than \f$H\f$. So we compute
  the surface gradient by

  \f[\nabla h = \frac{n}{(2n+2)} \eta^{(-n-2)/(2n+2)} \nabla \eta + \nabla b,\f]

  recalling that \f$h = H + b\f$. This method is only applied when \f$\eta >
  0\f$ at a given point; otherwise \f$\nabla h = \nabla b\f$.

  In all cases we are computing the gradient by finite differences onto a
  staggered grid. In the method with \f$\eta\f$ we apply centered differences
  using (roughly) the same method for \f$\eta\f$ and \f$b\f$ that applies
  directly to the surface elevation \f$h\f$ in the `mahaffy` and `haseloff`
  methods.

  \param[out] h_x the X-component of the surface gradient, on the staggered grid
  \param[out] h_y the Y-component of the surface gradient, on the staggered grid
*/
PetscErrorCode SIAFD::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode  ierr;

  const std::string method = config.get_string("surface_gradient_method");

  if (method == "eta") {

    ierr = surface_gradient_eta(h_x, h_y); CHKERRQ(ierr);

  } else if (method == "haseloff") {

    ierr = surface_gradient_haseloff(h_x, h_y); CHKERRQ(ierr);

  } else if (method == "mahaffy") {

    ierr = surface_gradient_mahaffy(h_x, h_y); CHKERRQ(ierr);

  } else {
    verbPrintf(1, grid.com,
               "PISM ERROR: value of surface_gradient_method, option '-gradient %s', not valid ... ending\n",
               method.c_str());
    PISMEnd();
  }

  return 0;
}

//! \brief Compute the ice surface gradient using the eta-transformation.
PetscErrorCode SIAFD::surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;

  const double n = flow_law->exponent(), // presumably 3.0
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0);
  const double dx = grid.dx, dy = grid.dy;  // convenience
  IceModelVec2S &eta = work_2d[0];

  // compute eta = H^{8/3}, which is more regular, on reg grid
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = eta.begin_access(); CHKERRQ(ierr);
  int GHOSTS = 2;
  for (int   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (int j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
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
  for (int o=0; o<2; o++) {

    GHOSTS = 1;
    for (int   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (int j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
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
  ierr = eta.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);

  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  return 0;
}


//! \brief Compute the ice surface gradient using the Mary Anne Mahaffy method;
//! see [\ref Mahaffy].
PetscErrorCode SIAFD::surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;
  const double dx = grid.dx, dy = grid.dy;  // convenience

  IceModelVec2S &h = *surface;
  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);
  ierr = surface->begin_access(); CHKERRQ(ierr);

  for (int o=0; o<2; o++) {
    int GHOSTS = 1;
    for (int   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (int j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
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

  ierr = surface->end_access(); CHKERRQ(ierr);
  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the ice surface gradient using a modification of Marianne Haseloff's approach.
/*!
 * The original code deals correctly with adjacent ice-free points with bed
 * elevations which are above the surface of the ice nearby. This is done by
 * setting surface gradient at the margin to zero at such locations.
 *
 * This code also deals with shelf fronts: sharp surface elevation change at
 * the ice shelf front would otherwise cause abnormally high diffusivity
 * values, which forces PISM to take shorter time-steps than necessary. (Note
 * that the mass continuity code does not use SIA fluxes in floating areas.)
 * This is done by assuming that the ice surface near shelf fronts is
 * horizontal (i.e. here the surface gradient is set to zero also).
 *
 * The code below uses an interpretation of the standard Mahaffy scheme. We
 * compute components of the surface gradient at staggered grid locations. The
 * field h_x stores the x-component on the i-offset and j-offset grids, h_y ---
 * the y-component.
 *
 * The Mahaffy scheme for the x-component at grid points on the i-offset grid
 * (offset in the x-direction) is just the centered finite difference using
 * adjacent regular-grid points. (Similarly for the y-component at j-offset
 * locations.)
 *
 * Mahaffy's prescription for computing the y-component on the i-offset can be
 * interpreted as:
 *
 * - compute the y-component at four surrounding j-offset staggered grid locations,
 * - compute the average of these four.
 *
 * The code below does just that.
 *
 * - The first double for-loop computes x-components at i-offset
 *   locations and y-components at j-offset locations. Each computed
 *   number is assigned a weight (w_i and w_j) that is used below
 *
 * - The second double for-loop computes x-components at j-offset
 *   locations and y-components at i-offset locations as averages of
 *   quantities computed earlier. The weight are used to keep track of
 *   the number of values used in the averaging process.
 *
 * This method communicates ghost values of h_x and h_y. They cannot be
 * computed locally because the first loop uses width=2 stencil of surface,
 * mask, and bed to compute values at all grid points including width=1 ghosts,
 * then the second loop uses width=1 stencil to compute local values. (In other
 * words, a purely local computation would require width=3 stencil of surface,
 * mask, and bed fields.)
 */
PetscErrorCode SIAFD::surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;
  const double dx = grid.dx, dy = grid.dy;  // convenience
  IceModelVec2S &h = *surface, &b = *bed,
    &w_i = work_2d[0], &w_j = work_2d[1]; // averaging weights

  MaskQuery m(*mask);

  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);
  ierr = w_i.begin_access(); CHKERRQ(ierr);
  ierr = w_j.begin_access(); CHKERRQ(ierr);

  ierr = h.begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = b.begin_access(); CHKERRQ(ierr);
  int GHOSTS = 1;
  for (int   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (int j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {

      // x-derivative, i-offset
      {
        if ((m.floating_ice(i,j) && m.ice_free_ocean(i+1,j)) ||
            (m.ice_free_ocean(i,j) && m.floating_ice(i+1,j))) {
          // marine margin
          h_x(i,j,0) = 0;
          w_i(i,j)   = 0;
        } else if ((m.icy(i,j) && m.ice_free(i+1,j) && b(i+1,j) > h(i,j)) ||
                   (m.ice_free(i,j) && m.icy(i+1,j) && b(i,j) > h(i+1,j))) {
          // ice next to a "cliff"
          h_x(i,j,0) = 0.0;
          w_i(i,j)   = 0;
        } else {
          // default case
          h_x(i,j,0) = (h(i+1,j) - h(i,j)) / dx;
          w_i(i,j)   = 1;
        }
      }

      // y-derivative, j-offset
      {
        if ((m.floating_ice(i,j) && m.ice_free_ocean(i,j+1)) ||
            (m.ice_free_ocean(i,j) && m.floating_ice(i,j+1))) {
          // marine margin
          h_y(i,j,1) = 0.0;
          w_j(i,j)   = 0.0;
        } else if ((m.icy(i,j) && m.ice_free(i,j+1) && b(i,j+1) > h(i,j)) ||
                   (m.ice_free(i,j) && m.icy(i,j+1) && b(i,j) > h(i,j+1))) {
          // ice next to a "cliff"
          h_y(i,j,1) = 0.0;
          w_j(i,j)   = 0.0;
        } else {
          // default case
          h_y(i,j,1) = (h(i,j+1) - h(i,j)) / dy;
          w_j(i,j)   = 1.0;
        }
      }
    } // inner loop (j)
  } // outer loop (i)
  ierr = b.end_access(); CHKERRQ(ierr);
  ierr = h.end_access(); CHKERRQ(ierr);

  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      // x-derivative, j-offset
      {
        if (w_j(i,j) > 0) {
          double W = w_i(i,j) + w_i(i-1,j) + w_i(i-1,j+1) + w_i(i,j+1);
          if (W > 0)
            h_x(i,j,1) = 1.0/W * (h_x(i,j,0) + h_x(i-1,j,0) + h_x(i-1,j+1,0) + h_x(i,j+1,0));
          else
            h_x(i,j,1) = 0.0;
        } else {
          if (m.icy(i,j)) {
            double W = w_i(i,j) + w_i(i-1,j);
            if (W > 0)
              h_x(i,j,1) = 1.0/W * (h_x(i,j,0) + h_x(i-1,j,0));
            else
              h_x(i,j,1) = 0.0;
          } else {
            double W = w_i(i,j+1) + w_i(i-1,j+1);
            if (W > 0)
              h_x(i,j,1) = 1.0/W * (h_x(i-1,j+1,0) + h_x(i,j+1,0));
            else
              h_x(i,j,1) = 0.0;
          }
        }
      } // end of "x-derivative, j-offset"

      // y-derivative, i-offset
      {
        if (w_i(i,j) > 0) {
          double W = w_j(i,j) + w_j(i,j-1) + w_j(i+1,j-1) + w_j(i+1,j);
          if (W > 0)
            h_y(i,j,0) = 1.0/W * (h_y(i,j,1) + h_y(i,j-1,1) + h_y(i+1,j-1,1) + h_y(i+1,j,1));
          else
            h_y(i,j,0) = 0.0;
        } else {
          if (m.icy(i,j)) {
            double W = w_j(i,j) + w_j(i,j-1);
            if (W > 0)
              h_y(i,j,0) = 1.0/W * (h_y(i,j,1) + h_y(i,j-1,1));
            else
              h_y(i,j,0) = 0.0;
          } else {
            double W = w_j(i+1,j-1) + w_j(i+1,j);
            if (W > 0)
              h_y(i,j,0) = 1.0/W * (h_y(i+1,j-1,1) + h_y(i+1,j,1));
            else
              h_y(i,j,0) = 0.0;
          }
        }
      } // end of "y-derivative, i-offset"
    }   // inner loop (j)
  }     // outer loop (i)


  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = w_j.end_access(); CHKERRQ(ierr);
  ierr = w_i.end_access(); CHKERRQ(ierr);
  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  ierr = h_x.update_ghosts(); CHKERRQ(ierr);
  ierr = h_y.update_ghosts(); CHKERRQ(ierr);

  return 0;
}


//! \brief Compute the SIA flux. If fast == false, also store delta on the staggered grid.
/*!
 * Recall that \f$ Q = -D \nabla h \f$ is the diffusive flux in the mass-continuity equation
 *
 * \f[ \frac{\partial H}{\partial t} = M - S - \nabla \cdot (Q + \mathbf{U}_b H),\f]
 *
 * where \f$h\f$ is the ice surface elevation, \f$M\f$ is the top surface
 * accumulation/ablation rate, \f$S\f$ is the basal mass balance and
 * \f$\mathbf{U}_b\f$ is the thickness-advective (in PISM: usually SSA) ice
 * velocity.
 *
 * Recall also that at any particular point in the map-plane (i.e. if \f$x\f$
 * and \f$y\f$ are fixed)
 *
 * \f[ D = 2\int_b^h F(z)P(z)(h-z)dz, \f]
 *
 * where \f$F(z)\f$ is a constitutive function and \f$P(z)\f$ is the pressure
 * at a level \f$z\f$.
 *
 * By defining
 *
 * \f[ \delta(z) = 2F(z)P(z) \f]
 *
 * one can write
 *
 * \f[D = \int_b^h\delta(z)(h-z)dz. \f]
 *
 * The advantage is that it is then possible to avoid re-evaluating
 * \f$F(z)\f$ (which is computationally expensive) in the horizontal ice
 * velocity (see compute_3d_horizontal_velocity()) computation.
 *
 * This method computes \f$Q\f$ and stores \f$\delta\f$ in delta[0,1] is fast == false.
 *
 * The trapezoidal rule is used to approximate the integral.
 *
 * \param[in]  h_x x-component of the surface gradient, on the staggered grid
 * \param[in]  h_y y-component of the surface gradient, on the staggered grid
 * \param[out] result diffusive ice flux
 * \param[in]  fast the boolean flag specitying if we're doing a "fast" update.
 */
PetscErrorCode SIAFD::compute_diffusive_flux(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y,
                                             IceModelVec2Stag &result, bool fast) {
  PetscErrorCode  ierr;
  IceModelVec2S &thk_smooth = work_2d[0],
    &theta = work_2d[1];

  bool full_update = (fast == false);

  ierr = result.set(0.0); CHKERRQ(ierr);

  double *delta_ij;
  delta_ij = new double[grid.Mz];

  const double enhancement_factor = flow_law->enhancement_factor();
  double ice_grain_size = config.get("ice_grain_size");

  bool compute_grain_size_using_age = config.get_flag("compute_grain_size_using_age");

  // some flow laws use grain size, and even need age to update grain size
  if (compute_grain_size_using_age && (!config.get_flag("do_age"))) {
    PetscPrintf(grid.com,
                "PISM ERROR in SIAFD::compute_diffusive_flux(): do_age not set but\n"
                "age is needed for grain-size-based flow law ...  ENDING! ...\n\n");
    PISMEnd();
  }

  const bool use_age = (IceFlowLawUsesGrainSize(flow_law) &&
                        compute_grain_size_using_age &&
                        config.get_flag("do_age"));

  // get "theta" from Schoof (2003) bed smoothness calculation and the
  // thickness relative to the smoothed bed; each IceModelVec2S involved must
  // have stencil width WIDE_GHOSTS for this too work
  ierr = bed_smoother->get_theta(*surface, &theta); CHKERRQ(ierr);

  ierr = bed_smoother->get_smoothed_thk(*surface, *thickness, *mask,
                                        &thk_smooth); CHKERRQ(ierr);

  ierr = theta.begin_access(); CHKERRQ(ierr);
  ierr = thk_smooth.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);

  double *age_ij, *age_offset;
  if (use_age) {
    ierr = age->begin_access(); CHKERRQ(ierr);
  }

  if (full_update) {
    ierr = delta[0].begin_access(); CHKERRQ(ierr);
    ierr = delta[1].begin_access(); CHKERRQ(ierr);
  }

  double *E_ij, *E_offset;
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);

  double my_D_max = 0.0;
  for (int o=0; o<2; o++) {
    int GHOSTS = 1;
    for (int   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (int j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        // staggered point: o=0 is i+1/2, o=1 is j+1/2, (i,j) and (i+oi,j+oj)
        //   are regular grid neighbors of a staggered point:
        const int oi = 1 - o, oj = o;

        const double
          thk = 0.5 * (thk_smooth(i,j) + thk_smooth(i+oi,j+oj));

        // zero thickness case:
        if (thk == 0.0) {
          result(i,j,o) = 0.0;
          if (full_update) {
            ierr = delta[o].setColumn(i, j, 0.0); CHKERRQ(ierr);
          }
          continue;
        }

        if (use_age) {
          ierr = age->getInternalColumn(i, j, &age_ij); CHKERRQ(ierr);
          ierr = age->getInternalColumn(i+oi, j+oj, &age_offset); CHKERRQ(ierr);
        }

        ierr = enthalpy->getInternalColumn(i, j, &E_ij); CHKERRQ(ierr);
        ierr = enthalpy->getInternalColumn(i+oi, j+oj, &E_offset); CHKERRQ(ierr);

        const double slope = (o==0) ? h_x(i,j,o) : h_y(i,j,o);
        const int      ks = grid.kBelowHeight(thk);
        const double   alpha =
          sqrt(PetscSqr(h_x(i,j,o)) + PetscSqr(h_y(i,j,o)));
        const double theta_local = 0.5 * (theta(i,j) + theta(i+oi,j+oj));

        double  Dfoffset = 0.0;  // diffusivity for deformational SIA flow
        for (int k = 0; k <= ks; ++k) {
          double depth = thk - grid.zlevels[k]; // FIXME issue #15
          // pressure added by the ice (i.e. pressure difference between the
          // current level and the top of the column)
          const double pressure = EC.getPressureFromDepth(depth);

          double flow;
          if (use_age) {
            ice_grain_size = grainSizeVostok(0.5 * (age_ij[k] + age_offset[k]));
          }
          // If the flow law does not use grain size, it will just ignore it,
          // no harm there
          double E = 0.5 * (E_ij[k] + E_offset[k]);
          flow = flow_law->flow(alpha * pressure, E, pressure, ice_grain_size);

          delta_ij[k] = enhancement_factor * theta_local * 2.0 * pressure * flow;

          if (k > 0) { // trapezoidal rule
            const double dz = grid.zlevels[k] - grid.zlevels[k-1];
            Dfoffset += 0.5 * dz * ((depth + dz) * delta_ij[k-1] + depth * delta_ij[k]);
          }
        }
        // finish off D with (1/2) dz (0 + (H-z[ks])*delta_ij[ks]), but dz=H-z[ks]:
        const double dz = thk - grid.zlevels[ks];
        Dfoffset += 0.5 * dz * dz * delta_ij[ks];

        // Override diffusivity at the edges of the domain. (At these
        // locations PISM uses ghost cells *beyond* the boundary of
        // the computational domain. This does not matter if the ice
        // does not extend all the way to the domain boundary, as in
        // whole-ice-sheet simulations. In a regional setup, though,
        // this adjustment lets us avoid taking very small time-steps
        // because of the possible thickness and bed elevation
        // "discontinuities" at the boundary.)
        if (i < 0 || i >= grid.Mx-1 || j < 0 || j >= grid.My-1)
          Dfoffset = 0.0;

        my_D_max = PetscMax(my_D_max, Dfoffset);

        // vertically-averaged SIA-only flux, sans sliding; note
        //   result(i,j,0) is  u  at E (east)  staggered point (i+1/2,j)
        //   result(i,j,1) is  v  at N (north) staggered point (i,j+1/2)
        result(i,j,o) = - Dfoffset * slope;

        // if doing the full update, fill the delta column above the ice and
        // store it:
        if (full_update) {
          for (unsigned int k = ks + 1; k < grid.Mz; ++k) {
            delta_ij[k] = 0.0;
          }
          ierr = delta[o].setInternalColumn(i,j,delta_ij); CHKERRQ(ierr);
        }
      } // o
    } // j
  } // i

  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = theta.end_access(); CHKERRQ(ierr);
  ierr = thk_smooth.end_access(); CHKERRQ(ierr);

  if (use_age) {
    ierr = age->end_access(); CHKERRQ(ierr);
  }

  ierr = enthalpy->end_access(); CHKERRQ(ierr);

  if (full_update) {
    ierr = delta[1].end_access(); CHKERRQ(ierr);
    ierr = delta[0].end_access(); CHKERRQ(ierr);
  }

  ierr = GlobalMax(&my_D_max, &D_max, grid.com); CHKERRQ(ierr);

  delete [] delta_ij;

  return 0;
}

//! \brief Compute diffusivity (diagnostically).
/*!
 * Computes \f$D\f$ as
 *
 * \f[D = \int_b^h\delta(z)(h-z)dz. \f]
 *
 * Uses the trapezoidal rule to approximate the integral.
 *
 * See compute_diffusive_flux() for the rationale and the definition of
 * \f$\delta\f$.
 * \param[out] result The diffusivity of the SIA flow.
 */
PetscErrorCode SIAFD::compute_diffusivity(IceModelVec2S &result) {
  PetscErrorCode ierr;
  IceModelVec2Stag &D_stag = work_2d_stag[0];

  ierr = this->compute_diffusivity_staggered(D_stag); CHKERRQ(ierr);

  ierr = D_stag.update_ghosts(); CHKERRQ(ierr);

  ierr = D_stag.staggered_to_regular(result); CHKERRQ(ierr);

  return 0;
}

/*!
 * \brief Computes the diffusivity of the SIA mass continuity equation on the
 * staggered grid (for debugging).
 */
PetscErrorCode SIAFD::compute_diffusivity_staggered(IceModelVec2Stag &D_stag) {
  PetscErrorCode ierr;

  // delta on the staggered grid:
  double *delta_ij;
  IceModelVec2S &thk_smooth = work_2d[0];

  ierr = bed_smoother->get_smoothed_thk(*surface, *thickness, *mask,
                                        &thk_smooth); CHKERRQ(ierr);

  ierr = thk_smooth.begin_access(); CHKERRQ(ierr);
  ierr = delta[0].begin_access(); CHKERRQ(ierr);
  ierr = delta[1].begin_access(); CHKERRQ(ierr);
  ierr = D_stag.begin_access(); CHKERRQ(ierr);
  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      for (int o = 0; o < 2; ++o) {
        const int oi = 1 - o, oj = o;

        ierr = delta[o].getInternalColumn(i,j,&delta_ij); CHKERRQ(ierr);

        const double
          thk = 0.5 * (thk_smooth(i,j) + thk_smooth(i+oi,j+oj));

        if (thk == 0) {
          D_stag(i,j,o) = 0.0;
          continue;
        }

        const unsigned int ks = grid.kBelowHeight(thk);
        double Dfoffset = 0.0;

        for (unsigned int k = 1; k <= ks; ++k) {
          double depth = thk - grid.zlevels[k];

          const double dz = grid.zlevels[k] - grid.zlevels[k-1];
          // trapezoidal rule
          Dfoffset += 0.5 * dz * ((depth + dz) * delta_ij[k-1] + depth * delta_ij[k]);
        }

        // finish off D with (1/2) dz (0 + (H-z[ks])*delta_ij[ks]), but dz=H-z[ks]:
        const double dz = thk - grid.zlevels[ks];
        Dfoffset += 0.5 * dz * dz * delta_ij[ks];

        D_stag(i,j,o) = Dfoffset;
      }
    }
  }
  ierr = D_stag.end_access(); CHKERRQ(ierr);
  ierr = delta[1].end_access(); CHKERRQ(ierr);
  ierr = delta[0].end_access(); CHKERRQ(ierr);
  ierr = thk_smooth.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Extend the grid vertically.
PetscErrorCode SIAFD::extend_the_grid(int old_Mz) {
  PetscErrorCode ierr;

  ierr = SSB_Modifier::extend_the_grid(old_Mz); CHKERRQ(ierr);

  ierr = delta[0].extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = delta[1].extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);

  ierr = work_3d[0].extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = work_3d[1].extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute I.
/*!
 * This computes
 * \f[ I(z) = \int_b^z\delta(s)ds.\f]
 *
 * Uses the trapezoidal rule to approximate the integral.
 *
 * See compute_diffusive_flux() for the definition of \f$\delta\f$.
 *
 * The result is stored in work_3d[0,1] and is used to compute the SIA component
 * of the 3D-distributed horizontal ice velocity.
 */
PetscErrorCode SIAFD::compute_I() {
  PetscErrorCode ierr;
  double *I_ij, *delta_ij;

  IceModelVec2S &thk_smooth = work_2d[0];
  IceModelVec3* I = work_3d;

  ierr = bed_smoother->get_smoothed_thk(*surface, *thickness, *mask,
                                        &thk_smooth); CHKERRQ(ierr);

  ierr = delta[0].begin_access(); CHKERRQ(ierr);
  ierr = delta[1].begin_access(); CHKERRQ(ierr);
  ierr = I[0].begin_access(); CHKERRQ(ierr);
  ierr = I[1].begin_access(); CHKERRQ(ierr);
  ierr = thk_smooth.begin_access(); CHKERRQ(ierr);

  for (int o = 0; o < 2; ++o) {
    int GHOSTS = 1;
    for (int   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
      for (int j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
        const int oi = 1-o, oj=o;
        const double
          thk = 0.5 * (thk_smooth(i,j) + thk_smooth(i+oi,j+oj));

        ierr = delta[o].getInternalColumn(i,j,&delta_ij); CHKERRQ(ierr);
        ierr = I[o].getInternalColumn(i,j,&I_ij); CHKERRQ(ierr);

        const unsigned int ks = grid.kBelowHeight(thk);

        // within the ice:
        I_ij[0] = 0.0;
        double I_current = 0.0;
        for (unsigned int k = 1; k <= ks; ++k) {
          const double dz = grid.zlevels[k] - grid.zlevels[k-1];
          // trapezoidal rule
          I_current += 0.5 * dz * (delta_ij[k-1] + delta_ij[k]);
          I_ij[k] = I_current;
        }
        // above the ice:
        for (unsigned int k = ks + 1; k < grid.Mz; ++k) {
          I_ij[k] = I_current;
        }
      }
    }
  }

  ierr = thk_smooth.end_access(); CHKERRQ(ierr);
  ierr = I[1].end_access(); CHKERRQ(ierr);
  ierr = I[0].end_access(); CHKERRQ(ierr);
  ierr = delta[1].end_access(); CHKERRQ(ierr);
  ierr = delta[0].end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute horizontal components of the SIA velocity (in 3D).
/*!
 * Recall that
 *
 * \f[ \mathbf{U}(z) = -2 \nabla h \int_b^z F(s)P(s)ds + \mathbf{U}_b,\f]
 *
 * which can be written in terms of \f$I(z)\f$ defined in compute_I():
 *
 * \f[ \mathbf{U}(z) = -I(z) \nabla h + \mathbf{U}_b. \f]
 *
 * \note This is one of the places where "hybridization" is done.
 *
 * \param[in] h_x the X-component of the surface gradient, on the staggered grid
 * \param[in] h_y the Y-component of the surface gradient, on the staggered grid
 * \param[in] vel_input the thickness-advective velocity from the underlying stress balance module
 * \param[out] u_out the X-component of the resulting horizontal velocity field
 * \param[out] v_out the Y-component of the resulting horizontal velocity field
 */
PetscErrorCode SIAFD::compute_3d_horizontal_velocity(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y,
                                                     IceModelVec2V *vel_input,
                                                     IceModelVec3 &u_out, IceModelVec3 &v_out) {
  PetscErrorCode ierr;

  ierr = compute_I(); CHKERRQ(ierr);
  // after the compute_I() call work_3d[0,1] contains I on the staggered grid
  IceModelVec3 *I = work_3d;

  double *u_ij, *v_ij, *IEAST, *IWEST, *INORTH, *ISOUTH;

  ierr = u_out.begin_access(); CHKERRQ(ierr);
  ierr = v_out.begin_access(); CHKERRQ(ierr);

  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);
  ierr = vel_input->begin_access(); CHKERRQ(ierr);

  ierr = I[0].begin_access(); CHKERRQ(ierr);
  ierr = I[1].begin_access(); CHKERRQ(ierr);

  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = I[0].getInternalColumn(i, j, &IEAST); CHKERRQ(ierr);
      ierr = I[0].getInternalColumn(i - 1, j, &IWEST); CHKERRQ(ierr);
      ierr = I[1].getInternalColumn(i, j, &INORTH); CHKERRQ(ierr);
      ierr = I[1].getInternalColumn(i, j - 1, &ISOUTH); CHKERRQ(ierr);

      ierr = u_out.getInternalColumn(i, j, &u_ij); CHKERRQ(ierr);
      ierr = v_out.getInternalColumn(i, j, &v_ij); CHKERRQ(ierr);

      // Fetch values from 2D fields *outside* of the k-loop:
      double h_x_w = h_x(i - 1, j, 0), h_x_e = h_x(i, j, 0),
        h_x_n = h_x(i, j, 1), h_x_s = h_x(i, j - 1, 1);

      double h_y_w = h_y(i - 1, j, 0), h_y_e = h_y(i, j, 0),
        h_y_n = h_y(i, j, 1), h_y_s = h_y(i, j - 1, 1);

      double vel_input_u = (*vel_input)(i, j).u,
        vel_input_v = (*vel_input)(i, j).v;

      for (unsigned int k = 0; k < grid.Mz; ++k) {
        u_ij[k] = - 0.25 * (IEAST[k]  * h_x_e + IWEST[k]  * h_x_w +
                             INORTH[k] * h_x_n + ISOUTH[k] * h_x_s);
        v_ij[k] = - 0.25 * (IEAST[k]  * h_y_e + IWEST[k]  * h_y_w +
                             INORTH[k] * h_y_n + ISOUTH[k] * h_y_s);

        // Add the "SSA" velocity:
        u_ij[k] += vel_input_u;
        v_ij[k] += vel_input_v;
      }
    }
  }

  ierr = I[1].end_access(); CHKERRQ(ierr);
  ierr = I[0].end_access(); CHKERRQ(ierr);

  ierr = vel_input->end_access(); CHKERRQ(ierr);
  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  ierr = u_out.end_access(); CHKERRQ(ierr);
  ierr = v_out.end_access(); CHKERRQ(ierr);

  // Communicate to get ghosts:
  ierr = u_out.update_ghosts(); CHKERRQ(ierr);
  ierr = v_out.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

//! Use the Vostok core as a source of a relationship between the age of the ice and the grain size.
/*! A data set is interpolated here. The intention is that the softness of the
  ice has nontrivial dependence on its age, through its grainsize, because of
  variable dustiness of the global climate. The grainsize is partly determined
  by at which point in the glacial cycle the given ice fell as snow.

  The data is from [\ref DeLaChapelleEtAl98] and [\ref LipenkovEtAl89]. In
  particular, Figure A2 in the former reference was hand-sampled with an
  attempt to include the ``wiggles'' in that figure. Ages of the oldest ice (>=
  300 ka) were estimated in a necessarily ad hoc way. The age value of 10000 ka
  was added simply to give interpolation for very old ice; ages beyond that get
  constant extrapolation. Linear interpolation is done between the samples.
 */
double SIAFD::grainSizeVostok(double age_seconds) const {
  const int numPoints = 22;
  const double ageAt[numPoints] = {  // ages in ka
    0.0000e+00, 5.0000e+01, 1.0000e+02, 1.2500e+02, 1.5000e+02,
    1.5800e+02, 1.6500e+02, 1.7000e+02, 1.8000e+02, 1.8800e+02,
    2.0000e+02, 2.2500e+02, 2.4500e+02, 2.6000e+02, 3.0000e+02,
    3.2000e+02, 3.5000e+02, 4.0000e+02, 5.0000e+02, 6.0000e+02,
    8.0000e+02, 1.0000e+04 };
  const double gsAt[numPoints] = {   // grain sizes in m
    1.8000e-03, 2.2000e-03, 3.0000e-03, 4.0000e-03, 4.3000e-03,
    3.0000e-03, 3.0000e-03, 4.6000e-03, 3.4000e-03, 3.3000e-03,
    5.9000e-03, 6.2000e-03, 5.4000e-03, 6.8000e-03, 3.5000e-03,
    6.0000e-03, 8.0000e-03, 8.3000e-03, 3.6000e-03, 3.8000e-03,
    9.5000e-03, 1.0000e-02 };
  const double a = age_seconds * second_to_kiloyear; // Age in ka
  int l = 0;               // Left end of the binary search
  int r = numPoints - 1;   // Right end

  // If we are out of range
  if (a < ageAt[l]) {
    return gsAt[l];
  } else if (a > ageAt[r]) {
    return gsAt[r];
  }
  // Binary search for the interval
  while (r > l + 1) {
    const int j = (r + l) / 2;
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

} // end of namespace pism
