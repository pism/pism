// Copyright (C) 2004--2018 Jed Brown, Craig Lingle, Ed Bueler and Constantine Khroulev
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

#include <cstdlib>
#include <cassert>

#include "SIAFD.hh"
#include "BedSmoother.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/stressbalance/StressBalance.hh"

#include "pism/util/Time.hh"

namespace pism {
namespace stressbalance {

SIAFD::SIAFD(IceGrid::ConstPtr g)
  : SSB_Modifier(g) {

  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

  // 2D temporary storage:
  for (int i = 0; i < 2; ++i) {
    char name[30];
    snprintf(name, sizeof(name), "work_vector_2d_%d", i);

    m_work_2d[i].create(m_grid, name, WITH_GHOSTS, WIDE_STENCIL);
  }

  m_h_x.create(m_grid, "h_x", WITH_GHOSTS);
  m_h_y.create(m_grid, "h_y", WITH_GHOSTS);
  m_D.create(m_grid, "diffusivity", WITH_GHOSTS);

  m_delta[0].create(m_grid, "delta_0", WITH_GHOSTS);
  m_delta[1].create(m_grid, "delta_1", WITH_GHOSTS);

  // 3D temporary storage:
  m_work_3d[0].create(m_grid, "work_3d_0", WITH_GHOSTS);
  m_work_3d[1].create(m_grid, "work_3d_1", WITH_GHOSTS);

  // bed smoother
  m_bed_smoother = new BedSmoother(m_grid, WIDE_STENCIL);

  m_second_to_kiloyear = units::convert(m_sys, 1, "second", "1000 years");

  {
    rheology::FlowLawFactory ice_factory("stress_balance.sia.", m_config, m_EC);
    m_flow_law = ice_factory.create();
  }

  const bool compute_grain_size_using_age = m_config->get_boolean("stress_balance.sia.grain_size_age_coupling");
  const bool age_model_enabled = m_config->get_boolean("age.enabled");
  const bool e_age_coupling = m_config->get_boolean("stress_balance.sia.e_age_coupling");

  if (compute_grain_size_using_age) {
    if (not FlowLawUsesGrainSize(*m_flow_law)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "flow law %s does not use grain size "
                                    "but sia.grain_size_age_coupling was set",
                                    m_flow_law->name().c_str());
    }

    if (not age_model_enabled) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "SIAFD: age model is not active but\n"
                                    "age is needed for grain-size-based flow law %s",
                                    m_flow_law->name().c_str());
    }
  }

  if (e_age_coupling and not age_model_enabled) {
      throw RuntimeError(PISM_ERROR_LOCATION, "SIAFD: age model is not active but\n"
                         "age is needed for age-dependent flow enhancement");
  }

  m_eemian_start   = m_config->get_double("time.eemian_start", "seconds");
  m_eemian_end     = m_config->get_double("time.eemian_end", "seconds");
  m_holocene_start = m_config->get_double("time.holocene_start", "seconds");
}

SIAFD::~SIAFD() {
  delete m_bed_smoother;
}

//! \brief Initialize the SIA module.
void SIAFD::init() {

  SSB_Modifier::init();

  m_log->message(2,
             "* Initializing the SIA stress balance modifier...\n");
  m_log->message(2,
             "  [using the %s flow law]\n", m_flow_law->name().c_str());


  // implements an option e.g. described in @ref Greve97Greenland that is the
  // enhancement factor is coupled to the age of the ice
  if (m_config->get_boolean("stress_balance.sia.e_age_coupling")) {
    m_log->message(2,
                   "  using age-dependent enhancement factor:\n"
                   "  e=%f for ice accumulated during interglacial periods\n"
                   "  e=%f for ice accumulated during glacial periods\n",
                   m_flow_law->enhancement_factor_interglacial(),
                   m_flow_law->enhancement_factor());
  }
}

//! \brief Do the update; if full_update == false skip the update of 3D velocities and strain
//! heating.
void SIAFD::update(const IceModelVec2V &sliding_velocity,
                   const Inputs &inputs,
                   bool full_update) {

  const Profiling &profiling = m_grid->ctx()->profiling();

  // Check if the smoothed bed computed by BedSmoother is out of date and
  // recompute if necessary.
  if (inputs.new_bed_elevation) {
    profiling.begin("sia.bed_smoother");
    m_bed_smoother->preprocess_bed(inputs.geometry->bed_elevation);
    profiling.end("sia.bed_smoother");
  }

  profiling.begin("sia.gradient");
  compute_surface_gradient(inputs, m_h_x, m_h_y);
  profiling.end("sia.gradient");

  profiling.begin("sia.flux");
  compute_diffusivity(full_update,
                      *inputs.geometry,
                      inputs.enthalpy,
                      inputs.age,
                      m_h_x, m_h_y, m_D);
  compute_diffusive_flux(m_h_x, m_h_y, m_D, m_diffusive_flux);
  profiling.end("sia.flux");

  if (full_update) {
    profiling.begin("sia.3d_velocity");
    compute_3d_horizontal_velocity(*inputs.geometry, m_h_x, m_h_y, sliding_velocity,
                                   m_u, m_v);
    profiling.end("sia.3d_velocity");
  }
}


//! \brief Compute the ice surface gradient for the SIA.
/*!
  There are three methods for computing the surface gradient. Which method is
  controlled by configuration parameter `sia.surface_gradient_method` which can
  have values `haseloff`, `mahaffy`, or `eta`.

  The most traditional method is to directly differentiate the surface
  elevation \f$h\f$ by the Mahaffy method [\ref Mahaffy]. The `haseloff` method,
  suggested by Marianne Haseloff, modifies the Mahaffy method only where
  ice-free adjacent bedrock points are above the ice surface, and in those
  cases the returned gradient component is zero.

  The alternative method, when `sia.surface_gradient_method` = `eta`, transforms
  the thickness to something more regular and differentiates that. We get back
  to the gradient of the surface by applying the chain rule. In particular, as
  shown in [\ref Vazquezetal2003] for the flat bed and \f$n=3\f$ case, if we define

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
void SIAFD::compute_surface_gradient(const Inputs &inputs,
                                     IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) const {

  const std::string method = m_config->get_string("stress_balance.sia.surface_gradient_method");

  if (method == "eta") {

    surface_gradient_eta(inputs, h_x, h_y);

  } else if (method == "haseloff") {

    surface_gradient_haseloff(inputs, h_x, h_y);

  } else if (method == "mahaffy") {

    surface_gradient_mahaffy(inputs, h_x, h_y);

  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "value of sia.surface_gradient_method, option '-gradient %s', is not valid",
                                  method.c_str());
  }
}

//! \brief Compute the ice surface gradient using the eta-transformation.
void SIAFD::surface_gradient_eta(const Inputs &inputs,
                                 IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) const {
  const double n = m_flow_law->exponent(), // presumably 3.0
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0);
  const double dx = m_grid->dx(), dy = m_grid->dy();  // convenience
  IceModelVec2S &eta = m_work_2d[0];

  // compute eta = H^{8/3}, which is more regular, on reg grid

  const IceModelVec2S
    &H = inputs.geometry->ice_thickness,
    &b = inputs.geometry->bed_elevation;

  IceModelVec::AccessList list{&eta, &H, &h_x, &h_y, &b};

  unsigned int GHOSTS = eta.stencil_width();
  assert(H.stencil_width() >= GHOSTS);

  for (PointsWithGhosts p(*m_grid, GHOSTS); p; p.next()) {
    const int i = p.i(), j = p.j();

    eta(i,j) = pow(H(i,j), etapow);
  }

  // now use Mahaffy on eta to get grad h on staggered;
  // note   grad h = (3/8) eta^{-5/8} grad eta + grad b  because  h = H + b

  assert(b.stencil_width()   >= 2);
  assert(eta.stencil_width() >= 2);
  assert(h_x.stencil_width() >= 1);
  assert(h_y.stencil_width() >= 1);

  for (int o=0; o<2; o++) {

    for (PointsWithGhosts p(*m_grid); p; p.next()) {
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
        h_x(i,j,o) += b.diff_x_stagE(i,j);
        h_y(i,j,o) += b.diff_y_stagE(i,j);
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
        h_y(i,j,o) += b.diff_y_stagN(i,j);
        h_x(i,j,o) += b.diff_x_stagN(i,j);
      }
    }
  }
}


//! \brief Compute the ice surface gradient using the Mary Anne Mahaffy method;
//! see [\ref Mahaffy].
void SIAFD::surface_gradient_mahaffy(const Inputs &inputs,
                                     IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) const {
  const double dx = m_grid->dx(), dy = m_grid->dy();  // convenience

  const IceModelVec2S &h = inputs.geometry->ice_surface_elevation;

  IceModelVec::AccessList list{&h_x, &h_y, &h};

  // h_x and h_y have to have ghosts
  assert(h_x.stencil_width() >= 1);
  assert(h_y.stencil_width() >= 1);
  // surface elevation needs more ghosts
  assert(h.stencil_width()   >= 2);

  for (PointsWithGhosts p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // I-offset
    h_x(i, j, 0) = (h(i + 1, j) - h(i, j)) / dx;
    h_y(i, j, 0) = (+ h(i + 1, j + 1) + h(i, j + 1)
                    - h(i + 1, j - 1) - h(i, j - 1)) / (4.0*dy);
    // J-offset
    h_y(i, j, 1) = (h(i, j + 1) - h(i, j)) / dy;
    h_x(i, j, 1) = (+ h(i + 1, j + 1) + h(i + 1, j)
                    - h(i - 1, j + 1) - h(i - 1, j)) / (4.0*dx);
  }
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
void SIAFD::surface_gradient_haseloff(const Inputs &inputs,
                                      IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) const {
  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();  // convenience
  const IceModelVec2S
    &h = inputs.geometry->ice_surface_elevation,
    &b = inputs.geometry->bed_elevation;
  IceModelVec2S
    &w_i = m_work_2d[0],
    &w_j = m_work_2d[1]; // averaging weights

  const IceModelVec2CellType &mask = inputs.geometry->cell_type;

  IceModelVec::AccessList list{&h_x, &h_y, &w_i, &w_j, &h, &mask, &b};

  assert(b.stencil_width()    >= 2);
  assert(mask.stencil_width() >= 2);
  assert(h.stencil_width()    >= 2);
  assert(h_x.stencil_width()  >= 1);
  assert(h_y.stencil_width()  >= 1);
  assert(w_i.stencil_width()  >= 1);
  assert(w_j.stencil_width()  >= 1);

  for (PointsWithGhosts p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // x-derivative, i-offset
    {
      if ((mask.floating_ice(i,j) && mask.ice_free_ocean(i+1,j)) ||
          (mask.ice_free_ocean(i,j) && mask.floating_ice(i+1,j))) {
        // marine margin
        h_x(i,j,0) = 0;
        w_i(i,j)   = 0;
      } else if ((mask.icy(i,j) && mask.ice_free(i+1,j) && b(i+1,j) > h(i,j)) ||
                 (mask.ice_free(i,j) && mask.icy(i+1,j) && b(i,j) > h(i+1,j))) {
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
      if ((mask.floating_ice(i,j) && mask.ice_free_ocean(i,j+1)) ||
          (mask.ice_free_ocean(i,j) && mask.floating_ice(i,j+1))) {
        // marine margin
        h_y(i,j,1) = 0.0;
        w_j(i,j)   = 0.0;
      } else if ((mask.icy(i,j) && mask.ice_free(i,j+1) && b(i,j+1) > h(i,j)) ||
                 (mask.ice_free(i,j) && mask.icy(i,j+1) && b(i,j) > h(i,j+1))) {
        // ice next to a "cliff"
        h_y(i,j,1) = 0.0;
        w_j(i,j)   = 0.0;
      } else {
        // default case
        h_y(i,j,1) = (h(i,j+1) - h(i,j)) / dy;
        w_j(i,j)   = 1.0;
      }
    }
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // x-derivative, j-offset
    {
      if (w_j(i,j) > 0) {
        double W = w_i(i,j) + w_i(i-1,j) + w_i(i-1,j+1) + w_i(i,j+1);
        if (W > 0) {
          h_x(i,j,1) = 1.0/W * (h_x(i,j,0) + h_x(i-1,j,0) + h_x(i-1,j+1,0) + h_x(i,j+1,0));
        } else {
          h_x(i,j,1) = 0.0;
        }
      } else {
        if (mask.icy(i,j)) {
          double W = w_i(i,j) + w_i(i-1,j);
          if (W > 0) {
            h_x(i,j,1) = 1.0/W * (h_x(i,j,0) + h_x(i-1,j,0));
          } else {
            h_x(i,j,1) = 0.0;
          }
        } else {
          double W = w_i(i,j+1) + w_i(i-1,j+1);
          if (W > 0) {
            h_x(i,j,1) = 1.0/W * (h_x(i-1,j+1,0) + h_x(i,j+1,0));
          } else {
            h_x(i,j,1) = 0.0;
          }
        }
      }
    } // end of "x-derivative, j-offset"

      // y-derivative, i-offset
    {
      if (w_i(i,j) > 0) {
        double W = w_j(i,j) + w_j(i,j-1) + w_j(i+1,j-1) + w_j(i+1,j);
        if (W > 0) {
          h_y(i,j,0) = 1.0/W * (h_y(i,j,1) + h_y(i,j-1,1) + h_y(i+1,j-1,1) + h_y(i+1,j,1));
        } else {
          h_y(i,j,0) = 0.0;
        }
      } else {
        if (mask.icy(i,j)) {
          double W = w_j(i,j) + w_j(i,j-1);
          if (W > 0) {
            h_y(i,j,0) = 1.0/W * (h_y(i,j,1) + h_y(i,j-1,1));
          } else {
            h_y(i,j,0) = 0.0;
          }
        } else {
          double W = w_j(i+1,j-1) + w_j(i+1,j);
          if (W > 0) {
            h_y(i,j,0) = 1.0/W * (h_y(i+1,j-1,1) + h_y(i+1,j,1));
          } else {
            h_y(i,j,0) = 0.0;
          }
        }
      }
    } // end of "y-derivative, i-offset"
  }

  h_x.update_ghosts();
  h_y.update_ghosts();
}


//! \brief Compute the SIA diffusivity. If full_update, also store delta on the staggered grid.
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
 * This method computes \f$D\f$ and stores \f$\delta\f$ in delta[0,1] if full_update is true.
 *
 * The trapezoidal rule is used to approximate the integral.
 *
 * \param[in]  full_update the boolean flag specitying if we're doing a "full" update.
 * \param[in]  h_x x-component of the surface gradient, on the staggered grid
 * \param[in]  h_y y-component of the surface gradient, on the staggered grid
 * \param[out] result diffusivity of the SIA flow
 */
void SIAFD::compute_diffusivity(bool full_update,
                                const Geometry &geometry,
                                const IceModelVec3 *enthalpy,
                                const IceModelVec3 *age,
                                const IceModelVec2Stag &h_x,
                                const IceModelVec2Stag &h_y,
                                IceModelVec2Stag &result) {
  IceModelVec2S
    &thk_smooth = m_work_2d[0],
    &theta      = m_work_2d[1];

  const IceModelVec2S
    &h = geometry.ice_surface_elevation,
    &H = geometry.ice_thickness;

  const IceModelVec2CellType &mask = geometry.cell_type;

  result.set(0.0);

  const double
    current_time                    = m_grid->ctx()->time()->current(),
    enhancement_factor              = m_flow_law->enhancement_factor(),
    enhancement_factor_interglacial = m_flow_law->enhancement_factor_interglacial(),
    D_limit                         = m_config->get_double("stress_balance.sia.max_diffusivity");

  const bool
    compute_grain_size_using_age = m_config->get_boolean("stress_balance.sia.grain_size_age_coupling"),
    e_age_coupling               = m_config->get_boolean("stress_balance.sia.e_age_coupling"),
    limit_diffusivity            = m_config->get_boolean("stress_balance.sia.limit_diffusivity"),
    use_age                      = compute_grain_size_using_age or e_age_coupling;

  // get "theta" from Schoof (2003) bed smoothness calculation and the
  // thickness relative to the smoothed bed; each IceModelVec2S involved must
  // have stencil width WIDE_GHOSTS for this too work
  m_bed_smoother->theta(h, theta);

  m_bed_smoother->smoothed_thk(h, H, mask, thk_smooth);

  IceModelVec::AccessList list{&result, &theta, &thk_smooth, &h_x, &h_y, enthalpy};

  if (use_age) {
    assert(age->stencil_width() >= 2);
    list.add(*age);
  }

  if (full_update) {
    list.add({&m_delta[0], &m_delta[1]});
    assert(m_delta[0].stencil_width()  >= 1);
    assert(m_delta[1].stencil_width()  >= 1);
  }

  assert(theta.stencil_width()      >= 2);
  assert(thk_smooth.stencil_width() >= 2);
  assert(result.stencil_width()     >= 1);
  assert(h_x.stencil_width()        >= 1);
  assert(h_y.stencil_width()        >= 1);
  assert(enthalpy->stencil_width()  >= 2);

  const std::vector<double> &z = m_grid->z();
  const unsigned int
    Mx = m_grid->Mx(),
    My = m_grid->My(),
    Mz = m_grid->Mz();

  std::vector<double> depth(Mz), stress(Mz), pressure(Mz), E(Mz), flow(Mz);
  std::vector<double> delta_ij(Mz);
  std::vector<double> A(Mz), ice_grain_size(Mz, m_config->get_double("constants.ice.grain_size", "m"));
  std::vector<double> e_factor(Mz, enhancement_factor);

  double D_max = 0.0;
  int high_diffusivity_counter = 0;
  for (int o=0; o<2; o++) {
    ParallelSection loop(m_grid->com);
    try {
      for (PointsWithGhosts p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        // staggered point: o=0 is i+1/2, o=1 is j+1/2, (i, j) and (i+oi, j+oj)
        //   are regular grid neighbors of a staggered point:
        const int oi = 1 - o, oj = o;

        const double
          thk = 0.5 * (thk_smooth(i, j) + thk_smooth(i+oi, j+oj));

        // zero thickness case:
        if (thk == 0.0) {
          result(i, j, o) = 0.0;
          if (full_update) {
            m_delta[o].set_column(i, j, 0.0);
          }
          continue;
        }

        const int ks = m_grid->kBelowHeight(thk);

        for (int k = 0; k <= ks; ++k) {
          depth[k] = thk - z[k];
        }

        // pressure added by the ice (i.e. pressure difference between the
        // current level and the top of the column)
        m_EC->pressure(depth, ks, pressure); // FIXME issue #15

        if (use_age) {
          const double
            *age_ij     = age->get_column(i, j),
            *age_offset = age->get_column(i+oi, j+oj);

          for (int k = 0; k <= ks; ++k) {
            A[k] = 0.5 * (age_ij[k] + age_offset[k]);
          }

          if (compute_grain_size_using_age) {
            for (int k = 0; k <= ks; ++k) {
              ice_grain_size[k] = grainSizeVostok(A[k]);
            }
          }

          if (e_age_coupling) {
            for (int k = 0; k <= ks; ++k) {
              const double accumulation_time = current_time - A[k];
              if (interglacial(accumulation_time)) {
                e_factor[k] = enhancement_factor_interglacial;
              } else {
                e_factor[k] = enhancement_factor;
              }
            }
          }
        }

        {
          const double
            *E_ij     = enthalpy->get_column(i, j),
            *E_offset = enthalpy->get_column(i+oi, j+oj);
          for (int k = 0; k <= ks; ++k) {
            E[k] = 0.5 * (E_ij[k] + E_offset[k]);
          }
        }

        const double alpha = sqrt(PetscSqr(h_x(i, j, o)) + PetscSqr(h_y(i, j, o)));
        for (int k = 0; k <= ks; ++k) {
          stress[k] = alpha * pressure[k];
        }

        m_flow_law->flow_n(&stress[0], &E[0], &pressure[0], &ice_grain_size[0], ks + 1,
                           &flow[0]);

        const double theta_local = 0.5 * (theta(i, j) + theta(i+oi, j+oj));
        for (int k = 0; k <= ks; ++k) {
          delta_ij[k] = e_factor[k] * theta_local * 2.0 * pressure[k] * flow[k];
        }

        double D = 0.0;  // diffusivity for deformational SIA flow
        {
          for (int k = 1; k <= ks; ++k) {
            // trapezoidal rule
            const double dz = z[k] - z[k-1];
            D += 0.5 * dz * ((depth[k] + dz) * delta_ij[k-1] + depth[k] * delta_ij[k]);
          }
          // finish off D with (1/2) dz (0 + (H-z[ks])*delta_ij[ks]), but dz=H-z[ks]:
          const double dz = thk - z[ks];
          D += 0.5 * dz * dz * delta_ij[ks];
        }

        // Override diffusivity at the edges of the domain. (At these
        // locations PISM uses ghost cells *beyond* the boundary of
        // the computational domain. This does not matter if the ice
        // does not extend all the way to the domain boundary, as in
        // whole-ice-sheet simulations. In a regional setup, though,
        // this adjustment lets us avoid taking very small time-steps
        // because of the possible thickness and bed elevation
        // "discontinuities" at the boundary.)
        if (i < 0 || i >= (int)Mx - 1 ||
            j < 0 || j >= (int)My - 1) {
          D = 0.0;
        }

        if (limit_diffusivity and D >= D_limit) {
          D = D_limit;
          high_diffusivity_counter += 1;
        }

        D_max = std::max(D_max, D);

        result(i, j, o) = D;

        // if doing the full update, fill the delta column above the ice and
        // store it:
        if (full_update) {
          for (unsigned int k = ks + 1; k < Mz; ++k) {
            delta_ij[k] = 0.0;
          }
          m_delta[o].set_column(i, j, &delta_ij[0]);
        }
      } // i, j-loop
    } catch (...) {
      loop.failed();
    }
    loop.check();
  } // o-loop

  m_D_max = GlobalMax(m_grid->com, D_max);

  high_diffusivity_counter = GlobalSum(m_grid->com, high_diffusivity_counter);

  if (m_D_max > D_limit) {

    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Maximum diffusivity of SIA flow (%f m2/s) is too high.\n"
                                  "This probably means that the bed elevation or the ice thickness is "
                                  "too rough.\n"
                                  "Increase stress_balance.sia.max_diffusivity to suppress this message.", m_D_max);

  } else if (high_diffusivity_counter > 0) {

    m_log->message(2, "  SIA diffusivity was capped at %.2f m2/s at %d locations.\n",
                   D_limit, high_diffusivity_counter);
  }
}

void SIAFD::compute_diffusive_flux(const IceModelVec2Stag &h_x, const IceModelVec2Stag &h_y,
                                   const IceModelVec2Stag &diffusivity,
                                   IceModelVec2Stag &result) {

  IceModelVec::AccessList list{&diffusivity, &h_x, &h_y, &result};

  for (int o = 0; o < 2; o++) {
    ParallelSection loop(m_grid->com);
    try {
      for (PointsWithGhosts p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        const double slope = (o == 0) ? h_x(i, j, o) : h_y(i, j, o);

        result(i, j, o) = - diffusivity(i, j, o) * slope;
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();
  } // o-loop
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
void SIAFD::compute_I(const Geometry &geometry) {

  IceModelVec2S &thk_smooth = m_work_2d[0];
  IceModelVec3* I = m_work_3d;

  const IceModelVec2S
    &h = geometry.ice_surface_elevation,
    &H = geometry.ice_thickness;

  const IceModelVec2CellType &mask = geometry.cell_type;

  m_bed_smoother->smoothed_thk(h, H, mask, thk_smooth);

  IceModelVec::AccessList list{&m_delta[0], &m_delta[1], &I[0], &I[1], &thk_smooth};

  assert(I[0].stencil_width()     >= 1);
  assert(I[1].stencil_width()     >= 1);
  assert(m_delta[0].stencil_width() >= 1);
  assert(m_delta[1].stencil_width() >= 1);
  assert(thk_smooth.stencil_width() >= 2);

  const unsigned int Mz = m_grid->Mz();

  std::vector<double> dz(Mz);
  for (unsigned int k = 1; k < Mz; ++k) {
    dz[k] = m_grid->z(k) - m_grid->z(k - 1);
  }

  for (int o = 0; o < 2; ++o) {
    ParallelSection loop(m_grid->com);
    try {
      for (PointsWithGhosts p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        const int oi = 1 - o, oj = o;
        const double
          thk = 0.5 * (thk_smooth(i, j) + thk_smooth(i + oi, j + oj));

        const double *delta_ij = m_delta[o].get_column(i, j);
        double       *I_ij     = I[o].get_column(i, j);

        const unsigned int ks = m_grid->kBelowHeight(thk);

        // within the ice:
        I_ij[0] = 0.0;
        double I_current = 0.0;
        for (unsigned int k = 1; k <= ks; ++k) {
          // trapezoidal rule
          I_current += 0.5 * dz[k] * (delta_ij[k - 1] + delta_ij[k]);
          I_ij[k] = I_current;
        }

        // above the ice:
        for (unsigned int k = ks + 1; k < Mz; ++k) {
          I_ij[k] = I_current;
        }
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();
  } // o-loop
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
 * \param[in] sliding_velocity the thickness-advective velocity from the underlying stress balance module
 * \param[out] u_out the X-component of the resulting horizontal velocity field
 * \param[out] v_out the Y-component of the resulting horizontal velocity field
 */
void SIAFD::compute_3d_horizontal_velocity(const Geometry &geometry,
                                           const IceModelVec2Stag &h_x,
                                           const IceModelVec2Stag &h_y,
                                           const IceModelVec2V &sliding_velocity,
                                           IceModelVec3 &u_out, IceModelVec3 &v_out) {

  compute_I(geometry);
  // after the compute_I() call work_3d[0,1] contains I on the staggered grid
  IceModelVec3 *I = m_work_3d;

  IceModelVec::AccessList list{&u_out, &v_out, &h_x, &h_y, &sliding_velocity, &I[0], &I[1]};

  const unsigned int Mz = m_grid->Mz();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      *I_e = I[0].get_column(i, j),
      *I_w = I[0].get_column(i - 1, j),
      *I_n = I[1].get_column(i, j),
      *I_s = I[1].get_column(i, j - 1);

    // Fetch values from 2D fields *outside* of the k-loop:
    const double
      h_x_w = h_x(i - 1, j, 0),
      h_x_e = h_x(i, j, 0),
      h_x_n = h_x(i, j, 1),
      h_x_s = h_x(i, j - 1, 1);

    const double
      h_y_w = h_y(i - 1, j, 0),
      h_y_e = h_y(i, j, 0),
      h_y_n = h_y(i, j, 1),
      h_y_s = h_y(i, j - 1, 1);

    const double
      sliding_velocity_u = sliding_velocity(i, j).u,
      sliding_velocity_v = sliding_velocity(i, j).v;

    double
      *u_ij = u_out.get_column(i, j),
      *v_ij = v_out.get_column(i, j);

    // split into two loops to encourage auto-vectorization
    for (unsigned int k = 0; k < Mz; ++k) {
      u_ij[k] = sliding_velocity_u - 0.25 * (I_e[k] * h_x_e + I_w[k] * h_x_w +
                                             I_n[k] * h_x_n + I_s[k] * h_x_s);
    }
    for (unsigned int k = 0; k < Mz; ++k) {
      v_ij[k] = sliding_velocity_v - 0.25 * (I_e[k] * h_y_e + I_w[k] * h_y_w +
                                             I_n[k] * h_y_n + I_s[k] * h_y_s);
    }
  }

  // Communicate to get ghosts:
  u_out.update_ghosts();
  v_out.update_ghosts();
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

  FIXME: Use GSL's interpolation code.
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
  const double a = age_seconds * m_second_to_kiloyear; // Age in ka
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
  if ((r == l) || (std::abs(r - l) > 1)) {
    throw RuntimeError(PISM_ERROR_LOCATION, "binary search in grainSizeVostok: oops");
  }
  // Linear interpolation on the interval
  return gsAt[l] + (a - ageAt[l]) * (gsAt[r] - gsAt[l]) / (ageAt[r] - ageAt[l]);
}

//! Determine if `accumulation_time` corresponds to an interglacial period.
bool SIAFD::interglacial(double accumulation_time) {
  if (accumulation_time < m_eemian_start) {
    return false;
  } else if (accumulation_time < m_eemian_end) {
    return true;
  } else if (accumulation_time < m_holocene_start) {
    return false;
  } else {
    return true;
  }
}

const IceModelVec2Stag& SIAFD::surface_gradient_x() const {
  return m_h_x;
}

const IceModelVec2Stag& SIAFD::surface_gradient_y() const {
  return m_h_y;
}

const IceModelVec2Stag& SIAFD::diffusivity() const {
  return m_D;
}

const BedSmoother& SIAFD::bed_smoother() const {
  return *m_bed_smoother;
}


} // end of namespace stressbalance
} // end of namespace pism
