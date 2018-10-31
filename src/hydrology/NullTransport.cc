// Copyright (C) 2012-2018 PISM Authors
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

#include "NullTransport.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/pism_utilities.hh" // clip

namespace pism {
namespace hydrology {

NullTransport::NullTransport(IceGrid::ConstPtr g)
  : Hydrology(g) {
  m_diffuse_tillwat    = m_config->get_boolean("hydrology.null_diffuse_till_water");
  m_diffusion_time     = m_config->get_double("hydrology.null_diffusion_time", "seconds");
  m_diffusion_distance = m_config->get_double("hydrology.null_diffusion_distance", "meters");
  m_tillwat_max        = m_config->get_double("hydrology.tillwat_max", "meters");
  m_tillwat_decay_rate = m_config->get_double("hydrology.tillwat_decay_rate", "m / second");

  if (m_tillwat_max < 0.0) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "hydrology::NullTransport: hydrology.tillwat_max is negative.\n"
                       "This is not allowed.");
  }

  if (m_diffuse_tillwat) {
    m_Wtill_old.create(m_grid, "Wtill_old", WITH_GHOSTS);
  }
}

NullTransport::~NullTransport() {
  // empty
}

void NullTransport::initialization_message() const {
  m_log->message(2,
                 "* Initializing the null-transport (till only) subglacial hydrology model ...\n");

  if (m_diffuse_tillwat) {
    m_log->message(2,
                   "  [using lateral diffusion of stored till water as in Bueler and Brown, 2009]\n");
  }
}

void NullTransport::restart_impl(const PIO &input_file, int record) {
  Hydrology::restart_impl(input_file, record);
}

void NullTransport::bootstrap_impl(const PIO &input_file,
                                   const IceModelVec2S &ice_thickness) {
  Hydrology::bootstrap_impl(input_file, ice_thickness);
}

void NullTransport::initialize_impl(const IceModelVec2S &W_till,
                                    const IceModelVec2S &W,
                                    const IceModelVec2S &P) {
  Hydrology::initialize_impl(W_till, W, P);
}

MaxTimestep NullTransport::max_timestep_impl(double t) const {
  (void) t;
  if (m_diffuse_tillwat) {
    const double
      dx2 = m_grid->dx() * m_grid->dx(),
      dy2 = m_grid->dy() * m_grid->dy(),
      L   = m_diffusion_distance,
      T   = m_diffusion_time,
      K   = L * L / (2.0 * T);

    return MaxTimestep(dx2 * dy2 / (2.0 * K * (dx2 + dy2)), "null-transport hydrology");
  } else {
    return MaxTimestep("null-transport hydrology");
  }
}

//! Update the till water thickness by simply integrating the melt input.
/*!
  Does a step of the trivial integration

  \f[ \frac{\partial W_{till}}{\partial t} = \frac{m}{\rho_w} - C\f]

  where \f$C=\f$`hydrology_tillwat_decay_rate`.  Enforces bounds
  \f$0 \le W_{till} \le W_{till}^{max}\f$ where the upper bound is
  `hydrology_tillwat_max`.  Here \f$m/\rho_w\f$ is `total_input`.

  Uses the current mass-continuity timestep `dt`.  (Compare
  hydrology::Routing::update_Wtill() which will generally be taking time steps
  determined by the evolving transportable water layer in that model.)

  There is no attempt to report on conservation errors because this
  hydrology::NullTransport model does not conserve water.

  There is no tranportable water thickness variable and no interaction with it.
*/
void NullTransport::update_impl(double t, double dt, const Inputs& inputs) {
  (void) t;

  // no transportable basal water
  m_W.set(0.0);

  m_input_change.add(dt, m_input_rate);

  const double
    water_density = m_config->get_double("constants.fresh_water.density"),
    kg_per_m      = m_grid->cell_area() * water_density; // kg m-1

  const IceModelVec2CellType &cell_type = *inputs.cell_type;

  IceModelVec::AccessList list{&cell_type, &m_Wtill, &m_input_rate,
      &m_conservation_error_change};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      W_old      = m_Wtill(i, j),
      dW_input   = dt * m_input_rate(i, j);

    if (W_old < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "negative subglacial water thickness of %f m at (%d, %d)",
                                    W_old, i, j);
    }

    // decay rate in areas under grounded ice
    double dW_decay = dt * (- m_tillwat_decay_rate);

    if (not cell_type.grounded_ice(i, j)) {
      dW_decay = 0.0;
    }

    m_Wtill(i, j) = (W_old + dW_input) + dW_decay;

    // dW_decay is a "conservation error"
    m_conservation_error_change(i, j) += dW_decay * kg_per_m;
  }

  if (m_diffuse_tillwat) {
    diffuse_till_water(dt);
  }

  // remove water in ice-free areas and account for changes
  enforce_bounds(*inputs.cell_type,
                 inputs.no_model_mask,
                 m_tillwat_max,
                 m_Wtill,
                 m_grounded_margin_change,
                 m_grounding_line_change,
                 m_conservation_error_change,
                 m_no_model_mask_change);
}

void NullTransport::diffuse_till_water(double dt) {
  // note: this call updates ghosts of m_Wtill_old
  m_Wtill_old.copy_from(m_Wtill);

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy(),
    L  = m_diffusion_distance,
    T  = m_diffusion_time,
    K  = L * L / (2.0 * T),
    Rx = K * dt / (dx * dx),
    Ry = K * dt / (dy * dy);

  IceModelVec::AccessList list{&m_Wtill, &m_Wtill_old, &m_flow_change};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto W = m_Wtill_old.star(i, j);

    double dWtill = (Rx * (W.w - 2.0 * W.ij + W.e) + Ry * (W.s - 2.0 * W.ij + W.n));

    m_Wtill(i, j) = W.ij + dWtill;

    m_flow_change(i, j) = dWtill;
  }
}

} // end of namespace hydrology
} // end of namespace pism
