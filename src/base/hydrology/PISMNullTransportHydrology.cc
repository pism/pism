// Copyright (C) 2012-2016 PISM Authors
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

#include "PISMHydrology.hh"
#include "base/util/Mask.hh"
#include "base/util/error_handling.hh"
#include "hydrology_diagnostics.hh"
#include "base/util/PISMVars.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {
namespace hydrology {

NullTransport::NullTransport(IceGrid::ConstPtr g)
  : Hydrology(g) {
  m_diffuse_tillwat    = m_config->get_boolean("hydrology_null_diffuse_till_water");
  m_diffusion_time     = m_config->get_double("hydrology_null_diffusion_time", "seconds");
  m_diffusion_distance = m_config->get_double("hydrology_null_diffusion_distance", "meters");
  m_tillwat_max        = m_config->get_double("hydrology_tillwat_max", "meters");
  m_tillwat_decay_rate = m_config->get_double("hydrology_tillwat_decay_rate");

  if (m_tillwat_max < 0.0) {
    throw RuntimeError("hydrology::NullTransport: hydrology_tillwat_max is negative.\n"
                       "This is not allowed.");
  }

  if (m_diffuse_tillwat) {
    m_Wtil_old.create(m_grid, "Wtil_old", WITH_GHOSTS);
  }
}

NullTransport::~NullTransport() {
}

void NullTransport::init() {
  m_log->message(2,
             "* Initializing the null-transport (till only) subglacial hydrology model ...\n");

  if (m_diffuse_tillwat) {
    m_log->message(2,
                   "  [using lateral diffusion of stored till water as in Bueler and Brown, 2009]\n");
  }

  Hydrology::init();
}

MaxTimestep NullTransport::max_timestep_impl(double t) {
  (void) t;
  if (m_diffuse_tillwat) {
    const double
      dx2 = m_grid->dx() * m_grid->dx(),
      dy2 = m_grid->dy() * m_grid->dy(),
      L   = m_diffusion_distance,
      T   = m_diffusion_time,
      K   = L * L / (2.0 * T);

    return MaxTimestep(dx2 * dy2 / (2.0 * K * (dx2 + dy2)));
  } else {
    return MaxTimestep();
  }
}

//! Set the transportable subglacial water thickness to zero; there is no tranport.
void NullTransport::subglacial_water_thickness(IceModelVec2S &result) {
  result.set(0.0);
}


//! Returns the (trivial) overburden pressure as the pressure of the non-existent transportable water, because this is the least harmful output if this is misused.
void NullTransport::subglacial_water_pressure(IceModelVec2S &result) {
  overburden_pressure(result);
}


//! Update the till water thickness by simply integrating the melt input.
/*!
Does a step of the trivial integration
  \f[ \frac{\partial W_{til}}{\partial t} = \frac{m}{\rho_w} - C\f]
where \f$C=\f$`hydrology_tillwat_decay_rate`.  Enforces bounds
\f$0 \le W_{til} \le W_{til}^{max}\f$ where the upper bound is
`hydrology_tillwat_max`.  Here \f$m/\rho_w\f$ is `total_input`.

Uses the current mass-continuity timestep `dt`.  (Compare
hydrology::Routing::raw_update_Wtil() which will generally be taking time steps
determined by the evolving transportable water layer in that model.)

There is no attempt to report on conservation errors because this
hydrology::NullTransport model does not conserve water.

There is no tranportable water thickness variable and no interaction with it.
 */
void NullTransport::update_impl(double t, double dt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(t - m_t) < 1e-6) and (fabs(dt - m_dt) < 1e-6)) {
    return;
  }
  m_t = t;
  m_dt = dt;

  get_input_rate(t, dt, m_total_input);

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(m_Wtil);
  list.add(m_total_input);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i,j) || mask.ice_free(i,j)) {
      m_Wtil(i, j) = 0.0;
    } else {
      m_Wtil(i, j) += dt * (m_total_input(i, j) - m_tillwat_decay_rate);
      m_Wtil(i, j) = std::min(std::max(0.0, m_Wtil(i, j)), m_tillwat_max);
    }
  }

  if (m_diffuse_tillwat) {
    diffuse_till_water(dt);
  }
}

void NullTransport::diffuse_till_water(double dt) {
  // note: this call updates ghosts of m_Wtil_old
  m_Wtil_old.copy_from(m_Wtil);

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy(),
    L  = m_diffusion_distance,
    T  = m_diffusion_time,
    K  = L * L / (2.0 * T),
    Rx = K * dt / (dx * dx),
    Ry = K * dt / (dy * dy);

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(m_Wtil);
  list.add(m_Wtil_old);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const StarStencil<double> W = m_Wtil_old.star(i, j);

    // diffuse Wtil
    m_Wtil(i, j) = ((1.0 - 2.0 * Rx - 2.0 * Ry) * W.ij + Rx * (W.w + W.e) + Ry * (W.s + W.n));

    // enforce bounds
    m_Wtil(i, j) = std::min(std::max(0.0, m_Wtil(i, j)), m_tillwat_max);

    // set to zero in
    if (mask.ocean(i, j) or mask.ice_free(i, j)) {
      m_Wtil(i, j) = 0.0;
    }
  }

}

} // end of namespace hydrology
} // end of namespace pism
