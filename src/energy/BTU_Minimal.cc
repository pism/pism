/* Copyright (C) 2016, 2017 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "BTU_Minimal.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace energy {

BTU_Minimal::BTU_Minimal(IceGrid::ConstPtr g)
  : BedThermalUnit(g) {
  // empty
}

void BTU_Minimal::init_impl(const InputOptions &opts) {
  m_log->message(2,
                 "* Initializing the minimal model for lithosphere:\n"
                 "  stored time-independent geothermal flux applied to ice base...\n");

  BedThermalUnit::init_impl(opts);

  // The flux through the top surface is the same as the flux through the bottom surface.
  m_top_surface_flux.copy_from(m_bottom_surface_flux);
}

double BTU_Minimal::vertical_spacing_impl() const {
  return 0.0;
}

double BTU_Minimal::depth_impl() const {
  return 0.0;
}

unsigned int BTU_Minimal::Mz_impl() const {
  return 0;
}

MaxTimestep BTU_Minimal::max_timestep_impl(double t) const {
  (void) t;
  // no time step restriction
  return MaxTimestep("minimal thermal bedrock layer");
}

void BTU_Minimal::update_impl(const array::Scalar &bedrock_top_temperature, double t, double dt) {
  (void) bedrock_top_temperature;
  (void) t;
  (void) dt;
  // empty
}

} // end of namespace energy
} // end of namespace pism
