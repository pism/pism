/* Copyright (C) 2016, 2017, 2018 PISM Authors
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

#ifndef POINITIALIZATION_H
#define POINITIALIZATION_H

#include "pism/coupler/OceanModel.hh"

namespace pism {
namespace ocean {

/*! Ocean model "modifier" that helps with initialization.
 *
 * This modifier saves *all* fields a ocean model provides as a part of the model state and re-loads
 * them during initialization so that they are available *before* the first time step in a
 * re-started run.
 *
 * It is
 *
 * - not visible to the user,
 * - is added automatically, and
 * - does not have a corresponding "keyword" in ocean::Factory.
 */
class InitializationHelper : public OceanModel {
public:
  InitializationHelper(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> in);

private:
  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;

  void update_impl(const Geometry &geometry, double t, double dt);
  void init_impl(const Geometry &geometry);

  const IceModelVec2S& shelf_base_temperature_impl() const;
  const IceModelVec2S& shelf_base_mass_flux_impl() const;
  const IceModelVec2S& melange_back_pressure_fraction_impl() const;

  // storage for melange_back_pressure_fraction is inherited from OceanModel
  IceModelVec2S::Ptr m_shelf_base_temperature;
  IceModelVec2S::Ptr m_shelf_base_mass_flux;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* POINITIALIZATION_H */
