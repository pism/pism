/* Copyright (C) 2016 PISM Authors
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

#include "POModifier.hh"

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
class InitializationHelper : public OceanModifier {
public:
  InitializationHelper(IceGrid::ConstPtr g, OceanModel* in);
protected:
  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;

  void update_impl(double t, double dt);
  void init_impl();

  void melange_back_pressure_fraction_impl(IceModelVec2S &result) const;
  void sea_level_elevation_impl(double &result) const;
  void shelf_base_temperature_impl(IceModelVec2S &result) const;
  void shelf_base_mass_flux_impl(IceModelVec2S &result) const;

private:
  IceModelVec2S m_melange_back_pressure_fraction;
  IceModelVec2S m_shelf_base_temperature;
  IceModelVec2S m_shelf_base_mass_flux;

  double m_sea_level_elevation;

  TimeseriesMetadata m_sea_level_metadata;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* POINITIALIZATION_H */
