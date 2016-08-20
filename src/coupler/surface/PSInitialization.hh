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

#ifndef PSINITIALIZATION_H
#define PSINITIALIZATION_H

#include "PSModifier.hh"

namespace pism {
namespace surface {

class Initialization : public SurfaceModifier {
public:
  Initialization(IceGrid::ConstPtr g, SurfaceModel* in);
protected:
  void attach_atmosphere_model_impl(atmosphere::AtmosphereModel *in);

  void update_impl(double my_t, double my_dt);
  void init_impl();

  void ice_surface_mass_flux_impl(IceModelVec2S &result);
  void ice_surface_temperature_impl(IceModelVec2S &result);
  void ice_surface_liquid_water_fraction_impl(IceModelVec2S &result);
  void mass_held_in_surface_layer_impl(IceModelVec2S &result);
  void surface_layer_thickness_impl(IceModelVec2S &result);

  void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  void define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                             IO_Type nctype);

private:
  bool m_update_called;
  // store pointers to fields to make it possible to iterate over them
  std::vector<const IceModelVec*> m_variables;
  // storage
  IceModelVec2S m_ice_surface_mass_flux;
  IceModelVec2S m_ice_surface_temperature;
  IceModelVec2S m_ice_surface_liquid_water_fraction;
  IceModelVec2S m_mass_held_in_surface_layer;
  IceModelVec2S m_surface_layer_thickness;
};

} // end of namespace surface
} // end of namespace pism


#endif /* PSINITIALIZATION_H */
