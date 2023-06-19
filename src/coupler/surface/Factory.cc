/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2023 PISM Authors
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

#include "pism/coupler/surface/Factory.hh"

// surface models:
#include "pism/coupler/surface/Anomaly.hh"
#include "pism/coupler/surface/Cache.hh"
#include "pism/coupler/surface/ConstantPIK.hh"
#include "pism/coupler/surface/Delta_T.hh"
#include "pism/coupler/surface/Elevation.hh"
#include "pism/coupler/surface/ElevationChange.hh"
#include "pism/coupler/surface/ForceThickness.hh"
#include "pism/coupler/surface/GivenClimate.hh"
#include "pism/coupler/surface/ISMIP6Climate.hh"
#include "pism/coupler/surface/NoGLRetreat.hh"
#include "pism/coupler/surface/DEBMSimple.hh"
#include "pism/coupler/surface/Simple.hh"
#include "pism/coupler/surface/TemperatureIndex.hh"

namespace pism {
namespace surface {

Factory::Factory(std::shared_ptr<const Grid> g, std::shared_ptr<atmosphere::AtmosphereModel> input)
  : PCFactory<SurfaceModel>(g, "surface.models"),
    m_input(input) {

  add_surface_model<Elevation>("elevation");
  add_surface_model<Given>("given");
  add_surface_model<ISMIP6>("ismip6");
  add_surface_model<TemperatureIndex>("pdd");
  add_surface_model<PIK>("pik");
  add_surface_model<Simple>("simple");
  add_surface_model<DEBMSimple>("debm_simple");

  add_modifier<Anomaly>("anomaly");
  add_modifier<Cache>("cache");
  add_modifier<Delta_T>("delta_T");
  add_modifier<ForceThickness>("forcing");
  add_modifier<ElevationChange>("elevation_change");
  add_modifier<NoGLRetreat>("no_gl_retreat");
}

std::shared_ptr<SurfaceModel> Factory::create(const std::string &type) {

  std::vector<std::string> choices = split(type, ',');

  // the first element has to be an *actual* model (not a modifier)
  auto j = choices.begin();

  auto result = surface_model(*j, m_input);

  ++j;

  // process remaining arguments:
  for (;j != choices.end(); ++j) {
    result = modifier(*j, result);
  }

  return result;
}

std::shared_ptr<SurfaceModel> Factory::surface_model(const std::string &type,
                                                     std::shared_ptr<InputModel> input) {
  if (m_surface_models.find(type) == m_surface_models.end()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "cannot allocate %s \"%s\".\n"
                                  "Available models:    %s\n",
                                  m_parameter.c_str(), type.c_str(),
                                  key_list(m_surface_models).c_str());
  }

  return m_surface_models[type]->create(m_grid, input);
}

} // end of namespace surface
} // end of namespace pism
