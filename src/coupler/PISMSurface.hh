// Copyright (C) 2008-2016 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

#ifndef __PISMSurfaceModel_hh
#define __PISMSurfaceModel_hh

/*!
 * This file should contain the class definition and nothing else.
 * Implementations should go in separate files.
 */

#include "base/util/PISMComponent.hh"

namespace pism {

namespace atmosphere {
class AtmosphereModel;
}

class IceModelVec2S;

//! @brief Surface models and modifiers: provide top-surface
//! temperature, mass flux, liquid water fraction, mass and thickness of the surface
//! layer.
namespace surface {

//! \brief The interface of PISM's surface models.
class SurfaceModel : public Component_TS {
public:
  SurfaceModel(IceGrid::ConstPtr g);
  virtual ~SurfaceModel();

  // the interface:
  void ice_surface_mass_flux(IceModelVec2S &result);

  void ice_surface_temperature(IceModelVec2S &result);
  void ice_surface_liquid_water_fraction(IceModelVec2S &result);

  void mass_held_in_surface_layer(IceModelVec2S &result);
  void surface_layer_thickness(IceModelVec2S &result);

  void attach_atmosphere_model(atmosphere::AtmosphereModel *input);

  void init();
protected:
  virtual void init_impl();
  virtual MaxTimestep max_timestep_impl(double my_t);

  virtual void attach_atmosphere_model_impl(atmosphere::AtmosphereModel *input);

  virtual void surface_layer_thickness_impl(IceModelVec2S &result);
  virtual void mass_held_in_surface_layer_impl(IceModelVec2S &result);

  virtual void ice_surface_temperature_impl(IceModelVec2S &result) = 0;  
  virtual void ice_surface_liquid_water_fraction_impl(IceModelVec2S &result);

  virtual void ice_surface_mass_flux_impl(IceModelVec2S &result) = 0;

  virtual void get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                    std::map<std::string, TSDiagnostic::Ptr> &ts_dict);
  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars,
                                     const PIO &nc, IO_Type nctype);
protected:
  atmosphere::AtmosphereModel *m_atmosphere;
};

} // end of namespace surface
} // end of namespace pism

#endif  // __PISMSurfaceModel_hh

