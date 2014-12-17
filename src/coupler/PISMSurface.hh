// Copyright (C) 2008-2014 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#include "PISMComponent.hh"

namespace pism {

class AtmosphereModel;
class IceModelVec2S;

//! \brief The interface of PISM's surface models.
class SurfaceModel : public Component_TS {
public:
  SurfaceModel(IceGrid &g);
  virtual ~SurfaceModel();

  // the interface:
  virtual void attach_atmosphere_model(AtmosphereModel *input);
  virtual void ice_surface_mass_flux(IceModelVec2S &result) = 0;
  virtual void ice_surface_temperature(IceModelVec2S &result) = 0;
  virtual void ice_surface_liquid_water_fraction(IceModelVec2S &result);
  virtual void mass_held_in_surface_layer(IceModelVec2S &result);
  virtual void surface_layer_thickness(IceModelVec2S &result);

  // provide default re-implementations of these parent's methods:
  virtual void init();
  virtual void get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                               std::map<std::string, TSDiagnostic*> &ts_dict);
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
  virtual void max_timestep(double my_t, double &my_dt, bool &restrict);
protected:
  AtmosphereModel *atmosphere;
};

} // end of namespace pism

#endif  // __PISMSurfaceModel_hh

