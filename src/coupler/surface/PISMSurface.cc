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

#include "PISMSurface.hh"
#include "PISMAtmosphere.hh"
#include "PIO.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "iceModelVec.hh"
#include <assert.h>

namespace pism {

///// PISMSurfaceModel base class:

SurfaceModel::SurfaceModel(IceGrid &g, const Config &conf)
  : Component_TS(g, conf) {
  atmosphere = NULL;
}

SurfaceModel::~SurfaceModel() {
  delete atmosphere;
}

void SurfaceModel::get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                                       std::map<std::string, TSDiagnostic*> &ts_dict) {
  if (atmosphere) {
    atmosphere->get_diagnostics(dict, ts_dict);
  }
}

void SurfaceModel::attach_atmosphere_model(AtmosphereModel *input) {
  if (atmosphere != NULL) {
    delete atmosphere;
  }
  atmosphere = input;
}

void SurfaceModel::init(Vars &vars) {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  assert(atmosphere != NULL);
  atmosphere->init(vars);
}

//! \brief Returns mass held in the surface layer.
/*!
 * Basic surface models currently implemented in PISM do not model the mass of
 * the surface layer.
 */
void SurfaceModel::mass_held_in_surface_layer(IceModelVec2S &result) {
  result.set(0.0);
}

//! \brief Returns thickness of the surface layer. Used to compute surface
//! elevation as a sum of elevation of the top surface of the ice and surface
//! layer (firn, etc) thickness.
/*!
 * Basic surface models currently implemented in PISM do not model surface
 * layer thickness.
 */
void SurfaceModel::surface_layer_thickness(IceModelVec2S &result) {
  result.set(0.0);
}

//! \brief Returns the liquid water fraction of the ice at the top ice surface.
/*!
 * Most PISM surface models return 0.
 */
void SurfaceModel::ice_surface_liquid_water_fraction(IceModelVec2S &result) {
  result.set(0.0);
}

void SurfaceModel::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  if (atmosphere != NULL) {
    atmosphere->define_variables(vars, nc, nctype);
  }
}

void SurfaceModel::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  if (atmosphere != NULL) {
    atmosphere->write_variables(vars, nc);
  }
}

void SurfaceModel::max_timestep(double my_t, double &my_dt, bool &restrict) {
  if (atmosphere != NULL) {
    atmosphere->max_timestep(my_t, my_dt, restrict);
  } else {
    my_dt = -1;
    restrict = false;
  }
}

void SurfaceModel::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (atmosphere != NULL) {
    atmosphere->add_vars_to_output(keyword, result);
  }
}

} // end of namespace pism
