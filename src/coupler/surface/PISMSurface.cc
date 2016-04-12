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

#include <cassert>
#include <gsl/gsl_math.h>

#include "coupler/PISMSurface.hh"
#include "coupler/PISMAtmosphere.hh"
#include "base/util/io/PIO.hh"
#include "base/util/PISMVars.hh"
#include "base/util/PISMTime.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_options.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/MaxTimestep.hh"

namespace pism {
namespace surface {

///// PISMSurfaceModel base class:

SurfaceModel::SurfaceModel(IceGrid::ConstPtr g)
  : Component_TS(g) {
  m_atmosphere = NULL;
}

SurfaceModel::~SurfaceModel() {
  delete m_atmosphere;
}

void SurfaceModel::ice_surface_mass_flux(IceModelVec2S &result) {
  this->ice_surface_mass_flux_impl(result);
}

void SurfaceModel::get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                       std::map<std::string, TSDiagnostic::Ptr> &ts_dict) {
  if (m_atmosphere) {
    m_atmosphere->get_diagnostics(dict, ts_dict);
  }
}

void SurfaceModel::attach_atmosphere_model(atmosphere::AtmosphereModel *input) {
  this->attach_atmosphere_model_impl(input);
}

void SurfaceModel::attach_atmosphere_model_impl(atmosphere::AtmosphereModel *input) {
  if (m_atmosphere != NULL) {
    delete m_atmosphere;
  }
  m_atmosphere = input;
}

void SurfaceModel::init() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock
  this->init_impl();
}

void SurfaceModel::init_impl() {
  assert(m_atmosphere != NULL);
  m_atmosphere->init();
}

//! \brief Returns mass held in the surface layer.
/*!
 * Basic surface models currently implemented in PISM do not model the mass of
 * the surface layer.
 */
void SurfaceModel::mass_held_in_surface_layer(IceModelVec2S &result) {
  this->mass_held_in_surface_layer_impl(result);
}

void SurfaceModel::mass_held_in_surface_layer_impl(IceModelVec2S &result) {
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
  this->surface_layer_thickness_impl(result);
}

void SurfaceModel::surface_layer_thickness_impl(IceModelVec2S &result) {
  result.set(0.0);
}

void SurfaceModel::ice_surface_temperature(IceModelVec2S &result) {
  this->ice_surface_temperature_impl(result);
}
//! \brief Returns the liquid water fraction of the ice at the top ice surface.
/*!
 * Most PISM surface models return 0.
 */
void SurfaceModel::ice_surface_liquid_water_fraction(IceModelVec2S &result) {
  this->ice_surface_liquid_water_fraction_impl(result);
}

void SurfaceModel::ice_surface_liquid_water_fraction_impl(IceModelVec2S &result) {
  result.set(0.0);
}

void SurfaceModel::define_variables_impl(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  if (m_atmosphere != NULL) {
    m_atmosphere->define_variables(vars, nc, nctype);
  }
}

void SurfaceModel::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {
  if (m_atmosphere != NULL) {
    m_atmosphere->write_variables(vars, nc);
  }
}

MaxTimestep SurfaceModel::max_timestep_impl(double my_t) {
  if (m_atmosphere != NULL) {
    return m_atmosphere->max_timestep(my_t);
  } else {
    return MaxTimestep();
  }
}

void SurfaceModel::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  if (m_atmosphere != NULL) {
    m_atmosphere->add_vars_to_output(keyword, result);
  }
}

} // end of namespace surface
} // end of namespace pism

