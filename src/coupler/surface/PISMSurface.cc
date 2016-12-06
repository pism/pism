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
#include "base/util/pism_utilities.hh"

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

void SurfaceModel::ice_surface_mass_flux(IceModelVec2S &result) const {
  this->ice_surface_mass_flux_impl(result);
}

std::map<std::string, Diagnostic::Ptr> SurfaceModel::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"climatic_mass_balance",             Diagnostic::Ptr(new PS_climatic_mass_balance(this))},
    {"ice_surface_temp",                  Diagnostic::Ptr(new PS_ice_surface_temp(this))},
    {"ice_surface_liquid_water_fraction", Diagnostic::Ptr(new PS_liquid_water_fraction(this))},
    {"surface_layer_mass",                Diagnostic::Ptr(new PS_surface_layer_mass(this))},
    {"surface_layer_thickness",           Diagnostic::Ptr(new PS_surface_layer_thickness(this))}
  };

  if (m_atmosphere) {
    result = pism::combine(result, m_atmosphere->diagnostics());
  }

  return result;
}

std::map<std::string, TSDiagnostic::Ptr> SurfaceModel::ts_diagnostics_impl() const {
  if (m_atmosphere) {
    return m_atmosphere->ts_diagnostics();
  } else {
    return {};
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
void SurfaceModel::mass_held_in_surface_layer(IceModelVec2S &result) const {
  this->mass_held_in_surface_layer_impl(result);
}

void SurfaceModel::mass_held_in_surface_layer_impl(IceModelVec2S &result) const {
  result.set(0.0);
}

//! \brief Returns thickness of the surface layer. Used to compute surface
//! elevation as a sum of elevation of the top surface of the ice and surface
//! layer (firn, etc) thickness.
/*!
 * Basic surface models currently implemented in PISM do not model surface
 * layer thickness.
 */
void SurfaceModel::surface_layer_thickness(IceModelVec2S &result) const {
  this->surface_layer_thickness_impl(result);
}

void SurfaceModel::surface_layer_thickness_impl(IceModelVec2S &result) const {
  result.set(0.0);
}

void SurfaceModel::ice_surface_temperature(IceModelVec2S &result) const {
  this->ice_surface_temperature_impl(result);
}
//! \brief Returns the liquid water fraction of the ice at the top ice surface.
/*!
 * Most PISM surface models return 0.
 */
void SurfaceModel::ice_surface_liquid_water_fraction(IceModelVec2S &result) const {
  this->ice_surface_liquid_water_fraction_impl(result);
}

void SurfaceModel::ice_surface_liquid_water_fraction_impl(IceModelVec2S &result) const {
  result.set(0.0);
}

void SurfaceModel::define_model_state_impl(const PIO &output) const {
  if (m_atmosphere != NULL) {
    m_atmosphere->define_model_state(output);
  }
}

void SurfaceModel::write_model_state_impl(const PIO &output) const {
  if (m_atmosphere != NULL) {
    m_atmosphere->write_model_state(output);
  }
}

MaxTimestep SurfaceModel::max_timestep_impl(double my_t) const {
  if (m_atmosphere != NULL) {
    return m_atmosphere->max_timestep(my_t);
  } else {
    return MaxTimestep("surface model");
  }
}

PS_climatic_mass_balance::PS_climatic_mass_balance(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "climatic_mass_balance")};

  set_attrs("surface mass balance (accumulation/ablation) rate",
            "land_ice_surface_specific_mass_balance_flux",
            "kg m-2 second-1", "kg m-2 year-1", 0);
}

IceModelVec::Ptr PS_climatic_mass_balance::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->ice_surface_mass_flux(*result);

  return result;
}

PS_ice_surface_temp::PS_ice_surface_temp(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "ice_surface_temp")};

  set_attrs("ice temperature at the ice surface", "",
            "Kelvin", "Kelvin", 0);
}

IceModelVec::Ptr PS_ice_surface_temp::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->ice_surface_temperature(*result);

  return result;
}

PS_liquid_water_fraction::PS_liquid_water_fraction(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "ice_surface_liquid_water_fraction")};

  set_attrs("ice liquid water fraction at the ice surface", "",
            "1", "1", 0);
}

IceModelVec::Ptr PS_liquid_water_fraction::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "ice_surface_liquid_water_fraction", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->ice_surface_liquid_water_fraction(*result);

  return result;
}

PS_surface_layer_mass::PS_surface_layer_mass(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "surface_layer_mass")};

  set_attrs("mass of the surface layer (snow and firn)", "",
            "kg", "kg", 0);
}

IceModelVec::Ptr PS_surface_layer_mass::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "surface_layer_mass", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->mass_held_in_surface_layer(*result);

  return result;
}

PS_surface_layer_thickness::PS_surface_layer_thickness(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "surface_layer_thickness")};

  set_attrs("thickness of the surface layer (snow and firn)", "",
            "meters", "meters", 0);
}

IceModelVec::Ptr PS_surface_layer_thickness::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "surface_layer_thickness", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->surface_layer_thickness(*result);

  return result;
}

} // end of namespace surface
} // end of namespace pism

