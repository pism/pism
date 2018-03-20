// Copyright (C) 2008-2018 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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
#include <gsl/gsl_math.h>       // GSL_NAN

#include "pism/coupler/SurfaceModel.hh"
#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/Vars.hh"
#include "pism/util/Time.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace surface {

// SurfaceModel diagnostics (these don't need to be in the header)

/*! @brief Climatic mass balance */
class PS_climatic_mass_balance : public Diag<SurfaceModel>
{
public:
  PS_climatic_mass_balance(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

/*! @brief Ice surface temperature. */
class PS_ice_surface_temp : public Diag<SurfaceModel>
{
public:
  PS_ice_surface_temp(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

/*! @brief Ice liquid water fraction at the ice surface. */
class PS_liquid_water_fraction : public Diag<SurfaceModel>
{
public:
  PS_liquid_water_fraction(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

/*! @brief Mass of the surface layer (snow and firn). */
class PS_layer_mass : public Diag<SurfaceModel>
{
public:
  PS_layer_mass(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

/*! @brief Surface layer (snow and firn) thickness. */
class PS_layer_thickness : public Diag<SurfaceModel>
{
public:
  PS_layer_thickness(const SurfaceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

///// PISMSurfaceModel base class:

SurfaceModel::SurfaceModel(IceGrid::ConstPtr g)
  : Component(g) {
  m_atmosphere = NULL;
}

SurfaceModel::~SurfaceModel() {
  // empty
}

void SurfaceModel::mass_flux(IceModelVec2S &result) const {
  this->mass_flux_impl(result);
}

DiagnosticList SurfaceModel::diagnostics_impl() const {
  DiagnosticList result = {
    {"climatic_mass_balance",             Diagnostic::Ptr(new PS_climatic_mass_balance(this))},
    {"ice_surface_temp",                  Diagnostic::Ptr(new PS_ice_surface_temp(this))},
    {"ice_surface_liquid_water_fraction", Diagnostic::Ptr(new PS_liquid_water_fraction(this))},
    {"surface_layer_mass",                Diagnostic::Ptr(new PS_layer_mass(this))},
    {"surface_layer_thickness",           Diagnostic::Ptr(new PS_layer_thickness(this))}
  };

  if (m_atmosphere) {
    result = pism::combine(result, m_atmosphere->diagnostics());
  }

  return result;
}

TSDiagnosticList SurfaceModel::ts_diagnostics_impl() const {
  if (m_atmosphere) {
    return m_atmosphere->ts_diagnostics();
  } else {
    return {};
  }
}

void SurfaceModel::attach_atmosphere_model(std::shared_ptr<atmosphere::AtmosphereModel> input) {
  this->attach_atmosphere_model_impl(input);
}

void SurfaceModel::attach_atmosphere_model_impl(std::shared_ptr<atmosphere::AtmosphereModel> input) {
  m_atmosphere = input;
}

void SurfaceModel::init(const Geometry &geometry) {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock
  this->init_impl(geometry);
}

void SurfaceModel::init_impl(const Geometry &geometry) {
  assert(m_atmosphere != NULL);
  m_atmosphere->init(geometry);
}

void SurfaceModel::update(const Geometry &geometry, double t, double dt) {
  this->update_impl(geometry, t, dt);
}

//! \brief Returns mass held in the surface layer.
/*!
 * Basic surface models currently implemented in PISM do not model the mass of
 * the surface layer.
 */
void SurfaceModel::layer_mass(IceModelVec2S &result) const {
  this->layer_mass_impl(result);
}

void SurfaceModel::layer_mass_impl(IceModelVec2S &result) const {
  result.set(0.0);
}

//! \brief Returns thickness of the surface layer. Used to compute surface
//! elevation as a sum of elevation of the top surface of the ice and surface
//! layer (firn, etc) thickness.
/*!
 * Basic surface models currently implemented in PISM do not model surface
 * layer thickness.
 */
void SurfaceModel::layer_thickness(IceModelVec2S &result) const {
  this->layer_thickness_impl(result);
}

void SurfaceModel::layer_thickness_impl(IceModelVec2S &result) const {
  result.set(0.0);
}

void SurfaceModel::temperature(IceModelVec2S &result) const {
  this->temperature_impl(result);
}
//! \brief Returns the liquid water fraction of the ice at the top ice surface.
/*!
 * Most PISM surface models return 0.
 */
void SurfaceModel::liquid_water_fraction(IceModelVec2S &result) const {
  this->liquid_water_fraction_impl(result);
}

void SurfaceModel::liquid_water_fraction_impl(IceModelVec2S &result) const {
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

IceModelVec2S::Ptr SurfaceModel::mass_flux() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "mass_flux", WITHOUT_GHOSTS));
  mass_flux(*result);
  return result;
}

IceModelVec2S::Ptr SurfaceModel::temperature() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "temperature", WITHOUT_GHOSTS));
  temperature(*result);
  return result;
}

IceModelVec2S::Ptr SurfaceModel::liquid_water_fraction() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "liquid_water_fraction", WITHOUT_GHOSTS));
  liquid_water_fraction(*result);
  return result;
}

IceModelVec2S::Ptr SurfaceModel::layer_mass() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "layer_mass", WITHOUT_GHOSTS));
  layer_mass(*result);
  return result;
}

IceModelVec2S::Ptr SurfaceModel::layer_thickness() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "layer_thickness", WITHOUT_GHOSTS));
  layer_thickness(*result);
  return result;
}

PS_climatic_mass_balance::PS_climatic_mass_balance(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "climatic_mass_balance")};

  set_attrs("surface mass balance (accumulation/ablation) rate",
            "land_ice_surface_specific_mass_balance_flux",
            "kg m-2 second-1", "kg m-2 year-1", 0);
}

IceModelVec::Ptr PS_climatic_mass_balance::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  model->mass_flux(*result);

  return result;
}

PS_ice_surface_temp::PS_ice_surface_temp(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "ice_surface_temp")};

  set_attrs("ice temperature at the ice surface", "",
            "Kelvin", "Kelvin", 0);
}

IceModelVec::Ptr PS_ice_surface_temp::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "ice_surface_temp", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  model->temperature(*result);

  return result;
}

PS_liquid_water_fraction::PS_liquid_water_fraction(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "ice_surface_liquid_water_fraction")};

  set_attrs("ice liquid water fraction at the ice surface", "",
            "1", "1", 0);
}

IceModelVec::Ptr PS_liquid_water_fraction::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "ice_surface_liquid_water_fraction", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  model->liquid_water_fraction(*result);

  return result;
}

PS_layer_mass::PS_layer_mass(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "surface_layer_mass")};

  set_attrs("mass of the surface layer (snow and firn)", "",
            "kg", "kg", 0);
}

IceModelVec::Ptr PS_layer_mass::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "surface_layer_mass", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  model->layer_mass(*result);

  return result;
}

PS_layer_thickness::PS_layer_thickness(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "surface_layer_thickness")};

  set_attrs("thickness of the surface layer (snow and firn)", "",
            "meters", "meters", 0);
}

IceModelVec::Ptr PS_layer_thickness::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "surface_layer_thickness", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  model->layer_thickness(*result);

  return result;
}

} // end of namespace surface
} // end of namespace pism

