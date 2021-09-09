// Copyright (C) 2008-2020 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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
#include "pism/util/io/File.hh"
#include "pism/util/Vars.hh"
#include "pism/util/Time.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace surface {

IceModelVec2S::Ptr SurfaceModel::allocate_layer_mass(IceGrid::ConstPtr grid) {
  IceModelVec2S::Ptr result(new IceModelVec2S(grid, "surface_layer_mass", WITHOUT_GHOSTS));

  result->set_attrs("climate_forcing", "mass held in the surface layer",
                    "kg", "kg", "", 0);

  result->metadata()["valid_min"] = {0.0};

  return result;
}

IceModelVec2S::Ptr SurfaceModel::allocate_layer_thickness(IceGrid::ConstPtr grid) {

  IceModelVec2S::Ptr result(new IceModelVec2S(grid, "surface_layer_thickness", WITHOUT_GHOSTS));

  result->set_attrs("climate_forcing",
                    "thickness of the surface process layer at the top surface of the ice",
                    "m", "m", "", 0);

  result->metadata()["valid_min"] = {0.0};

  return result;
}

IceModelVec2S::Ptr SurfaceModel::allocate_liquid_water_fraction(IceGrid::ConstPtr grid) {

  IceModelVec2S::Ptr result(new IceModelVec2S(grid,
                                              "ice_surface_liquid_water_fraction", WITHOUT_GHOSTS));

  result->set_attrs("climate_forcing",
                    "liquid water fraction of the ice at the top surface",
                    "1", "1", "", 0);

  result->metadata()["valid_range"] = {0.0, 1.0};

  return result;
}

IceModelVec2S::Ptr SurfaceModel::allocate_mass_flux(IceGrid::ConstPtr grid) {

  IceModelVec2S::Ptr result(new IceModelVec2S(grid, "climatic_mass_balance", WITHOUT_GHOSTS));

  result->set_attrs("climate_forcing",
                    "surface mass balance (accumulation/ablation) rate",
                    "kg m-2 second-1", "kg m-2 year-1",
                    "land_ice_surface_specific_mass_balance_flux", 0);

  Config::ConstPtr config = grid->ctx()->config();
  const double smb_max = config->get_number("surface.given.smb_max", "kg m-2 second-1");

  result->metadata()["valid_range"] = {-smb_max, smb_max};

  return result;
}

IceModelVec2S::Ptr SurfaceModel::allocate_temperature(IceGrid::ConstPtr grid) {

  IceModelVec2S::Ptr result(new IceModelVec2S(grid, "ice_surface_temp", WITHOUT_GHOSTS));

  result->set_attrs("climate_forcing",
                    "temperature of the ice at the ice surface but below firn processes",
                    "Kelvin", "Kelvin", "", 0);

  result->metadata()["valid_range"] = {0.0, 323.15}; // [0C, 50C]

  return result;
}

IceModelVec2S::Ptr SurfaceModel::allocate_accumulation(IceGrid::ConstPtr grid) {

  IceModelVec2S::Ptr result(new IceModelVec2S(grid, "surface_accumulation_flux", WITHOUT_GHOSTS));

  result->set_attrs("diagnostic",
                    "surface accumulation (precipitation minus rain)",
                    "kg m-2", "kg m-2", "", 0);

  return result;
}

IceModelVec2S::Ptr SurfaceModel::allocate_melt(IceGrid::ConstPtr grid) {

  IceModelVec2S::Ptr result(new IceModelVec2S(grid, "surface_melt_flux", WITHOUT_GHOSTS));

  result->set_attrs("diagnostic",
                    "surface melt",
                    "kg m-2", "kg m-2", "", 0);

  return result;
}

IceModelVec2S::Ptr SurfaceModel::allocate_runoff(IceGrid::ConstPtr grid) {

  IceModelVec2S::Ptr result(new IceModelVec2S(grid, "surface_runoff_flux", WITHOUT_GHOSTS));

  result->set_attrs("diagnostic",
                    "surface meltwater runoff",
                    "kg m-2", "kg m-2", "", 0);

  return result;
}

SurfaceModel::SurfaceModel(IceGrid::ConstPtr grid)
  : Component(grid) {

  m_liquid_water_fraction = allocate_liquid_water_fraction(grid);
  m_layer_mass            = allocate_layer_mass(grid);
  m_layer_thickness       = allocate_layer_thickness(grid);
  m_accumulation          = allocate_accumulation(grid);
  m_melt                  = allocate_melt(grid);
  m_runoff                = allocate_runoff(grid);

  // default values
  m_layer_thickness->set(0.0);
  m_layer_mass->set(0.0);
  m_liquid_water_fraction->set(0.0);
  m_accumulation->set(0.0);
  m_melt->set(0.0);
  m_runoff->set(0.0);
}

SurfaceModel::SurfaceModel(IceGrid::ConstPtr g, std::shared_ptr<SurfaceModel> input)
  : Component(g) {
  m_input_model = input;
  // this is a modifier: allocate storage only if necessary (in derived classes)
}

SurfaceModel::SurfaceModel(IceGrid::ConstPtr grid,
                           std::shared_ptr<atmosphere::AtmosphereModel> atmosphere)
  : SurfaceModel(grid) {        // this constructor will allocate storage

  m_atmosphere = atmosphere;
}

SurfaceModel::~SurfaceModel() {
  // empty
}


//! \brief Returns accumulation
/*!
 * Basic surface models currently implemented in PISM do not model accumulation
 */
const IceModelVec2S& SurfaceModel::accumulation() const {
  return accumulation_impl();
}

//! \brief Returns melt
/*!
 * Basic surface models currently implemented in PISM do not model melt
 */
const IceModelVec2S& SurfaceModel::melt() const {
  return melt_impl();
}

//! \brief Returns runoff
/*!
 * Basic surface models currently implemented in PISM do not model runoff
 */
const IceModelVec2S& SurfaceModel::runoff() const {
  return runoff_impl();
}

const IceModelVec2S& SurfaceModel::mass_flux() const {
  return mass_flux_impl();
}

const IceModelVec2S& SurfaceModel::temperature() const {
  return temperature_impl();
}

//! \brief Returns the liquid water fraction of the ice at the top ice surface.
/*!
 * Most PISM surface models return 0.
 */
const IceModelVec2S& SurfaceModel::liquid_water_fraction() const {
  return liquid_water_fraction_impl();
}

//! \brief Returns mass held in the surface layer.
/*!
 * Basic surface models currently implemented in PISM do not model the mass of
 * the surface layer.
 */
const IceModelVec2S& SurfaceModel::layer_mass() const {
  return layer_mass_impl();
}

//! \brief Returns thickness of the surface layer. Could be used to compute surface
//! elevation as a sum of elevation of the top surface of the ice and surface layer (firn,
//! etc) thickness.
/*!
 * Basic surface models currently implemented in PISM do not model surface
 * layer thickness.
 */
const IceModelVec2S& SurfaceModel::layer_thickness() const {
  return layer_thickness_impl();
}

const IceModelVec2S& SurfaceModel::accumulation_impl() const {
  if (m_input_model) {
    return m_input_model->accumulation();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

const IceModelVec2S& SurfaceModel::melt_impl() const {
  if (m_input_model) {
    return m_input_model->melt();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

const IceModelVec2S& SurfaceModel::runoff_impl() const {
  if (m_input_model) {
    return m_input_model->runoff();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

const IceModelVec2S& SurfaceModel::mass_flux_impl() const {
  if (m_input_model) {
    return m_input_model->mass_flux();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

const IceModelVec2S& SurfaceModel::temperature_impl() const {
  if (m_input_model) {
    return m_input_model->temperature();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

const IceModelVec2S& SurfaceModel::liquid_water_fraction_impl() const {
  if (m_input_model) {
    return m_input_model->liquid_water_fraction();
  } else {
    return *m_liquid_water_fraction;
  }
}

const IceModelVec2S& SurfaceModel::layer_mass_impl() const {
  if (m_input_model) {
    return m_input_model->layer_mass();
  } else {
    return *m_layer_mass;
  }
}

const IceModelVec2S& SurfaceModel::layer_thickness_impl() const {
  if (m_input_model) {
    return m_input_model->layer_thickness();
  } else {
    return *m_layer_thickness;
  }
}

TSDiagnosticList SurfaceModel::ts_diagnostics_impl() const {
  if (m_atmosphere) {
    return m_atmosphere->ts_diagnostics();
  }

  if (m_input_model) {
    return m_input_model->ts_diagnostics();
  }

  return {};
}

void SurfaceModel::init(const Geometry &geometry) {
  this->init_impl(geometry);
}

void SurfaceModel::init_impl(const Geometry &geometry) {
  if (m_atmosphere) {
    m_atmosphere->init(geometry);
  }

  if (m_input_model) {
    m_input_model->init(geometry);
  }
}

void SurfaceModel::update(const Geometry &geometry, double t, double dt) {
  this->update_impl(geometry, t, dt);
}

void SurfaceModel::update_impl(const Geometry &geometry, double t, double dt) {
  if (m_atmosphere) {
    m_atmosphere->update(geometry, t, dt);
  }

  if (m_input_model) {
    m_input_model->update(geometry, t, dt);
  }
}

void SurfaceModel::define_model_state_impl(const File &output) const {
  if (m_atmosphere) {
    m_atmosphere->define_model_state(output);
  }

  if (m_input_model) {
    m_input_model->define_model_state(output);
  }
}

void SurfaceModel::write_model_state_impl(const File &output) const {
  if (m_atmosphere) {
    m_atmosphere->write_model_state(output);
  }

  if (m_input_model) {
    m_input_model->write_model_state(output);
  }
}

MaxTimestep SurfaceModel::max_timestep_impl(double t) const {
  if (m_atmosphere) {
    return m_atmosphere->max_timestep(t);
  }

  if (m_input_model) {
    return m_input_model->max_timestep(t);
  }

  return MaxTimestep("surface model");
}

/*!
 * Use the surface mass balance to compute dummy accumulation.
 *
 * This is used by surface models that compute the SMB but do not provide accumulation,
 * melt, and runoff.
 *
 * We assume that the positive part of the SMB is accumulation and the negative part is
 * runoff. This ensures that outputs of PISM's surface models satisfy "SMB = accumulation
 * - runoff".
 */
void SurfaceModel::dummy_accumulation(const IceModelVec2S& smb, IceModelVec2S& result) {

  IceModelVec::AccessList list{&result, &smb};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    result(i,j) = std::max(smb(i,j), 0.0);
  }
}

/*!
 * Use the surface mass balance to compute dummy runoff.
 *
 * This is used by surface models that compute the SMB but do not provide accumulation,
 * melt, and runoff.
 *
 * We assume that the positive part of the SMB is accumulation and the negative part is
 * runoff. This ensures that outputs of PISM's surface models satisfy "SMB = accumulation
 * - runoff".
 */
void SurfaceModel::dummy_runoff(const IceModelVec2S& smb, IceModelVec2S& result) {

  IceModelVec::AccessList list{&result, &smb};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    result(i,j) = std::max(-smb(i,j), 0.0);
  }
}

/*!
 * Use the surface mass balance to compute dummy runoff.
 *
 * This is used by surface models that compute the SMB but do not provide accumulation,
 * melt, and runoff.
 *
 * We assume that all melt runs off, i.e. runoff = melt, but treat melt as a "derived"
 * quantity.
 */
void SurfaceModel::dummy_melt(const IceModelVec2S& smb, IceModelVec2S& result) {
  dummy_runoff(smb, result);
}

namespace diagnostics {

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

  result->copy_from(model->mass_flux());

  return result;
}

PS_ice_surface_temp::PS_ice_surface_temp(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {


  auto ismip6 = m_config->get_flag("output.ISMIP6");

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys,
                                    ismip6 ? "litemptop" : "ice_surface_temp")};

  set_attrs("ice temperature at the top ice surface",
            "temperature_at_top_of_ice_sheet_model",
            "K", "K", 0);
}

IceModelVec::Ptr PS_ice_surface_temp::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "ice_surface_temp", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->temperature());

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

  result->copy_from(model->liquid_water_fraction());

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

  result->copy_from(model->layer_mass());

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

  result->copy_from(model->layer_thickness());

  return result;
}
} // end of namespace diagnostics

DiagnosticList SurfaceModel::diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    {"climatic_mass_balance",             Diagnostic::Ptr(new PS_climatic_mass_balance(this))},
    {"ice_surface_temp",                  Diagnostic::Ptr(new PS_ice_surface_temp(this))},
    {"ice_surface_liquid_water_fraction", Diagnostic::Ptr(new PS_liquid_water_fraction(this))},
    {"surface_layer_mass",                Diagnostic::Ptr(new PS_layer_mass(this))},
    {"surface_layer_thickness",           Diagnostic::Ptr(new PS_layer_thickness(this))}
  };

  if (m_config->get_flag("output.ISMIP6")) {
    result["litemptop"] = Diagnostic::Ptr(new PS_ice_surface_temp(this));
  }

  if (m_atmosphere) {
    result = pism::combine(result, m_atmosphere->diagnostics());
  }

  if (m_input_model) {
    result = pism::combine(result, m_input_model->diagnostics());
  }

  return result;
}

} // end of namespace surface
} // end of namespace pism

