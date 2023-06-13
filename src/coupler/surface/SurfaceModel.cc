// Copyright (C) 2008-2023 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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
#include "pism/util/IceGrid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace surface {

array::Scalar::Ptr SurfaceModel::allocate_layer_mass(IceGrid::ConstPtr grid) {
  array::Scalar::Ptr result(new array::Scalar(grid, "surface_layer_mass"));

  result->set_attrs("climate_forcing", "mass held in the surface layer",
                    "kg", "kg", "", 0);

  result->metadata()["valid_min"] = {0.0};

  return result;
}

array::Scalar::Ptr SurfaceModel::allocate_layer_thickness(IceGrid::ConstPtr grid) {

  array::Scalar::Ptr result(new array::Scalar(grid, "surface_layer_thickness"));

  result->set_attrs("climate_forcing",
                    "thickness of the surface process layer at the top surface of the ice",
                    "m", "m", "", 0);

  result->metadata()["valid_min"] = {0.0};

  return result;
}

array::Scalar::Ptr SurfaceModel::allocate_liquid_water_fraction(IceGrid::ConstPtr grid) {

  array::Scalar::Ptr result(new array::Scalar(grid,
                                              "ice_surface_liquid_water_fraction"));

  result->set_attrs("climate_forcing",
                    "liquid water fraction of the ice at the top surface",
                    "1", "1", "", 0);

  result->metadata()["valid_range"] = {0.0, 1.0};

  return result;
}

array::Scalar::Ptr SurfaceModel::allocate_mass_flux(IceGrid::ConstPtr grid) {

  array::Scalar::Ptr result(new array::Scalar(grid, "climatic_mass_balance"));

  result->set_attrs("climate_forcing",
                    "surface mass balance (accumulation/ablation) rate",
                    "kg m-2 second-1", "kg m-2 year-1",
                    "land_ice_surface_specific_mass_balance_flux", 0);

  Config::ConstPtr config = grid->ctx()->config();
  const double smb_max = config->get_number("surface.given.smb_max", "kg m-2 second-1");

  result->metadata()["valid_range"] = {-smb_max, smb_max};

  return result;
}

array::Scalar::Ptr SurfaceModel::allocate_temperature(IceGrid::ConstPtr grid) {

  array::Scalar::Ptr result(new array::Scalar(grid, "ice_surface_temp"));

  result->set_attrs("climate_forcing",
                    "temperature of the ice at the ice surface but below firn processes",
                    "Kelvin", "Kelvin", "", 0);

  result->metadata()["valid_range"] = {0.0, 323.15}; // [0C, 50C]

  return result;
}

array::Scalar::Ptr SurfaceModel::allocate_accumulation(IceGrid::ConstPtr grid) {

  array::Scalar::Ptr result(new array::Scalar(grid, "surface_accumulation_flux"));

  result->set_attrs("diagnostic",
                    "surface accumulation (precipitation minus rain)",
                    "kg m-2", "kg m-2", "", 0);

  return result;
}

array::Scalar::Ptr SurfaceModel::allocate_melt(IceGrid::ConstPtr grid) {

  array::Scalar::Ptr result(new array::Scalar(grid, "surface_melt_flux"));

  result->set_attrs("diagnostic",
                    "surface melt",
                    "kg m-2", "kg m-2", "", 0);

  return result;
}

array::Scalar::Ptr SurfaceModel::allocate_runoff(IceGrid::ConstPtr grid) {

  array::Scalar::Ptr result(new array::Scalar(grid, "surface_runoff_flux"));

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


//! \brief Returns accumulation
/*!
 * Basic surface models currently implemented in PISM do not model accumulation
 */
const array::Scalar& SurfaceModel::accumulation() const {
  return accumulation_impl();
}

//! \brief Returns melt
/*!
 * Basic surface models currently implemented in PISM do not model melt
 */
const array::Scalar& SurfaceModel::melt() const {
  return melt_impl();
}

//! \brief Returns runoff
/*!
 * Basic surface models currently implemented in PISM do not model runoff
 */
const array::Scalar& SurfaceModel::runoff() const {
  return runoff_impl();
}

const array::Scalar& SurfaceModel::mass_flux() const {
  return mass_flux_impl();
}

const array::Scalar& SurfaceModel::temperature() const {
  return temperature_impl();
}

//! \brief Returns the liquid water fraction of the ice at the top ice surface.
/*!
 * Most PISM surface models return 0.
 */
const array::Scalar& SurfaceModel::liquid_water_fraction() const {
  return liquid_water_fraction_impl();
}

//! \brief Returns mass held in the surface layer.
/*!
 * Basic surface models currently implemented in PISM do not model the mass of
 * the surface layer.
 */
const array::Scalar& SurfaceModel::layer_mass() const {
  return layer_mass_impl();
}

//! \brief Returns thickness of the surface layer. Could be used to compute surface
//! elevation as a sum of elevation of the top surface of the ice and surface layer (firn,
//! etc) thickness.
/*!
 * Basic surface models currently implemented in PISM do not model surface
 * layer thickness.
 */
const array::Scalar& SurfaceModel::layer_thickness() const {
  return layer_thickness_impl();
}

const array::Scalar& SurfaceModel::accumulation_impl() const {
  if (m_input_model) {
    return m_input_model->accumulation();
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
}

const array::Scalar& SurfaceModel::melt_impl() const {
  if (m_input_model) {
    return m_input_model->melt();
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
}

const array::Scalar& SurfaceModel::runoff_impl() const {
  if (m_input_model) {
    return m_input_model->runoff();
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
}

const array::Scalar& SurfaceModel::mass_flux_impl() const {
  if (m_input_model) {
    return m_input_model->mass_flux();
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
}

const array::Scalar& SurfaceModel::temperature_impl() const {
  if (m_input_model) {
    return m_input_model->temperature();
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
}

const array::Scalar& SurfaceModel::liquid_water_fraction_impl() const {
  if (m_input_model) {
    return m_input_model->liquid_water_fraction();
  }

  return *m_liquid_water_fraction;
}

const array::Scalar& SurfaceModel::layer_mass_impl() const {
  if (m_input_model) {
    return m_input_model->layer_mass();
  }

  return *m_layer_mass;
}

const array::Scalar& SurfaceModel::layer_thickness_impl() const {
  if (m_input_model) {
    return m_input_model->layer_thickness();
  }

  return *m_layer_thickness;
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
void SurfaceModel::dummy_accumulation(const array::Scalar& smb, array::Scalar& result) {

  array::AccessScope list{&result, &smb};

  for (auto p = m_grid->points(); p; p.next()) {
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
void SurfaceModel::dummy_runoff(const array::Scalar& smb, array::Scalar& result) {

  array::AccessScope list{&result, &smb};

  for (auto p = m_grid->points(); p; p.next()) {
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
void SurfaceModel::dummy_melt(const array::Scalar& smb, array::Scalar& result) {
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
  array::Array::Ptr compute_impl() const;
};

/*! @brief Ice surface temperature. */
class PS_ice_surface_temp : public Diag<SurfaceModel>
{
public:
  PS_ice_surface_temp(const SurfaceModel *m);
protected:
  array::Array::Ptr compute_impl() const;
};

/*! @brief Ice liquid water fraction at the ice surface. */
class PS_liquid_water_fraction : public Diag<SurfaceModel>
{
public:
  PS_liquid_water_fraction(const SurfaceModel *m);
protected:
  array::Array::Ptr compute_impl() const;
};

/*! @brief Mass of the surface layer (snow and firn). */
class PS_layer_mass : public Diag<SurfaceModel>
{
public:
  PS_layer_mass(const SurfaceModel *m);
protected:
  array::Array::Ptr compute_impl() const;
};

/*! @brief Surface layer (snow and firn) thickness. */
class PS_layer_thickness : public Diag<SurfaceModel>
{
public:
  PS_layer_thickness(const SurfaceModel *m);
protected:
  array::Array::Ptr compute_impl() const;
};

PS_climatic_mass_balance::PS_climatic_mass_balance(const SurfaceModel *m)
  : Diag<SurfaceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "climatic_mass_balance")};

  set_attrs("surface mass balance (accumulation/ablation) rate",
            "land_ice_surface_specific_mass_balance_flux",
            "kg m-2 second-1", "kg m-2 year-1", 0);
}

array::Array::Ptr PS_climatic_mass_balance::compute_impl() const {

  array::Scalar::Ptr result(new array::Scalar(m_grid, "climatic_mass_balance"));
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

array::Array::Ptr PS_ice_surface_temp::compute_impl() const {

  array::Scalar::Ptr result(new array::Scalar(m_grid, "ice_surface_temp"));
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

array::Array::Ptr PS_liquid_water_fraction::compute_impl() const {

  array::Scalar::Ptr result(new array::Scalar(m_grid, "ice_surface_liquid_water_fraction"));
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

array::Array::Ptr PS_layer_mass::compute_impl() const {

  array::Scalar::Ptr result(new array::Scalar(m_grid, "surface_layer_mass"));
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

array::Array::Ptr PS_layer_thickness::compute_impl() const {

  array::Scalar::Ptr result(new array::Scalar(m_grid, "surface_layer_thickness"));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->layer_thickness());

  return result;
}

enum AmountKind {AMOUNT, MASS};

/*! @brief Report surface melt, averaged over the reporting interval */
class SurfaceMelt : public DiagAverageRate<SurfaceModel>
{
public:
  SurfaceMelt(const SurfaceModel *m, AmountKind kind)
    : DiagAverageRate<SurfaceModel>(m,
                                        kind == AMOUNT
                                        ? "surface_melt_flux"
                                        : "surface_melt_rate",
                                        TOTAL_CHANGE),
      m_melt_mass(m_grid, "melt_mass"),
      m_kind(kind)
  {

    std::string
      name              = "surface_melt_flux",
      long_name         = "surface melt, averaged over the reporting interval",
      standard_name     = "surface_snow_and_ice_melt_flux",
      accumulator_units = "kg m-2",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "surface_melt_rate";
      standard_name     = "";
      accumulator_units = "kg",
      internal_units    = "kg second-1";
      external_units    = "Gt year-1" ;
    }

    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata()["units"] = accumulator_units;

    set_attrs(long_name, standard_name, internal_units, external_units, 0);
    m_vars[0]["cell_methods"] = "time: mean";

    m_vars[0]["_FillValue"] = {to_internal(m_fill_value)};
  }

protected:
  const array::Scalar& model_input() {
    const array::Scalar &melt_amount = model->melt();

    if (m_kind == MASS) {
      double cell_area = m_grid->cell_area();

      array::AccessScope list{&m_melt_mass, &melt_amount};

      for (auto p = m_grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();
        m_melt_mass(i, j) = melt_amount(i, j) * cell_area;
      }
      return m_melt_mass;
    }

    return melt_amount;
  }
private:
  array::Scalar m_melt_mass;
  AmountKind m_kind;
};

/*! @brief Report surface runoff, averaged over the reporting interval */
class SurfaceRunoff : public DiagAverageRate<SurfaceModel>
{
public:
  SurfaceRunoff(const SurfaceModel *m, AmountKind kind)
    : DiagAverageRate<SurfaceModel>(m,
                                        kind == AMOUNT
                                        ? "surface_runoff_flux"
                                        : "surface_runoff_rate",
                                        TOTAL_CHANGE),
      m_kind(kind),
      m_runoff_mass(m_grid, "runoff_mass") {

    std::string
      name              = "surface_runoff_flux",
      long_name         = "surface runoff, averaged over the reporting interval",
      standard_name     = "surface_runoff_flux",
      accumulator_units = "kg m-2",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "surface_runoff_rate";
      standard_name     = "",
      accumulator_units = "kg",
      internal_units    = "kg second-1";
      external_units    = "Gt year-1" ;
    }

    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata()["units"] = accumulator_units;

    set_attrs(long_name, standard_name, internal_units, external_units, 0);
    m_vars[0]["cell_methods"] = "time: mean";

    m_vars[0]["_FillValue"] = {to_internal(m_fill_value)};
  }

protected:
  const array::Scalar& model_input() {
    const array::Scalar &runoff_amount = model->runoff();

    if (m_kind == MASS) {
      double cell_area = m_grid->cell_area();

      array::AccessScope list{&m_runoff_mass, &runoff_amount};

      for (auto p = m_grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();
        m_runoff_mass(i, j) = runoff_amount(i, j) * cell_area;
      }
      return m_runoff_mass;
    }

    return runoff_amount;
  }
private:
  AmountKind m_kind;
  array::Scalar m_runoff_mass;
};

/*! @brief Report accumulation (precipitation minus rain), averaged over the reporting interval */
class Accumulation : public DiagAverageRate<SurfaceModel>
{
public:
  Accumulation(const SurfaceModel *m, AmountKind kind)
    : DiagAverageRate<SurfaceModel>(m,
                                        kind == AMOUNT
                                        ? "surface_accumulation_flux"
                                        : "surface_accumulation_rate",
                                        TOTAL_CHANGE),
      m_kind(kind),
      m_accumulation_mass(m_grid, "accumulation_mass") {

    // possible standard name: surface_accumulation_flux
    std::string
      name              = "surface_accumulation_flux",
      long_name         = "accumulation (precipitation minus rain), averaged over the reporting interval",
      accumulator_units = "kg m-2",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "surface_accumulation_rate";
      accumulator_units = "kg",
      internal_units    = "kg second-1";
      external_units    = "Gt year-1" ;
    }


    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata()["units"] = accumulator_units;

    set_attrs(long_name, "", internal_units, external_units, 0);
    m_vars[0]["cell_methods"] = "time: mean";

    m_vars[0]["_FillValue"] = {to_internal(m_fill_value)};
  }

protected:
  const array::Scalar& model_input() {
    const array::Scalar &accumulation_amount = model->accumulation();

    if (m_kind == MASS) {
      double cell_area = m_grid->cell_area();

      array::AccessScope list{&m_accumulation_mass, &accumulation_amount};

      for (auto p = m_grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();
        m_accumulation_mass(i, j) = accumulation_amount(i, j) * cell_area;
      }
      return m_accumulation_mass;
    }

    return accumulation_amount;
  }
private:
  AmountKind m_kind;
  array::Scalar m_accumulation_mass;
};

/*!
 * Integrate a field over the computational domain.
 *
 * If the input has units kg/m^2, the output will be in kg.
 */
static double integrate(const array::Scalar &input) {
  IceGrid::ConstPtr grid = input.grid();

  double cell_area = grid->cell_area();

  array::AccessScope list{&input};

  double result = 0.0;

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    result += input(i, j) * cell_area;
  }

  return GlobalSum(grid->com, result);
}


//! \brief Reports the total accumulation rate.
class TotalSurfaceAccumulation : public TSDiag<TSFluxDiagnostic, SurfaceModel>
{
public:
  TotalSurfaceAccumulation(const SurfaceModel *m)
    : TSDiag<TSFluxDiagnostic, SurfaceModel>(m, "surface_accumulation_rate") {

    set_units("kg s-1", "kg year-1");
    m_variable["long_name"] = "surface accumulation rate (PDD model)";
  }

  double compute() {
    return integrate(model->accumulation());
  }
};


//! \brief Reports the total melt rate.
class TotalSurfaceMelt : public TSDiag<TSFluxDiagnostic, SurfaceModel>
{
public:
  TotalSurfaceMelt(const SurfaceModel *m)
    : TSDiag<TSFluxDiagnostic, SurfaceModel>(m, "surface_melt_rate") {

    set_units("kg s-1", "kg year-1");
    m_variable["long_name"] = "surface melt rate (PDD model)";
  }

  double compute() {
    return integrate(model->melt());
  }
};


//! \brief Reports the total top surface ice flux.
class TotalSurfaceRunoff : public TSDiag<TSFluxDiagnostic, SurfaceModel>
{
public:
  TotalSurfaceRunoff(const SurfaceModel *m)
    : TSDiag<TSFluxDiagnostic, SurfaceModel>(m, "surface_runoff_rate") {

    set_units("kg s-1", "kg year-1");
    m_variable["long_name"] = "surface runoff rate (PDD model)";
  }

  double compute() {
    return integrate(model->runoff());
  }
};

} // end of namespace diagnostics

DiagnosticList SurfaceModel::diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    {"surface_accumulation_flux",         Diagnostic::Ptr(new Accumulation(this, AMOUNT))},
    {"surface_accumulation_rate",         Diagnostic::Ptr(new Accumulation(this, MASS))},
    {"surface_melt_flux",                 Diagnostic::Ptr(new SurfaceMelt(this, AMOUNT))},
    {"surface_melt_rate",                 Diagnostic::Ptr(new SurfaceMelt(this, MASS))},
    {"surface_runoff_flux",               Diagnostic::Ptr(new SurfaceRunoff(this, AMOUNT))},
    {"surface_runoff_rate",               Diagnostic::Ptr(new SurfaceRunoff(this, MASS))},
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

TSDiagnosticList SurfaceModel::ts_diagnostics_impl() const {
  using namespace diagnostics;

  TSDiagnosticList result = {
    {"surface_accumulation_rate", TSDiagnostic::Ptr(new TotalSurfaceAccumulation(this))},
    {"surface_melt_rate",         TSDiagnostic::Ptr(new TotalSurfaceMelt(this))},
    {"surface_runoff_rate",       TSDiagnostic::Ptr(new TotalSurfaceRunoff(this))},
  };

  if (m_atmosphere) {
    return pism::combine(result, m_atmosphere->ts_diagnostics());
  }

  if (m_input_model) {
    return pism::combine(result, m_input_model->ts_diagnostics());
  }

  return result;
}

} // end of namespace surface
} // end of namespace pism

