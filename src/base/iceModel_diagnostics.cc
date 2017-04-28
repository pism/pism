// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017 Constantine Khroulev
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

#include "iceModel_diagnostics.hh"

#include "base/rheology/FlowLaw.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/stressbalance/SSB_Modifier.hh"
#include "base/stressbalance/ShallowStressBalance.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec3Custom.hh"
#include "enthalpyConverter.hh"
#include "base/util/PISMVars.hh"
#include "base/util/pism_utilities.hh"
#include "coupler/PISMOcean.hh"
#include "earth/PISMBedDef.hh"

#include "base/grounded_cell_fraction.hh"
#include "base/part_grid_threshold_thickness.hh"
#include "base/util/projection.hh"
#include "base/energy/utilities.hh"
#include "base/energy/EnergyModel.hh"

#if (PISM_USE_PROJ4==1)
#include "base/util/Proj.hh"
#endif

namespace pism {

// Horrendous names used by InitMIP (and ISMIP6, and CMIP5). Ugh.
static const char* land_ice_area_fraction_name           = "sftgif";
static const char* grounded_ice_sheet_area_fraction_name = "sftgrf";
static const char* floating_ice_sheet_area_fraction_name = "sftflf";

namespace diagnostics {

enum AreaType {GROUNDED, FLOATING, BOTH};

enum SurfaceType {TOP, BOTTOM};

/*! @brief Ocean pressure difference at calving fronts. Used to debug CF boundary conditins. */
class CalvingFrontPressureDifference : public Diag<IceModel>
{
public:
  CalvingFrontPressureDifference(IceModel *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};


CalvingFrontPressureDifference::CalvingFrontPressureDifference(IceModel *m)
  : Diag<IceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "ocean_pressure_difference")};
  m_vars[0].set_double("_FillValue", m_fill_value);

  set_attrs("ocean pressure difference at calving fronts", "",
            "", "", 0);
}

IceModelVec::Ptr CalvingFrontPressureDifference::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "ocean_pressure_difference", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  IceModelVec2CellType mask;
  mask.create(m_grid, "mask", WITH_GHOSTS);

  const IceModelVec2S &H   = model->geometry().ice_thickness;
  const IceModelVec2S &bed = model->geometry().bed_elevation;
  const double sea_level = model->ocean_model()->sea_level_elevation(); // FIXME: use 2D sea level

  {
    const double H_threshold = m_config->get_double("stress_balance.ice_free_thickness_standard");
    GeometryCalculator gc(*m_config);
    gc.set_icefree_thickness(H_threshold);

    gc.compute_mask(sea_level, bed, H, mask);
  }

  const double
    rho_ice   = m_config->get_double("constants.ice.density"),
    rho_ocean = m_config->get_double("constants.sea_water.density"),
    g         = m_config->get_double("constants.standard_gravity");

  const bool dry_mode = m_config->get_boolean("ocean.always_grounded");

  IceModelVec::AccessList list{&H, &bed, &mask, result.get()};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.icy(i, j) and mask.next_to_ice_free_ocean(i, j)) {
        (*result)(i, j) = stressbalance::ocean_pressure_difference(mask.ocean(i, j), dry_mode,
                                                                   H(i, j), bed(i, j), sea_level,
                                                                   rho_ice, rho_ocean, g);
      } else {
        (*result)(i, j) = m_fill_value;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return result;
}


/*! @brief Report average basal mass balance flux over the reporting interval (grounded or floating
  areas) */
class BMBSplit : public DiagAverageRate<IceModel>
{
public:
  BMBSplit(const IceModel *m, AreaType flag)
    : DiagAverageRate<IceModel>(m,
                            flag == GROUNDED
                            ? "basal_grounded_mass_flux"
                            : "basal_floating_mass_flux",
                            TOTAL_CHANGE), m_kind(flag) {
    assert(flag != BOTH);

    std::string name, description;
    if (m_kind == GROUNDED) {
      name        = "basal_grounded_mass_flux";
      description = "average basal mass flux over the reporting interval (grounded areas)";
    } else {
      name        = "basal_floating_mass_flux";
      description = "average basal mass flux over the reporting interval (floating areas)";
    }

    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata().set_string("units", "kg m-2");

    set_attrs(description, "", "kg m-2 s-1", "kg m-2 year-1", 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value, "year-1", "second-1");
    m_vars[0].set_double("_FillValue", fill_value);
    m_vars[0].set_string("comment", "positive flux corresponds to ice gain");
  }

protected:
  AreaType m_kind;
  void update_impl(double dt) {
    const IceModelVec2S &input = model->geometry_evolution().bottom_surface_mass_balance();
    const IceModelVec2CellType &cell_type = model->geometry().cell_type;

    IceModelVec::AccessList list{&input, &cell_type, &m_accumulator};

    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        if (m_kind == GROUNDED and cell_type.grounded(i, j)) {
          m_accumulator(i, j) += input(i, j);
        } else if (m_kind == FLOATING and cell_type.ocean(i, j)) {
          m_accumulator(i, j) += input(i, j);
        } else {
          m_accumulator(i, j) = 0.0;
        }
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();

    m_interval_length += dt;
  }
};

HardnessAverage::HardnessAverage(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "hardav")};

  // choice to use SSA power; see #285
  const double power = 1.0 / m_config->get_double("stress_balance.ssa.Glen_exponent");
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);

  set_attrs("vertical average of ice hardness", "",
            unitstr, unitstr, 0);

  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

//! \brief Computes vertically-averaged ice hardness.
IceModelVec::Ptr HardnessAverage::compute_impl() const {

  const rheology::FlowLaw *flow_law = model->stress_balance()->shallow()->flow_law();
  if (flow_law == NULL) {
    flow_law = model->stress_balance()->modifier()->flow_law();
    if (flow_law == NULL) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Can't compute vertically-averaged hardness: no flow law is used.");
    }
  }

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "hardav", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  const IceModelVec3& ice_enthalpy = model->energy_balance_model()->enthalpy();
  const IceModelVec2S& ice_thickness = model->geometry().ice_thickness;

  IceModelVec::AccessList list{&cell_type, &ice_enthalpy, &ice_thickness, result.get()};
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *Eij = ice_enthalpy.get_column(i,j);
      const double H = ice_thickness(i,j);
      if (cell_type.icy(i, j)) {
        (*result)(i,j) = rheology::averaged_hardness(*flow_law,
                                                     H, m_grid->kBelowHeight(H),
                                                     &(m_grid->z()[0]), Eij);
      } else { // put negative value below valid range
        (*result)(i,j) = m_fill_value;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


Rank::Rank(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "rank")};
  set_attrs("processor rank", "", "1", "", 0);
  m_vars[0].set_time_independent(true);
  m_vars[0].set_output_type(PISM_INT);
}

IceModelVec::Ptr Rank::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "rank", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  IceModelVec::AccessList list{result.get()};

  for (Points p(*m_grid); p; p.next()) {
    (*result)(p.i(),p.j()) = m_grid->rank();
  }

  return result;
}


CTS::CTS(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "cts", m_grid->z())};

  set_attrs("cts = E/E_s(p), so cold-temperate transition surface is at cts = 1", "",
            "", "", 0);
}

IceModelVec::Ptr CTS::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "cts", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  energy::compute_cts(model->energy_balance_model()->enthalpy(),
                      model->geometry().ice_thickness, *result);

  return result;
}

Temperature::Temperature(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "temp", m_grid->z())};

  set_attrs("ice temperature", "land_ice_temperature", "K", "K", 0);
  m_vars[0].set_double("valid_min", 0);
}

IceModelVec::Ptr Temperature::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3(m_grid, "temp", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  const IceModelVec2S &thickness = model->geometry().ice_thickness;
  const IceModelVec3 &enthalpy = model->energy_balance_model()->enthalpy();

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  double *Tij;
  const double *Enthij; // columns of these values

  IceModelVec::AccessList list{result.get(), &enthalpy, &thickness};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      Tij = result->get_column(i,j);
      Enthij = enthalpy.get_column(i,j);
      for (unsigned int k=0; k <m_grid->Mz(); ++k) {
        const double depth = thickness(i,j) - m_grid->z(k);
        Tij[k] = EC->temperature(Enthij[k], EC->pressure(depth));
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


TemperaturePA::TemperaturePA(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "temp_pa", m_grid->z())};

  set_attrs("pressure-adjusted ice temperature (degrees above pressure-melting point)", "",
            "deg_C", "deg_C", 0);
  m_vars[0].set_double("valid_max", 0);
}

IceModelVec::Ptr TemperaturePA::compute_impl() const {
  bool cold_mode = m_config->get_boolean("energy.temperature_based");
  double melting_point_temp = m_config->get_double("constants.fresh_water.melting_point_temperature");

  IceModelVec3::Ptr result(new IceModelVec3(m_grid, "temp_pa", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  const IceModelVec2S &thickness = model->geometry().ice_thickness;
  const IceModelVec3  &enthalpy  = model->energy_balance_model()->enthalpy();

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  double *Tij;
  const double *Enthij; // columns of these values

  IceModelVec::AccessList list{result.get(), &enthalpy, &thickness};

  ParallelSection loop(m_grid->com);
  try {
    for (Points pt(*m_grid); pt; pt.next()) {
      const int i = pt.i(), j = pt.j();

      Tij = result->get_column(i,j);
      Enthij = enthalpy.get_column(i,j);
      for (unsigned int k=0; k < m_grid->Mz(); ++k) {
        const double depth = thickness(i,j) - m_grid->z(k),
          p = EC->pressure(depth);
        Tij[k] = EC->pressure_adjusted_temperature(Enthij[k], p);

        if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
          // is 273.15
          if (EC->is_temperate_relaxed(Enthij[k],p) && (thickness(i,j) > 0)) {
            Tij[k] = melting_point_temp;
          }
        }

      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->shift(-melting_point_temp);

  return result;
}

TemperaturePABasal::TemperaturePABasal(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "temppabase")};

  set_attrs("pressure-adjusted ice temperature at the base of ice", "",
            "Celsius", "Celsius", 0);
}

IceModelVec::Ptr TemperaturePABasal::compute_impl() const {

  bool cold_mode = m_config->get_boolean("energy.temperature_based");
  double melting_point_temp = m_config->get_double("constants.fresh_water.melting_point_temperature");

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "temp_pa_base", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  const IceModelVec2S &thickness = model->geometry().ice_thickness;
  const IceModelVec3 &enthalpy = model->energy_balance_model()->enthalpy();

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  const double *Enthij; // columns of these values

  IceModelVec::AccessList list{result.get(), &enthalpy, &thickness};

  ParallelSection loop(m_grid->com);
  try {
    for (Points pt(*m_grid); pt; pt.next()) {
      const int i = pt.i(), j = pt.j();

      Enthij = enthalpy.get_column(i,j);

      const double depth = thickness(i,j),
        p = EC->pressure(depth);
      (*result)(i,j) = EC->pressure_adjusted_temperature(Enthij[0], p);

      if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
        // is 273.15
        if (EC->is_temperate_relaxed(Enthij[0],p) && (thickness(i,j) > 0)) {
          (*result)(i,j) = melting_point_temp;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->shift(-melting_point_temp);

  return result;
}

EnthalpySurface::EnthalpySurface(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "enthalpysurf")};

  set_attrs("ice enthalpy at 1m below the ice surface", "",
            "J kg-1", "J kg-1", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr EnthalpySurface::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "enthalpysurf", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  // compute levels corresponding to 1 m below the ice surface:

  const IceModelVec3& ice_enthalpy = model->energy_balance_model()->enthalpy();
  const IceModelVec2S& ice_thickness = model->geometry().ice_thickness;

  IceModelVec::AccessList list{&ice_thickness, result.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) = std::max(ice_thickness(i,j) - 1.0, 0.0);
  }

  ice_enthalpy.getSurfaceValues(*result, *result);  // z=0 slice

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ice_thickness(i,j) <= 1.0) {
      (*result)(i,j) = m_fill_value;
    }
  }

  return result;
}

EnthalpyBasal::EnthalpyBasal(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "enthalpybase")};

  set_attrs("ice enthalpy at the base of ice", "",
            "J kg-1", "J kg-1", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr EnthalpyBasal::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "enthalpybase", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  model->energy_balance_model()->enthalpy().getHorSlice(*result, 0.0);  // z=0 slice

  result->mask_by(model->geometry().ice_thickness, m_fill_value);

  return result;
}


TemperatureBasal::TemperatureBasal(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "tempbase")};

  set_attrs("ice temperature at the base of ice",
            "land_ice_basal_temperature", // InitMIP "standard" name
            "K", "K", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr TemperatureBasal::compute_impl() const {

  const IceModelVec2S &thickness = model->geometry().ice_thickness;

  IceModelVec::Ptr enth = EnthalpyBasal(model).compute();

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(enth);

  // result contains basal enthalpy; note that it is allocated by
  // EnthalpyBasal::compute().

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  IceModelVec::AccessList list{&cell_type, result.get(), &thickness};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double depth = thickness(i,j),
        pressure = EC->pressure(depth);
      if (cell_type.icy(i, j)) {
        (*result)(i,j) = EC->temperature((*result)(i,j), pressure);
      } else {
        (*result)(i,j) = m_fill_value;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->metadata(0) = m_vars[0];
  return result;
}

TemperatureSurface::TemperatureSurface(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "tempsurf")};

  set_attrs("ice temperature at 1m below the ice surface",
            "temperature_at_ground_level_in_snow_or_firn", // InitMIP "standard" name
            "K", "K", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr TemperatureSurface::compute_impl() const {

  const IceModelVec2S &thickness = model->geometry().ice_thickness;

  IceModelVec::Ptr enth = EnthalpySurface(model).compute();
  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(enth);

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  // result contains surface enthalpy; note that it is allocated by
  // EnthalpySurface::compute().

  IceModelVec::AccessList list{result.get(), &thickness};

  double depth = 1.0,
    pressure = EC->pressure(depth);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (thickness(i,j) > 1) {
        (*result)(i,j) = EC->temperature((*result)(i,j), pressure);
      } else {
        (*result)(i,j) = m_fill_value;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->metadata(0) = m_vars[0];
  return result;
}


LiquidFraction::LiquidFraction(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "liqfrac", m_grid->z())};

  set_attrs("liquid water fraction in ice (between 0 and 1)", "",
            "1", "1", 0);
  m_vars[0].set_doubles("valid_range", {0.0, 1.0});
}

IceModelVec::Ptr LiquidFraction::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "liqfrac", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  bool cold_mode = m_config->get_boolean("energy.temperature_based");

  if (cold_mode) {
    result->set(0.0);
  } else {
    energy::compute_liquid_water_fraction(model->energy_balance_model()->enthalpy(),
                                          model->geometry().ice_thickness,
                                          *result);
  }

  return result;
}

TemperateIceThickness::TemperateIceThickness(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys,
                                    "tempicethk")};

  set_attrs("temperate ice thickness (total column content)", "",
            "m", "m", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr TemperateIceThickness::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "tempicethk", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;
  const IceModelVec3& ice_enthalpy = model->energy_balance_model()->enthalpy();
  const IceModelVec2S& ice_thickness = model->geometry().ice_thickness;

  IceModelVec::AccessList list{&cell_type, result.get(), &ice_enthalpy, &ice_thickness};

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.icy(i, j)) {
        const double *Enth = ice_enthalpy.get_column(i,j);
        double H_temperate = 0.0;
        const double H = ice_thickness(i,j);
        const unsigned int ks = m_grid->kBelowHeight(H);

        for (unsigned int k=0; k<ks; ++k) { // FIXME issue #15
          double pressure = EC->pressure(H - m_grid->z(k));

          if (EC->is_temperate_relaxed(Enth[k], pressure)) {
            H_temperate += m_grid->z(k+1) - m_grid->z(k);
          }
        }

        double pressure = EC->pressure(H - m_grid->z(ks));
        if (EC->is_temperate_relaxed(Enth[ks], pressure)) {
          H_temperate += H - m_grid->z(ks);
        }

        (*result)(i,j) = H_temperate;
      } else {
        // ice-free
        (*result)(i,j) = m_fill_value;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

TemperateIceThicknessBasal::TemperateIceThicknessBasal(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys,
                                    "tempicethk_basal")};

  set_attrs("thickness of the basal layer of temperate ice", "",
            "m", "m", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

/*!
 * Uses linear interpolation to go beyond vertical grid resolution.
 */
IceModelVec::Ptr TemperateIceThicknessBasal::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "tempicethk_basal", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;
  const IceModelVec3& ice_enthalpy = model->energy_balance_model()->enthalpy();
  const IceModelVec2S& ice_thickness = model->geometry().ice_thickness;

  IceModelVec::AccessList list{&cell_type, result.get(), &ice_thickness, &ice_enthalpy};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double H = ice_thickness(i,j);

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (cell_type.ice_free(i, j)) {
        (*result)(i,j) = m_fill_value;
        continue;
      }

      const double *Enth = ice_enthalpy.get_column(i,j);

      unsigned int ks = m_grid->kBelowHeight(H);

      unsigned int k = 0;
      double pressure = EC->pressure(H - m_grid->z(k));
      while (k <= ks) {         // FIXME issue #15
        pressure = EC->pressure(H - m_grid->z(k));

        if (EC->is_temperate_relaxed(Enth[k],pressure)) {
          k++;
        } else {
          break;
        }
      }
      // after this loop 'pressure' is equal to the pressure at the first level
      // that is cold

      // no temperate ice at all; go to the next grid point
      if (k == 0) {
        (*result)(i,j) = 0.0;
        continue;
      }

      // the whole column is temperate (except, possibly, some ice between
      // z(ks) and the total thickness; we ignore it)
      if (k == ks + 1) {
        (*result)(i,j) = m_grid->z(ks);
        continue;
      }

      double
        pressure_0 = EC->pressure(H - m_grid->z(k-1)),
        dz         = m_grid->z(k) - m_grid->z(k-1),
        slope1     = (Enth[k] - Enth[k-1]) / dz,
        slope2     = (EC->enthalpy_cts(pressure) - EC->enthalpy_cts(pressure_0)) / dz;

      if (slope1 != slope2) {
        (*result)(i,j) = m_grid->z(k-1) +
          (EC->enthalpy_cts(pressure_0) - Enth[k-1]) / (slope1 - slope2);

        // check if the resulting thickness is valid:
        (*result)(i,j) = std::max((*result)(i,j), m_grid->z(k-1));
        (*result)(i,j) = std::min((*result)(i,j), m_grid->z(k));
      } else {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Linear interpolation of the thickness of"
                                      " the basal temperate layer failed:\n"
                                      "(i=%d, j=%d, k=%d, ks=%d)\n",
                                      i, j, k, ks);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return result;
}

VolumeGlacierized::VolumeGlacierized(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "volume_glacierized") {

  m_ts.variable().set_string("units", "m3");
  m_ts.variable().set_string("long_name", "volume of the ice in glacierized areas");
  m_ts.variable().set_double("valid_min", 0.0);
}

double VolumeGlacierized::compute() {
  return model->ice_volume(m_config->get_double("output.ice_free_thickness_standard"));
}

VolumeNonGlacierized::VolumeNonGlacierized(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "volume_nonglacierized") {

  m_ts.variable().set_string("units", "m3");
  m_ts.variable().set_string("long_name", "volume of the ice, including seasonal cover");
  m_ts.variable().set_double("valid_min", 0.0);
}

double VolumeNonGlacierized::compute() {
  return model->ice_volume(0.0);
}

SeaLevelVolume::SeaLevelVolume(const IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "slvol") {

  m_ts.variable().set_string("units", "m");
  m_ts.variable().set_string("long_name", "total sea-level relevant ice IN SEA-LEVEL EQUIVALENT");
  m_ts.variable().set_double("valid_min", 0.0);
}

double SeaLevelVolume::compute() {
  return model->sealevel_volume(m_config->get_double("output.ice_free_thickness_standard"));
}

VolumeRateOfChangeGlacierized::VolumeRateOfChangeGlacierized(IceModel *m)
  : TSDiag<TSRateDiagnostic, IceModel>(m, "volume_rate_of_change_glacierized") {

  m_ts.variable().set_string("units", "m3 s-1");
  m_ts.variable().set_string("glaciological_units", "m3 year-1");
  m_ts.variable().set_string("long_name", "rate of change of the ice volume in glacierized areas");
}

double VolumeRateOfChangeGlacierized::compute() {
  return model->ice_volume(m_config->get_double("output.ice_free_thickness_standard"));
}

VolumeRateOfChangeNonGlacierized::VolumeRateOfChangeNonGlacierized(IceModel *m)
  : TSDiag<TSRateDiagnostic, IceModel>(m, "volume_rate_of_change_nonglacierized") {

  m_ts.variable().set_string("units", "m3 s-1");
  m_ts.variable().set_string("glaciological_units", "m3 year-1");
  m_ts.variable().set_string("long_name",
                             "rate of change of the ice volume, including seasonal cover");
}

double VolumeRateOfChangeNonGlacierized::compute() {
  return model->ice_volume(0.0);
}


AreaGlacierized::AreaGlacierized(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "area_glacierized") {

  m_ts.variable().set_string("units", "m2");
  m_ts.variable().set_string("long_name", "glacierized area");
  m_ts.variable().set_double("valid_min", 0.0);
}

double AreaGlacierized::compute() {
  return model->ice_area(m_config->get_double("output.ice_free_thickness_standard"));
}

MassGlacierized::MassGlacierized(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "mass_glacierized") {

  m_ts.variable().set_string("units", "kg");
  m_ts.variable().set_string("long_name", "mass of the ice in glacierized areas");
  m_ts.variable().set_double("valid_min", 0.0);
}

double MassGlacierized::compute() {
  return (model->ice_volume(m_config->get_double("output.ice_free_thickness_standard")) *
          m_grid->ctx()->config()->get_double("constants.ice.density"));
}

MassNonGlacierized::MassNonGlacierized(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "mass_nonglacierized") {

  m_ts.variable().set_string("units", "kg");
  m_ts.variable().set_string("long_name", "mass of the ice, including seasonal cover");
  m_ts.variable().set_double("valid_min", 0.0);
}

double MassNonGlacierized::compute() {
  return (model->ice_volume(0.0) *
          m_grid->ctx()->config()->get_double("constants.ice.density"));
}

MassRateOfChangeGlacierized::MassRateOfChangeGlacierized(IceModel *m)
  : TSDiag<TSRateDiagnostic, IceModel>(m, "mass_rate_of_change_glacierized") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("glaciological_units", "kg year-1");
  m_ts.variable().set_string("long_name", "rate of change of the mass of ice in glacierized areas");
}

double MassRateOfChangeGlacierized::compute() {

  const double
    ice_density = m_grid->ctx()->config()->get_double("constants.ice.density"),
    ice_volume  = model->ice_volume(m_config->get_double("output.ice_free_thickness_standard"));

  return ice_volume * ice_density;
}

MassRateOfChangeNonGlacierized::MassRateOfChangeNonGlacierized(IceModel *m)
  : TSDiag<TSRateDiagnostic, IceModel>(m, "mass_rate_of_change_nonglacierized") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("glaciological_units", "kg year-1");
  m_ts.variable().set_string("long_name",
                             "rate of change of the mass of ice, including seasonal cover");
}

double MassRateOfChangeNonGlacierized::compute() {
  const double ice_density = m_grid->ctx()->config()->get_double("constants.ice.density");
  return model->ice_volume(0.0) * ice_density;
}


VolumeGlacierizedTemperate::VolumeGlacierizedTemperate(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "volume_glacierized_temperate") {

  m_ts.variable().set_string("units", "m3");
  m_ts.variable().set_string("long_name", "volume of temperate ice in glacierized areas");
  m_ts.variable().set_double("valid_min", 0.0);
}

double VolumeGlacierizedTemperate::compute() {
  return model->ice_volume_temperate(m_config->get_double("output.ice_free_thickness_standard"));
}

VolumeNonGlacierizedTemperate::VolumeNonGlacierizedTemperate(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "volume_nonglacierized_temperate") {

  m_ts.variable().set_string("units", "m3");
  m_ts.variable().set_string("long_name", "volume of temperate ice, including seasonal cover");
  m_ts.variable().set_double("valid_min", 0.0);
}

double VolumeNonGlacierizedTemperate::compute() {
  return model->ice_volume_temperate(0.0);
}


VolumeGlacierizedCold::VolumeGlacierizedCold(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "volume_glacierized_cold") {

  m_ts.variable().set_string("units", "m3");
  m_ts.variable().set_string("long_name", "volume of cold ice in glacierized areas");
  m_ts.variable().set_double("valid_min", 0.0);
}

double VolumeGlacierizedCold::compute() {
  return model->ice_volume_cold(m_config->get_double("output.ice_free_thickness_standard"));
}

VolumeNonGlacierizedCold::VolumeNonGlacierizedCold(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "volume_nonglacierized_cold") {

  m_ts.variable().set_string("units", "m3");
  m_ts.variable().set_string("long_name", "volume of cold ice, including seasonal cover");
  m_ts.variable().set_double("valid_min", 0.0);
}

double VolumeNonGlacierizedCold::compute() {
  return model->ice_volume_cold(0.0);
}

AreaGlacierizedTemperateBase::AreaGlacierizedTemperateBase(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "area_glacierized_temperate_base") {

  m_ts.variable().set_string("units", "m2");
  m_ts.variable().set_string("long_name", "glacierized area where basal ice is temperate");
  m_ts.variable().set_double("valid_min", 0.0);
}

double AreaGlacierizedTemperateBase::compute() {
  return model->ice_area_temperate(m_config->get_double("output.ice_free_thickness_standard"));
}

AreaGlacierizedColdBase::AreaGlacierizedColdBase(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "area_glacierized_cold_base") {

  m_ts.variable().set_string("units", "m2");
  m_ts.variable().set_string("long_name", "glacierized area where basal ice is cold");
  m_ts.variable().set_double("valid_min", 0.0);
}

double AreaGlacierizedColdBase::compute() {
  return model->ice_area_cold(m_config->get_double("output.ice_free_thickness_standard"));
}

EnthalpyGlacierized::EnthalpyGlacierized(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "enthalpy_glacierized") {

  m_ts.variable().set_string("units", "J");
  m_ts.variable().set_string("long_name", "enthalpy of the ice in glacierized areas");
  m_ts.variable().set_double("valid_min", 0.0);
}

double EnthalpyGlacierized::compute() {
  return energy::total_ice_enthalpy(m_config->get_double("output.ice_free_thickness_standard"),
                                    model->energy_balance_model()->enthalpy(),
                                    model->geometry().ice_thickness,
                                    model->geometry().cell_area);
}

EnthalpyNonGlacierized::EnthalpyNonGlacierized(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "enthalpy_nonglacierized") {

  m_ts.variable().set_string("units", "J");
  m_ts.variable().set_string("long_name", "enthalpy of the ice, including seasonal cover");
  m_ts.variable().set_double("valid_min", 0.0);
}

double EnthalpyNonGlacierized::compute() {
  return energy::total_ice_enthalpy(0.0,
                                    model->energy_balance_model()->enthalpy(),
                                    model->geometry().ice_thickness,
                                    model->geometry().cell_area);
}

AreaGlacierizedGrounded::AreaGlacierizedGrounded(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "area_glacierized_grounded") {

  m_ts.variable().set_string("units", "m2");
  m_ts.variable().set_string("long_name", "area of grounded ice in glacierized areas");
  m_ts.variable().set_double("valid_min", 0.0);
}

double AreaGlacierizedGrounded::compute() {
  return model->ice_area_grounded(m_config->get_double("output.ice_free_thickness_standard"));
}

AreaGlacierizedShelf::AreaGlacierizedShelf(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "area_glacierized_shelf") {

  m_ts.variable().set_string("units", "m2");
  m_ts.variable().set_string("long_name", "area of ice shelves in glacierized areas");
  m_ts.variable().set_double("valid_min", 0.0);
}

double AreaGlacierizedShelf::compute() {
  return model->ice_area_floating(m_config->get_double("output.ice_free_thickness_standard"));
}

TimeStepLength::TimeStepLength(const IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "dt") {

  m_ts.variable().set_string("units", "second");
  m_ts.variable().set_string("glaciological_units", "year");
  m_ts.variable().set_string("long_name", "mass continuity time step");
  m_ts.variable().set_double("valid_min", 0.0);
}

double TimeStepLength::compute() {
  return model->dt();
}

MaxDiffusivity::MaxDiffusivity(const IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "max_diffusivity") {

  m_ts.variable().set_string("units", "m2 s-1");
  m_ts.variable().set_string("long_name", "maximum diffusivity");
  m_ts.variable().set_double("valid_min", 0.0);
}

double MaxDiffusivity::compute() {
  return model->stress_balance()->max_diffusivity();
}

double mass_change(const IceModel *model, SurfaceType surface, AreaType area) {
  const IceGrid &grid = *model->grid();
  const Config &config = *grid.ctx()->config();

  const double ice_density = config.get_double("constants.ice.density");

  const IceModelVec2S &cell_area = model->geometry().cell_area;
  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  const IceModelVec2S &thickness_change = (surface == TOP) ?
    model->geometry_evolution().top_surface_mass_balance() :
    model->geometry_evolution().bottom_surface_mass_balance();

  IceModelVec::AccessList list{&cell_area, &cell_type, &thickness_change};

  double volume_change = 0.0;
  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if ((area == BOTH) or
          (area == GROUNDED and cell_type.grounded(i, j)) or
          (area == FLOATING and cell_type.ocean(i, j))) {

        volume_change += cell_area(i, j) * thickness_change(i, j);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // (kg / m3) * m3 = kg
  return ice_density * GlobalSum(grid.com, volume_change);
}

MassFluxSurface::MassFluxSurface(const IceModel *m)
  : TSDiag<TSFluxDiagnostic, IceModel>(m, "surface_ice_flux") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("glaciological_units", "kg year-1");
  m_ts.variable().set_string("long_name", "total over ice domain of top surface ice mass flux");
  m_ts.variable().set_string("comment", "positive means ice gain");
}

double MassFluxSurface::compute() {
  return mass_change(model, TOP, BOTH);
}

MassFluxBasalGrounded::MassFluxBasalGrounded(const IceModel *m)
  : TSDiag<TSFluxDiagnostic, IceModel>(m, "grounded_basal_ice_flux") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("glaciological_units", "kg year-1");
  m_ts.variable().set_string("long_name", "total over grounded ice domain of basal mass flux");
  m_ts.variable().set_string("comment", "positive means ice gain");
}

double MassFluxBasalGrounded::compute() {
  return mass_change(model, BOTTOM, BOTH);
}

MassFluxBasalFloating::MassFluxBasalFloating(const IceModel *m)
  : TSDiag<TSFluxDiagnostic, IceModel>(m, "sub_shelf_ice_flux") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("glaciological_units", "kg year-1");
  m_ts.variable().set_string("long_name", "total sub-shelf ice flux");
  m_ts.variable().set_string("comment", "positive means ice gain");
}

double MassFluxBasalFloating::compute() {
  return mass_change(model, BOTTOM, FLOATING);
}

MassFluxDischarge::MassFluxDischarge(const IceModel *m)
  : TSDiag<TSFluxDiagnostic, IceModel>(m, "discharge_flux") {

  m_ts.variable().set_string("units", "kg s-1");
  m_ts.variable().set_string("glaciological_units", "kg year-1");
  m_ts.variable().set_string("long_name", "discharge (calving & icebergs) flux");
  m_ts.variable().set_string("comment", "positive means ice gain");
}

double MassFluxDischarge::compute() {
  return 0.0;                   // FIXME_ (not zero in general)
}

//! \brief Computes dHdt, the ice thickness rate of change.
class ThicknessRateOfChange : public Diag<IceModel>
{
public:
  ThicknessRateOfChange(const IceModel *m)
    : Diag<IceModel>(m),
    m_last_thickness(m_grid, "last_ice_thickness", WITHOUT_GHOSTS),
    m_interval_length(0.0) {

    // set metadata:
    m_vars = {SpatialVariableMetadata(m_sys, "dHdt")};

    set_attrs("ice thickness rate of change",
              "tendency_of_land_ice_thickness",
              "m second-1", "m year-1", 0);

    m_fill_value = units::convert(m_sys, m_fill_value,
                                  "m year-1", "m second-1");

    const double valid_range = units::convert(m_sys, -1e6, "m year-1", "m second-1");

    m_vars[0].set_doubles("valid_range",  {-valid_range, valid_range});
    m_vars[0].set_double("_FillValue", m_fill_value);
    m_vars[0].set_string("cell_methods", "time: mean");

    m_last_thickness.set_attrs("internal",
                               "ice thickness at the time of the last report of dHdt",
                               "m", "land_ice_thickness");
  }
protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "dHdt", WITHOUT_GHOSTS));
    result->metadata() = m_vars[0];

    if (m_interval_length > 0.0) {
      model->geometry().ice_thickness.add(-1.0, m_last_thickness, *result);
      result->scale(1.0 / m_interval_length);
    } else {
      result->set(m_fill_value);
    }

    return result;
  }

  void reset_impl() {
    m_interval_length = 0.0;
    m_last_thickness.copy_from(model->geometry().ice_thickness);
  }

  void update_impl(double dt) {
    m_interval_length += dt;
  }

protected:
  IceModelVec2S m_last_thickness;
  double m_interval_length;
};

VolumeGlacierizedGrounded::VolumeGlacierizedGrounded(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "volume_glacierized_grounded") {

  m_ts.variable().set_string("units", "m3");
  m_ts.variable().set_string("long_name", "volume of grounded ice in glacierized areas");
  m_ts.variable().set_double("valid_min", 0.0);
}

double VolumeGlacierizedGrounded::compute() {
  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  const IceModelVec2S
    &cell_area     = model->geometry().cell_area,
    &ice_thickness = model->geometry().ice_thickness;

  const double thickness_threshold = m_config->get_double("output.ice_free_thickness_standard");

  IceModelVec::AccessList list{&ice_thickness, &cell_type, &cell_area};

  double volume = 0.0;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double H = ice_thickness(i, j);

    if (cell_type.grounded(i, j) and H >= thickness_threshold) {
      volume += cell_area(i, j) * H;
    }
  }

  return GlobalSum(m_grid->com, volume);
}

VolumeGlacierizedShelf::VolumeGlacierizedShelf(IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "volume_glacierized_shelf") {

  m_ts.variable().set_string("units", "m3");
  m_ts.variable().set_string("long_name", "volume of ice shelves in glacierized areas");
  m_ts.variable().set_double("valid_min", 0.0);
}

double VolumeGlacierizedShelf::compute() {
  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  const IceModelVec2S
    &cell_area     = model->geometry().cell_area,
    &ice_thickness = model->geometry().ice_thickness;

  const double thickness_threshold = m_config->get_double("output.ice_free_thickness_standard");

  IceModelVec::AccessList list{&ice_thickness, &cell_type, &cell_area};

  double volume = 0.0;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double H = ice_thickness(i, j);

    if (cell_type.ocean(i, j) and H >= thickness_threshold) {
      volume += cell_area(i, j) * H;
    }
  }

  return GlobalSum(m_grid->com, volume);
}

//! \brief Reports the maximum horizontal absolute velocity component over the grid.
/*!
 * This is the value used by the adaptive time-stepping code in the CFL condition
 * for horizontal advection (i.e. in energy and mass conservation time steps).
 *
 * This is not the maximum horizontal speed, but rather the maximum of components.
 *
 * Note that this picks up the value computed during the time-step taken at a
 * reporting time. (It is not the "average over the reporting interval computed using
 * differencing in time", as other rate-of-change diagnostics.)
 */
MaxHorizontalVelocity::MaxHorizontalVelocity(const IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "max_hor_vel") {

  m_ts.variable().set_string("units", "m second-1");
  m_ts.variable().set_string("glaciological_units", "m year-1");
  m_ts.variable().set_string("long_name",
                             "maximum abs component of horizontal ice velocity"
                             " over grid in last time step during time-series reporting interval");
  m_ts.variable().set_double("valid_min", 0.0);
}

double MaxHorizontalVelocity::compute() {
  CFLData cfl = model->stress_balance()->max_timestep_cfl_3d();

  return std::max(cfl.u_max, cfl.v_max);
}

MassNotDisplacingSeaWater::MassNotDisplacingSeaWater(const IceModel *m)
  : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "limnsw") {

  m_ts.variable().set_string("units", "kg");
  m_ts.variable().set_string("long_name", "mass of the ice not displacing sea water");
  m_ts.variable().set_double("valid_min", 0.0);
}

double MassNotDisplacingSeaWater::compute() {

  const double
    ice_density = m_config->get_double("constants.ice.density"),
    ice_volume  = model->ice_volume_not_displacing_seawater(m_config->get_double("output.ice_free_thickness_standard")),
    ice_mass    = ice_volume * ice_density;

  return ice_mass;
}

LatLonBounds::LatLonBounds(const IceModel *m,
                           const std::string &var_name,
                           const std::string &proj_string)
  : Diag<IceModel>(m) {
  assert(var_name == "lat" || var_name == "lon");
  m_var_name = var_name;

  // set metadata:
  std::vector<double> levels(4);
  for (int k = 0; k < 4; ++k) {
    levels[k] = k;
  }

  m_vars = {SpatialVariableMetadata(m_sys, m_var_name + "_bnds", levels)};
  m_vars[0].get_z().set_name("nv4");
  m_vars[0].get_z().clear_all_strings();
  m_vars[0].get_z().clear_all_doubles();
  m_vars[0].set_time_independent(true);

  if (m_var_name == "lon") {
    set_attrs("longitude bounds", "", "degree_east", "degree_east", 0);
    m_vars[0].set_double("valid_min", -180);
    m_vars[0].set_double("valid_max", 180);
  } else {
    set_attrs("latitude bounds", "", "degree_north", "degree_north", 0);
    m_vars[0].set_double("valid_min", -90);
    m_vars[0].set_double("valid_max", 90);
  }
  m_vars[0].set_string("coordinates", "");

  m_proj_string = proj_string;

#if (PISM_USE_PROJ4==1)
  // create PROJ.4 objects to check if proj_string is OK.
  Proj lonlat("+proj=latlong +datum=WGS84 +ellps=WGS84");
  Proj pism(m_proj_string);
#endif
  // If PISM_USE_PROJ4 is not 1 we don't need to check validity of m_proj_string: this diagnostic
  // will not be available and so this code will not run.
}

IceModelVec::Ptr LatLonBounds::compute_impl() const {
  std::map<std::string,std::string> attrs;
  std::vector<double> indices(4);

  IceModelVec3Custom::Ptr result(new IceModelVec3Custom);
  result->create(m_grid, m_var_name + "_bnds", "nv4",
                 indices, attrs);
  result->metadata(0) = m_vars[0];

  bool latitude = true;
  if (m_var_name == "lon") {
    latitude = false;
  }

  if (latitude) {
    compute_lat_bounds(m_proj_string, *result);
  } else {
    compute_lon_bounds(m_proj_string, *result);
  }

  return result;
}

IceAreaFraction::IceAreaFraction(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, land_ice_area_fraction_name)};
  set_attrs("fraction of a grid cell covered by ice (grounded or floating)",
            "land_ice_area_fraction", // InitMIP "standard" name
            "1", "1", 0);
}

IceModelVec::Ptr IceAreaFraction::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, land_ice_area_fraction_name, WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  const IceModelVec2S
    &thickness         = model->geometry().ice_thickness,
    &surface_elevation = model->geometry().ice_surface_elevation,
    &bed_topography    = model->geometry().bed_elevation;

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  IceModelVec::AccessList list{&thickness, &surface_elevation, &bed_topography, &cell_type,
      result.get()};

  const bool do_part_grid = m_config->get_boolean("geometry.part_grid.enabled");
  const IceModelVec2S &Href = model->geometry().ice_area_specific_volume;;
  if (do_part_grid) {
    list.add(Href);
  }

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.icy(i, j)) {
        // an "icy" cell: the area fraction is one
        (*result)(i, j) = 1.0;
      } else if (cell_type.ice_free_ocean(i, j)) {
        // an ice-free ocean cell may be "partially-filled", in which case we need to compute its
        // ice area fraction by dividing Href by the threshold thickness.

        double H_reference = do_part_grid ? Href(i, j) : 0.0;

        if (H_reference > 0.0) {
          const double H_threshold = part_grid_threshold_thickness(cell_type.int_star(i, j),
                                                                   thickness.star(i, j),
                                                                   surface_elevation.star(i, j),
                                                                   bed_topography(i,j));
          // protect from a division by zero
          if (H_threshold > 0.0) {
            (*result)(i, j) = H_reference / H_threshold;
          } else {
            (*result)(i, j) = 1.0;
          }
        } else {
          // H_reference is zero
          (*result)(i, j) = 0.0;
        }
      } else {
        // an ice-free-ground cell: the area fraction is zero
        (*result)(i, j) = 0.0;
      }
    } // end of the loop over grid points
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

IceAreaFractionGrounded::IceAreaFractionGrounded(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, grounded_ice_sheet_area_fraction_name)};
  set_attrs("fraction of a grid cell covered by grounded ice",
            "grounded_ice_sheet_area_fraction", // InitMIP "standard" name
            "1", "1", 0);
}

IceModelVec::Ptr IceAreaFractionGrounded::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, grounded_ice_sheet_area_fraction_name, WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  const double sea_level = model->ocean_model()->sea_level_elevation();

  const double
    ice_density   = m_config->get_double("constants.ice.density"),
    ocean_density = m_config->get_double("constants.sea_water.density");

  const IceModelVec2S
    &ice_thickness  = model->geometry().ice_thickness,
    &bed_topography = model->geometry().bed_elevation;

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  compute_grounded_cell_fraction(ice_density, ocean_density, sea_level,
                                 ice_thickness, bed_topography, cell_type,
                                 *result, NULL, NULL);

  // All grounded areas have the grounded cell fraction of one, so now we make sure that ice-free
  // areas get the value of 0 (they are grounded but not covered by a grounded ice sheet).

  IceModelVec::AccessList list{&cell_type, result.get()};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      if (cell_type.ice_free(i, j)) {
        (*result)(i, j) = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

IceAreaFractionFloating::IceAreaFractionFloating(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, floating_ice_sheet_area_fraction_name)};
  set_attrs("fraction of a grid cell covered by floating ice",
            "floating_ice_sheet_area_fraction", // InitMIP "standard" name
            "1", "1", 0);
}

IceModelVec::Ptr IceAreaFractionFloating::compute_impl() const {

  IceAreaFraction land_ice_area_fraction(model);
  IceModelVec::Ptr ice_area_fraction = land_ice_area_fraction.compute();

  IceAreaFractionGrounded grounded_ice_sheet_area_fraction(model);
  IceModelVec::Ptr grounded_area_fraction = grounded_ice_sheet_area_fraction.compute();

  IceModelVec::Ptr result = ice_area_fraction;
  result->metadata() = m_vars[0];

  // Floating area fraction is total area fraction minus grounded area fraction.
  result->add(-1.0, *grounded_area_fraction);

  return result;
}

HeightAboveFloatation::HeightAboveFloatation(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "height_above_flotation")};

  set_attrs("the height above flotation", "",
            "m", "m", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr HeightAboveFloatation::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "height_above_flotation", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  const double
    ice_density   = m_config->get_double("constants.ice.density"),
    ocean_density = m_config->get_double("constants.sea_water.density"),
    sea_level = model->ocean_model()->sea_level_elevation();

  const IceModelVec2S
    &ice_thickness  = model->geometry().ice_thickness,
    &bed_topography = model->geometry().bed_elevation;

  IceModelVec::AccessList list{&cell_type, result.get(), &ice_thickness, &bed_topography};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double thk = ice_thickness(i,j);
      double bed = bed_topography(i,j);

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (cell_type.ice_free(i, j)) {
        (*result)(i,j) = m_fill_value;
        continue;
      }
      (*result)(i,j) = ((ice_density / ocean_density) * thk) + (bed - sea_level);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return result;
}

IceMass::IceMass(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "ice_mass")};

  set_attrs("mass per cell",
            "",                 // no standard name
            "kg", "kg", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceMass::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "ice_mass", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  const double
    ice_density = m_config->get_double("constants.ice.density");

  const IceModelVec2S
    &ice_thickness = model->geometry().ice_thickness,
    &cell_area     = model->geometry().cell_area;

  IceModelVec::AccessList list{&cell_type, result.get(), &ice_thickness, &cell_area};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // count all ice, including cells which have so little they
      // are considered "ice-free"
      if (ice_thickness(i, j) > 0.0) {
        (*result)(i,j) = ice_density * ice_thickness(i, j) * cell_area(i, j);
      } else {
        (*result)(i,j) = m_fill_value;
      }
    } // end of loop over grid points

  } catch (...) {
    loop.failed();
  }
  loop.check();

  // Add the mass of ice in Href:
  if (m_config->get_boolean("geometry.part_grid.enabled")) {
    const IceModelVec2S &Href = model->geometry().ice_area_specific_volume;
    list.add(Href);
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (ice_thickness(i, j) <= 0.0 and Href(i, j) > 0.0) {
        (*result)(i, j) = ice_density * Href(i, j) * cell_area(i,j);
      }
    }
  }

  return result;
}

BedTopographySeaLevelAdjusted::BedTopographySeaLevelAdjusted(const IceModel *m)
  : Diag<IceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "topg_sl_adjusted")};

  set_attrs("sea-level adjusted bed topography (zero is at sea level)", "",
            "meters", "meters", 0);
}

IceModelVec::Ptr BedTopographySeaLevelAdjusted::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "topg_sl_adjusted", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->bed_model()->bed_elevation());
  // result = topg - sea_level
  result->shift(-model->ocean_model()->sea_level_elevation());

  return result;
}

Hardness::Hardness(const IceModel *m)
  : Diag<IceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "hardness", m_grid->z())};

  const double power = 1.0 / m_config->get_double("stress_balance.sia.Glen_exponent");
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);

  set_attrs("ice hardness computed using the SIA flow law", "",
            unitstr, unitstr, 0);
}

IceModelVec::Ptr Hardness::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "hardness", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  EnthalpyConverter::Ptr EC = m_grid->ctx()->enthalpy_converter();

  const IceModelVec3  &ice_enthalpy  = model->energy_balance_model()->enthalpy();
  const IceModelVec2S &ice_thickness = model->geometry().ice_thickness;

  const rheology::FlowLaw *flow_law = model->stress_balance()->modifier()->flow_law();

  IceModelVec::AccessList list{&ice_enthalpy, &ice_thickness, result.get()};

  const unsigned int Mz = m_grid->Mz();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      const double *E = ice_enthalpy.get_column(i, j);
      const double H = ice_thickness(i, j);

      double *hardness = result->get_column(i, j);

      for (unsigned int k = 0; k < Mz; ++k) {
        const double depth = H - m_grid->z(k);

        // EC->pressure() handles negative depths correctly
        const double pressure = EC->pressure(depth);

        hardness[k] = flow_law->hardness(E[k], pressure);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

Viscosity::Viscosity(IceModel *m)
  : Diag<IceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "effective_viscosity", m_grid->z())};

  set_attrs("effective viscosity of ice", "",
            "Pascal second", "kPascal second", 0);
  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

static inline double square(double x) {
  return x * x;
}

IceModelVec::Ptr Viscosity::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "effective_viscosity", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  IceModelVec3 W;
  W.create(m_grid, "wvel", WITH_GHOSTS);

  using mask::ice_free;

  EnthalpyConverter::Ptr EC = m_grid->ctx()->enthalpy_converter();

  const rheology::FlowLaw *flow_law = model->stress_balance()->modifier()->flow_law();

  const IceModelVec2S &ice_thickness = model->geometry().ice_thickness;

  const IceModelVec3
    &ice_enthalpy     = model->energy_balance_model()->enthalpy(),
    &U                = model->stress_balance()->velocity_u(),
    &V                = model->stress_balance()->velocity_v(),
    &W_without_ghosts = model->stress_balance()->velocity_w();

  W_without_ghosts.update_ghosts(W);

  const unsigned int Mz = m_grid->Mz();
  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();
  const std::vector<double> &z = m_grid->z();

  const IceModelVec2CellType &mask = model->geometry().cell_type;

  IceModelVec::AccessList list{&U, &V, &W, &ice_enthalpy, &ice_thickness, &mask, result.get()};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *E = ice_enthalpy.get_column(i, j);
      const double H = ice_thickness(i, j);

      const double
        *u   = U.get_column(i, j),
        *u_n = U.get_column(i, j + 1),
        *u_e = U.get_column(i + 1, j),
        *u_s = U.get_column(i, j - 1),
        *u_w = U.get_column(i - 1, j);

      const double
        *v   = V.get_column(i, j),
        *v_n = V.get_column(i, j + 1),
        *v_e = V.get_column(i + 1, j),
        *v_s = V.get_column(i, j - 1),
        *v_w = V.get_column(i - 1, j);

      const double
        *w   = W.get_column(i, j),
        *w_n = W.get_column(i, j + 1),
        *w_e = W.get_column(i + 1, j),
        *w_s = W.get_column(i, j - 1),
        *w_w = W.get_column(i - 1, j);

      StarStencil<int> m = mask.int_star(i, j);
      const unsigned int
        east  = ice_free(m.e) ? 0 : 1,
        west  = ice_free(m.w) ? 0 : 1,
        south = ice_free(m.s) ? 0 : 1,
        north = ice_free(m.n) ? 0 : 1;

      double *viscosity = result->get_column(i, j);

      if (ice_free(m.ij)) {
        result->set_column(i, j, m_fill_value);
        continue;
      }

      for (unsigned int k = 0; k < Mz; ++k) {
        const double depth = H - z[k];

        if (depth < 0.0) {
          viscosity[k] = m_fill_value;
          continue;
        }

        // EC->pressure() handles negative depths correctly
        const double pressure = EC->pressure(depth);

        const double hardness = flow_law->hardness(E[k], pressure);

        double u_x = 0.0, v_x = 0.0, w_x = 0.0;
        if (west + east > 0) {
          const double D = 1.0 / (dx * (west + east));
          u_x = D * (west * (u[k] - u_w[k]) + east * (u_e[k] - u[k]));
          v_x = D * (west * (v[k] - v_w[k]) + east * (v_e[k] - v[k]));
          w_x = D * (west * (w[k] - w_w[k]) + east * (w_e[k] - w[k]));
        }

        double u_y = 0.0, v_y = 0.0, w_y = 0.0;
        if (south + north > 0) {
          const double D = 1.0 / (dy * (south + north));
          u_y = D * (south * (u[k] - u_s[k]) + north * (u_n[k] - u[k]));
          v_y = D * (south * (v[k] - v_s[k]) + north * (v_n[k] - v[k]));
          w_y = D * (south * (w[k] - w_s[k]) + north * (w_n[k] - w[k]));
        }

        double
          u_z = 0.0,
          v_z = 0.0,
          w_z = 0.0;

        if (k == 0) {
          const double dz = z[1] - z[0];
          u_z = (u[1] - u[0]) / dz;
          v_z = (v[1] - v[0]) / dz;
          w_z = (w[1] - w[0]) / dz;
        } else if (k == Mz - 1) {
          const double dz = z[Mz - 1] - z[Mz - 2];
          u_z = (u[Mz - 1] - u[Mz - 2]) / dz;
          v_z = (v[Mz - 1] - v[Mz - 2]) / dz;
          w_z = (w[Mz - 1] - w[Mz - 2]) / dz;
        } else {
          const double
            dz_p = z[k + 1] - z[k],
            dz_m = z[k] - z[k - 1];
          u_z = 0.5 * ((u[k + 1] - u[k]) / dz_p + (u[k] - u[k - 1]) / dz_m);
          v_z = 0.5 * ((v[k + 1] - v[k]) / dz_p + (v[k] - v[k - 1]) / dz_m);
          w_z = 0.5 * ((w[k + 1] - w[k]) / dz_p + (w[k] - w[k - 1]) / dz_m);
        }

        // These should be "epsilon dot", but that's just too long.
        const double
          eps_xx = u_x,
          eps_yy = v_y,
          eps_zz = w_z,
          eps_xy = 0.5 * (u_y + v_x),
          eps_xz = 0.5 * (u_z + w_x),
          eps_yz = 0.5 * (v_z + w_y);

        // The second invariant of the 3D strain rate tensor; see equation 4.8 in [@ref
        // GreveBlatter2009]. Unlike secondInvariant_2D(), this code does not make assumptions about
        // the input velocity field: we do not ignore w_x and w_y and do not assume that u_z and v_z
        // are zero.
        const double
          gamma = (square(eps_xx) + square(eps_yy) + square(eps_zz) +
                   2.0 * (square(eps_xy) + square(eps_xz) + square(eps_yz)));

        double nu = 0.0;
        // Note: in PISM gamma has an extra factor of 1/2; compare to
        flow_law->effective_viscosity(hardness, 0.5 * gamma, &nu, NULL);

        viscosity[k] = nu;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

} // end of namespace diagnostics

void IceModel::init_diagnostics() {

  using namespace diagnostics;

  typedef Diagnostic::Ptr f;   // "f" for "field"
  m_diagnostics = {
    {"bedtoptemp",                          Diagnostic::wrap(m_bedtoptemp)},
    {"cts",                                 f(new CTS(this))},
    {"enthalpybase",                        f(new EnthalpyBasal(this))},
    {"enthalpysurf",                        f(new EnthalpySurface(this))},
    {"hardav",                              f(new HardnessAverage(this))},
    {"hardness",                            f(new Hardness(this))},
    {"liqfrac",                             f(new LiquidFraction(this))},
    {"rank",                                f(new Rank(this))},
    {"temp",                                f(new Temperature(this))},
    {"temp_pa",                             f(new TemperaturePA(this))},
    {"tempbase",                            f(new TemperatureBasal(this))},
    {"tempicethk",                          f(new TemperateIceThickness(this))},
    {"tempicethk_basal",                    f(new TemperateIceThicknessBasal(this))},
    {"temppabase",                          f(new TemperaturePABasal(this))},
    {"tempsurf",                            f(new TemperatureSurface(this))},
    {"dHdt",                                f(new ThicknessRateOfChange(this))},
    {"effective_viscosity",                 f(new Viscosity(this))},
    {"basal_grounded_mass_flux",            f(new BMBSplit(this, GROUNDED))},
    {"basal_floating_mass_flux",            f(new BMBSplit(this, FLOATING))},
    {land_ice_area_fraction_name,           f(new IceAreaFraction(this))},
    {grounded_ice_sheet_area_fraction_name, f(new IceAreaFractionGrounded(this))},
    {floating_ice_sheet_area_fraction_name, f(new IceAreaFractionFloating(this))},
    {"height_above_flotation",              f(new HeightAboveFloatation(this))},
    {"ice_mass",                            f(new IceMass(this))},
    {"topg_sl_adjusted",                    f(new BedTopographySeaLevelAdjusted(this))},
    {"bmelt",                               Diagnostic::wrap(m_basal_melt_rate)},
    {"cell_area",                           Diagnostic::wrap(m_geometry.cell_area)},
    {"lat",                                 Diagnostic::wrap(m_geometry.latitude)},
    {"lon",                                 Diagnostic::wrap(m_geometry.longitude)},
    {"thk",                                 Diagnostic::wrap(m_geometry.ice_thickness)},
    {"ice_area_specific_volume",            Diagnostic::wrap(m_geometry.ice_area_specific_volume)},
    {"mask",                                Diagnostic::wrap(m_geometry.cell_type)},
    {"cell_grounded_fraction",              Diagnostic::wrap(m_geometry.cell_grounded_fraction)},
    {"usurf",                               Diagnostic::wrap(m_geometry.ice_surface_elevation)},
    {"ssa_bc_mask",                         Diagnostic::wrap(m_ssa_dirichlet_bc_mask)},
    {"ssa_bc_vel",                          Diagnostic::wrap(m_ssa_dirichlet_bc_values)},
    {"ocean_pressure_difference",           f(new CalvingFrontPressureDifference(this))},
  };

#if (PISM_USE_PROJ4==1)
  std::string proj4 = m_grid->get_mapping_info().proj4;
  if (not proj4.empty()) {
    m_diagnostics["lat_bnds"] = f(new LatLonBounds(this, "lat", proj4));
    m_diagnostics["lon_bnds"] = f(new LatLonBounds(this, "lon", proj4));
  }
#endif

  typedef TSDiagnostic::Ptr s; // "s" for "scalar"
  m_ts_diagnostics = {
    {"volume_glacierized",                   s(new VolumeGlacierized(this))},
    {"volume_nonglacierized",                s(new VolumeNonGlacierized(this))},
    {"slvol",                                s(new SeaLevelVolume(this))},
    {"volume_rate_of_change_glacierized",    s(new VolumeRateOfChangeGlacierized(this))},
    {"volume_rate_of_change_nonglacierized", s(new VolumeRateOfChangeNonGlacierized(this))},
    {"area_glacierized",                     s(new AreaGlacierized(this))},
    {"mass_glacierized",                     s(new MassGlacierized(this))},
    {"mass_nonglacierized",                  s(new MassNonGlacierized(this))},
    {"mass_rate_of_change_glacierized",      s(new MassRateOfChangeGlacierized(this))},
    {"mass_rate_of_change_nonglacierized",   s(new MassRateOfChangeNonGlacierized(this))},
    {"volume_glacierized_temperate",         s(new VolumeGlacierizedTemperate(this))},
    {"volume_nonglacierized_temperate",      s(new VolumeNonGlacierizedTemperate(this))},
    {"volume_glacierized_cold",              s(new VolumeGlacierizedCold(this))},
    {"volume_nonglacierized_cold",           s(new VolumeNonGlacierizedCold(this))},
    {"volume_glacierized_grounded",          s(new VolumeGlacierizedGrounded(this))},
    {"volume_glacierized_shelf",             s(new VolumeGlacierizedShelf(this))},
    {"area_glacierized_temperate_base",      s(new AreaGlacierizedTemperateBase(this))},
    {"area_glacierized_cold_base",           s(new AreaGlacierizedColdBase(this))},
    {"area_glacierized_grounded",            s(new AreaGlacierizedGrounded(this))},
    {"area_glacierized_shelf",               s(new AreaGlacierizedShelf(this))},
    {"dt",                                   s(new TimeStepLength(this))},
    {"max_diffusivity",                      s(new MaxDiffusivity(this))},
    {"enthalpy_glacierized",                 s(new EnthalpyGlacierized(this))},
    {"enthalpy_nonglacierized",              s(new EnthalpyNonGlacierized(this))},
    {"max_hor_vel",                          s(new MaxHorizontalVelocity(this))},
    {"limnsw",                               s(new MassNotDisplacingSeaWater(this))},
    {"surface_ice_flux",                     s(new MassFluxSurface(this))},
    {"grounded_basal_ice_flux",              s(new MassFluxBasalGrounded(this))},
    {"sub_shelf_ice_flux",                   s(new MassFluxBasalFloating(this))},
    {"discharge_flux",                       s(new MassFluxDischarge(this))},
  };

  // get diagnostics from submodels
  for (auto m : m_submodels) {
    m_diagnostics = pism::combine(m_diagnostics, m.second->diagnostics());
    m_ts_diagnostics = pism::combine(m_ts_diagnostics, m.second->ts_diagnostics());
  }
}

void IceModel::list_diagnostics() {

  m_log->message(1, "\n");

  // 2D and 3D diagnostics
  for (unsigned int d = 3; d > 1; --d) {

    m_log->message(1,
                   "======== Available %dD diagnostic quantities ========\n",
                   d);

    for (auto f : m_diagnostics) {
      Diagnostic::Ptr diag = f.second;

      std::string
        name                = f.first,
        units               = diag->metadata().get_string("units"),
        glaciological_units = diag->metadata().get_string("glaciological_units");

      if (not glaciological_units.empty()) {
        units = glaciological_units;
      }

      if (diag->metadata().get_n_spatial_dimensions() == d) {

        m_log->message(1, "   Name: %s [%s]\n", name.c_str(), units.c_str());

        for (unsigned int k = 0; k < diag->n_variables(); ++k) {
          SpatialVariableMetadata var = diag->metadata(k);

          std::string long_name = var.get_string("long_name");

          m_log->message(1, "      -  %s\n", long_name.c_str());
        }

        m_log->message(1, "\n");
      }
    } // end of the loop over diagnostics
  }

  // scalar time-series
  m_log->message(1, "======== Available time-series ========\n");

  for (auto d : m_ts_diagnostics) {
    const VariableMetadata &m = d.second->metadata();

    std::string
      name                = d.first,
      long_name           = m.get_string("long_name"),
      units               = m.get_string("units"),
      glaciological_units = m.get_string("glaciological_units");

    if (not glaciological_units.empty()) {
      units = glaciological_units;
    }

    m_log->message(1,
                   "   Name: %s [%s]\n"
                   "      -  %s\n\n",
                   name.c_str(), units.c_str(), long_name.c_str());
  }
}

} // end of namespace pism
