// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Constantine Khroulev
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

#include "diagnostics.hh"

#include "pism/age/AgeModel.hh"
#include "pism/energy/EnergyModel.hh"
#include "pism/energy/utilities.hh"
#include "pism/geometry/grounded_cell_fraction.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/SSB_Modifier.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec3Custom.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"
#include "pism/earth/BedDef.hh"

#if (PISM_USE_PROJ4==1)
#include "pism/util/Proj.hh"
#endif

#include "flux_balance.hh"

namespace pism {

// Horrendous names used by InitMIP (and ISMIP6, and CMIP5). Ugh.
static const char* land_ice_area_fraction_name           = "sftgif";
static const char* grounded_ice_sheet_area_fraction_name = "sftgrf";
static const char* floating_ice_sheet_area_fraction_name = "sftflf";

namespace diagnostics {

enum AreaType {GROUNDED, SHELF, BOTH};

enum TermType {SMB, BMB, FLOW, ERROR};

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

  auto
    &H         = model->geometry().ice_thickness,
    &bed       = model->geometry().bed_elevation,
    &sea_level = model->geometry().sea_level_elevation;

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
                                                                   H(i, j), bed(i, j), sea_level(i, j),
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
                            ? "basal_mass_flux_grounded"
                            : "basal_mass_flux_floating",
                            TOTAL_CHANGE), m_kind(flag) {
    assert(flag != BOTH);

    std::string name, description;
    if (m_kind == GROUNDED) {
      name        = "basal_mass_flux_grounded";
      description = "average basal mass flux over the reporting interval (grounded areas)";
    } else {
      name        = "basal_mass_flux_floating";
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

    double ice_density = m_config->get_double("constants.ice.density");

    // the accumulator has the units of kg/m^2, computed as
    //
    // accumulator += BMB (m) * ice_density (kg / m^3)

    IceModelVec::AccessList list{&input, &cell_type, &m_accumulator};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_kind == GROUNDED and cell_type.grounded(i, j)) {
        m_accumulator(i, j) += input(i, j) * ice_density;
      } else if (m_kind == SHELF and cell_type.ocean(i, j)) {
        m_accumulator(i, j) += input(i, j) * ice_density;
      } else {
        m_accumulator(i, j) = 0.0;
      }
    }

    m_interval_length += dt;
  }
};

HardnessAverage::HardnessAverage(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "hardav")};

  // choice to use SSA power; see #285
  const double power = 1.0 / m_config->get_double("stress_balance.ssa.Glen_exponent");
  auto unitstr = pism::printf("Pa s%f", power);

  set_attrs("vertical average of ice hardness", "",
            unitstr, unitstr, 0);

  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

//! \brief Computes vertically-averaged ice hardness.
IceModelVec::Ptr HardnessAverage::compute_impl() const {

  const rheology::FlowLaw *flow_law = model->stress_balance()->shallow()->flow_law().get();
  if (flow_law == NULL) {
    flow_law = model->stress_balance()->modifier()->flow_law().get();
    if (flow_law == NULL) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "Can't compute vertically-averaged hardness: no flow law is used.");
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
            "1", "1", 0);
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

IceEnthalpySurface::IceEnthalpySurface(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "enthalpysurf")};

  set_attrs("ice enthalpy at 1m below the ice surface", "",
            "J kg-1", "J kg-1", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceEnthalpySurface::compute_impl() const {

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

IceEnthalpyBasal::IceEnthalpyBasal(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "enthalpybase")};

  set_attrs("ice enthalpy at the base of ice", "",
            "J kg-1", "J kg-1", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceEnthalpyBasal::compute_impl() const {

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

  IceModelVec::Ptr enth = IceEnthalpyBasal(model).compute();

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(enth);

  // result contains basal enthalpy; note that it is allocated by
  // IceEnthalpyBasal::compute().

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

  IceModelVec::Ptr enth = IceEnthalpySurface(model).compute();
  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(enth);

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  // result contains surface enthalpy; note that it is allocated by
  // IceEnthalpySurface::compute().

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

namespace scalar {

//! \brief Computes the total ice volume in glacierized areas.
class IceVolumeGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceVolumeGlacierized(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized") {

    m_ts.variable().set_string("units", "m3");
    m_ts.variable().set_string("long_name", "volume of the ice in glacierized areas");
    m_ts.variable().set_double("valid_min", 0.0);
  }
  double compute() {
    return model->ice_volume(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total ice volume.
class IceVolume : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceVolume(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume") {

    m_ts.variable().set_string("units", "m3");
    m_ts.variable().set_string("long_name", "volume of the ice, including seasonal cover");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_volume(0.0);
  }
};

//! \brief Computes the total ice volume which is relevant for sea-level
class SeaLevelRisePotential : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  SeaLevelRisePotential(const IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "sea_level_rise_potential") {

    m_ts.variable().set_string("units", "m");
    m_ts.variable().set_string("long_name", "the sea level rise that would result if all the ice were melted");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->sea_level_rise_potential(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the rate of change of the total ice volume in glacierized areas.
class IceVolumeRateOfChangeGlacierized : public TSDiag<TSRateDiagnostic, IceModel>
{
public:
  IceVolumeRateOfChangeGlacierized(IceModel *m)
    : TSDiag<TSRateDiagnostic, IceModel>(m, "tendency_of_ice_volume_glacierized") {

    m_ts.variable().set_string("units", "m3 s-1");
    m_ts.variable().set_string("glaciological_units", "m3 year-1");
    m_ts.variable().set_string("long_name", "rate of change of the ice volume in glacierized areas");
  }

  double compute() {
    return model->ice_volume(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the rate of change of the total ice volume.
class IceVolumeRateOfChange : public TSDiag<TSRateDiagnostic, IceModel>
{
public:
  IceVolumeRateOfChange(IceModel *m)
    : TSDiag<TSRateDiagnostic, IceModel>(m, "tendency_of_ice_volume") {

    m_ts.variable().set_string("units", "m3 s-1");
    m_ts.variable().set_string("glaciological_units", "m3 year-1");
    m_ts.variable().set_string("long_name",
                               "rate of change of the ice volume, including seasonal cover");
  }

  double compute() {
    return model->ice_volume(0.0);
  }
};

//! \brief Computes the total ice area.
class IceAreaGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceAreaGlacierized(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized") {

    m_ts.variable().set_string("units", "m2");
    m_ts.variable().set_string("long_name", "glacierized area");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_area(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total mass of the ice not displacing sea water.
class IceMassNotDisplacingSeaWater : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceMassNotDisplacingSeaWater(const IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "limnsw") {

    m_ts.variable().set_string("units", "kg");
    m_ts.variable().set_string("long_name", "mass of the ice not displacing sea water");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {

    const double
      thickness_standard = m_config->get_double("output.ice_free_thickness_standard"),
      ice_density        = m_config->get_double("constants.ice.density"),
      ice_volume         = model->ice_volume_not_displacing_seawater(thickness_standard),
      ice_mass           = ice_volume * ice_density;

    return ice_mass;
  }
};

//! \brief Computes the total ice mass in glacierized areas.
class IceMassGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceMassGlacierized(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_mass_glacierized") {

    m_ts.variable().set_string("units", "kg");
    m_ts.variable().set_string("long_name", "mass of the ice in glacierized areas");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    double
      ice_density        = m_config->get_double("constants.ice.density"),
      thickness_standard = m_config->get_double("output.ice_free_thickness_standard");
    return model->ice_volume(thickness_standard) * ice_density;
  }
};

//! \brief Computes the total ice mass.
class IceMass : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceMass(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_mass") {

    m_ts.variable().set_string("units", "kg");
    m_ts.variable().set_string("long_name", "mass of the ice, including seasonal cover");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return (model->ice_volume(0.0) *
            m_config->get_double("constants.ice.density"));
  }
};

//! \brief Computes the rate of change of the total ice mass in glacierized areas.
class IceMassRateOfChangeGlacierized : public TSDiag<TSRateDiagnostic, IceModel>
{
public:
  IceMassRateOfChangeGlacierized(IceModel *m)
    : TSDiag<TSRateDiagnostic, IceModel>(m, "tendency_of_ice_mass_glacierized") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "Gt year-1");
    m_ts.variable().set_string("long_name", "rate of change of the ice mass in glacierized areas");
  }

  double compute() {
    double
      ice_density         = m_config->get_double("constants.ice.density"),
      thickness_threshold = m_config->get_double("output.ice_free_thickness_standard");
    return model->ice_volume(thickness_threshold) * ice_density;
  }
};

//! \brief Computes the rate of change of the total ice mass due to flow (influx due to
//! prescribed constant-in-time ice thickness).
/*!
 * This is the change in mass resulting from prescribing (fixing) ice thickness.
 */
class IceMassRateOfChangeDueToFlow : public TSDiag<TSFluxDiagnostic, IceModel>
{
public:
  IceMassRateOfChangeDueToFlow(IceModel *m)
    : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_flow") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "Gt year-1");
    m_ts.variable().set_string("long_name", "rate of change of the mass of ice due to flow"
                               " (i.e. prescribed ice thickness)");
  }

  double compute() {

    const double
      ice_density = m_config->get_double("constants.ice.density");

    const IceModelVec2S
      &dH = model->geometry_evolution().thickness_change_due_to_flow(),
      &dV = model->geometry_evolution().area_specific_volume_change_due_to_flow();

    auto cell_area = m_grid->cell_area();

    IceModelVec::AccessList list{&dH, &dV};

    double volume_change = 0.0;
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      // m * m^2 = m^3
      volume_change += (dH(i, j) + dV(i, j)) * cell_area;
    }

    // (kg/m^3) * m^3 = kg
    return ice_density * GlobalSum(m_grid->com, volume_change);
  }
};

//! \brief Computes the rate of change of the total ice mass.
class IceMassRateOfChange : public TSDiag<TSRateDiagnostic, IceModel>
{
public:
  IceMassRateOfChange(IceModel *m)
    : TSDiag<TSRateDiagnostic, IceModel>(m, "tendency_of_ice_mass") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "Gt year-1");
    m_ts.variable().set_string("long_name",
                               "rate of change of the mass of ice, including seasonal cover");
  }

  double compute() {
    const double ice_density = m_config->get_double("constants.ice.density");
    return model->ice_volume(0.0) * ice_density;
  }
};


//! \brief Computes the total volume of the temperate ice in glacierized areas.
class IceVolumeGlacierizedTemperate : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceVolumeGlacierizedTemperate(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized_temperate") {

    m_ts.variable().set_string("units", "m3");
    m_ts.variable().set_string("long_name", "volume of temperate ice in glacierized areas");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_volume_temperate(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total volume of the temperate ice.
class IceVolumeTemperate : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceVolumeTemperate(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_temperate") {

    m_ts.variable().set_string("units", "m3");
    m_ts.variable().set_string("long_name", "volume of temperate ice, including seasonal cover");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_volume_temperate(0.0);
  }
};

//! \brief Computes the total volume of the cold ice in glacierized areas.
class IceVolumeGlacierizedCold : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceVolumeGlacierizedCold(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized_cold") {

    m_ts.variable().set_string("units", "m3");
    m_ts.variable().set_string("long_name", "volume of cold ice in glacierized areas");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_volume_cold(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total volume of the cold ice.
class IceVolumeCold : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceVolumeCold(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_cold") {

    m_ts.variable().set_string("units", "m3");
    m_ts.variable().set_string("long_name", "volume of cold ice, including seasonal cover");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_volume_cold(0.0);
  }
};

//! \brief Computes the total area of the temperate ice.
class IceAreaGlacierizedTemperateBase : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceAreaGlacierizedTemperateBase(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized_temperate_base") {

    m_ts.variable().set_string("units", "m2");
    m_ts.variable().set_string("long_name", "glacierized area where basal ice is temperate");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_area_temperate(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total area of the cold ice.
class IceAreaGlacierizedColdBase : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceAreaGlacierizedColdBase(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized_cold_base") {

    m_ts.variable().set_string("units", "m2");
    m_ts.variable().set_string("long_name", "glacierized area where basal ice is cold");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_area_cold(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total ice enthalpy in glacierized areas.
class IceEnthalpyGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceEnthalpyGlacierized(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_enthalpy_glacierized") {

    m_ts.variable().set_string("units", "J");
    m_ts.variable().set_string("long_name", "enthalpy of the ice in glacierized areas");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return energy::total_ice_enthalpy(m_config->get_double("output.ice_free_thickness_standard"),
                                      model->energy_balance_model()->enthalpy(),
                                      model->geometry().ice_thickness);
  }
};

//! \brief Computes the total ice enthalpy.
class IceEnthalpy : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceEnthalpy(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_enthalpy") {

    m_ts.variable().set_string("units", "J");
    m_ts.variable().set_string("long_name", "enthalpy of the ice, including seasonal cover");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return energy::total_ice_enthalpy(0.0,
                                      model->energy_balance_model()->enthalpy(),
                                      model->geometry().ice_thickness);
  }
};

//! \brief Computes the total grounded ice area.
class IceAreaGlacierizedGrounded : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceAreaGlacierizedGrounded(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized_grounded") {

    m_ts.variable().set_string("units", "m2");
    m_ts.variable().set_string("long_name", "area of grounded ice in glacierized areas");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_area_grounded(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total floating ice area.
class IceAreaGlacierizedShelf : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceAreaGlacierizedShelf(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized_floating") {

    m_ts.variable().set_string("units", "m2");
    m_ts.variable().set_string("long_name", "area of ice shelves in glacierized areas");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->ice_area_floating(m_config->get_double("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total grounded ice volume.
class IceVolumeGlacierizedGrounded : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceVolumeGlacierizedGrounded(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized_grounded") {

    m_ts.variable().set_string("units", "m3");
    m_ts.variable().set_string("long_name", "volume of grounded ice in glacierized areas");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    const IceModelVec2CellType &cell_type = model->geometry().cell_type;

    const IceModelVec2S &ice_thickness = model->geometry().ice_thickness;

    const double
      thickness_threshold = m_config->get_double("output.ice_free_thickness_standard"),
      cell_area           = m_grid->cell_area();

    IceModelVec::AccessList list{&ice_thickness, &cell_type};

    double volume = 0.0;
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double H = ice_thickness(i, j);

      if (cell_type.grounded(i, j) and H >= thickness_threshold) {
        volume += cell_area * H;
      }
    }

    return GlobalSum(m_grid->com, volume);
  }
};

//! \brief Computes the total floating ice volume.
class IceVolumeGlacierizedShelf : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  IceVolumeGlacierizedShelf(IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized_floating") {

    m_ts.variable().set_string("units", "m3");
    m_ts.variable().set_string("long_name", "volume of ice shelves in glacierized areas");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    const IceModelVec2CellType &cell_type = model->geometry().cell_type;

    const IceModelVec2S &ice_thickness = model->geometry().ice_thickness;

    const double
      thickness_threshold = m_config->get_double("output.ice_free_thickness_standard"),
      cell_area           = m_grid->cell_area();

    IceModelVec::AccessList list{&ice_thickness, &cell_type};

    double volume = 0.0;
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double H = ice_thickness(i, j);

      if (cell_type.ocean(i, j) and H >= thickness_threshold) {
        volume += cell_area * H;
      }
    }

    return GlobalSum(m_grid->com, volume);
  }
};

//! \brief Reports the mass continuity time step.
class TimeStepLength : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  TimeStepLength(const IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "dt") {

    m_ts.variable().set_string("units", "second");
    m_ts.variable().set_string("glaciological_units", "year");
    m_ts.variable().set_string("long_name", "mass continuity time step");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->dt();
  }
};

//! \brief Reports maximum diffusivity.
class MaxDiffusivity : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  MaxDiffusivity(const IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "max_diffusivity") {

    m_ts.variable().set_string("units", "m2 s-1");
    m_ts.variable().set_string("long_name", "maximum diffusivity");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    return model->stress_balance()->max_diffusivity();
  }
};

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
class MaxHorizontalVelocity : public TSDiag<TSSnapshotDiagnostic, IceModel>
{
public:
  MaxHorizontalVelocity(const IceModel *m)
    : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "max_hor_vel") {

    m_ts.variable().set_string("units", "m second-1");
    m_ts.variable().set_string("glaciological_units", "m year-1");
    m_ts.variable().set_string("long_name",
                               "maximum abs component of horizontal ice velocity"
                               " over grid in last time step during time-series reporting interval");
    m_ts.variable().set_double("valid_min", 0.0);
  }

  double compute() {
    CFLData cfl = model->stress_balance()->max_timestep_cfl_3d();

    return std::max(cfl.u_max, cfl.v_max);
  }
};

/*!
 * Return total mass change due to one of the terms in the mass continuity equation.
 *
 * Possible terms are
 *
 * - SMB: surface mass balance
 * - BMB: basal mass balance
 * - FLOW: ice flow
 * - ERROR: numerical flux needed to preserve non-negativity of thickness
 *
 * This computation can be restricted to grounded and floating areas
 * using the `area` argument.
 *
 * - BOTH: include all contributions
 * - GROUNDED: include grounded areas only
 * - SHELF: include floating areas only
 *
 * When computing mass changes due to flow it is important to remember
 * that ice mass in a cell can be represented by its thickness *or* an
 * "area specific volume". Transferring mass from one representation
 * to the other does not change the mass in a cell. This explains the
 * special case used when `term == FLOW`. (Note that surface and basal
 * mass balances do not affect the area specific volume field.)
 */
double mass_change(const IceModel *model, TermType term, AreaType area) {
  const IceGrid &grid = *model->grid();
  const Config &config = *grid.ctx()->config();

  const double
    ice_density = config.get_double("constants.ice.density"),
    cell_area   = grid.cell_area();

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  const IceModelVec2S *thickness_change = nullptr;

  switch (term) {
  case FLOW:
    thickness_change = &model->geometry_evolution().thickness_change_due_to_flow();
    break;
  case SMB:
    thickness_change = &model->geometry_evolution().top_surface_mass_balance();
    break;
  case BMB:
    thickness_change = &model->geometry_evolution().bottom_surface_mass_balance();
    break;
  case ERROR:
    thickness_change = &model->geometry_evolution().conservation_error();
    break;
  default:
    // can't happen
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid term type");
  }

  const IceModelVec2S &dV_flow = model->geometry_evolution().area_specific_volume_change_due_to_flow();

  IceModelVec::AccessList list{&cell_type, thickness_change};

  if (term == FLOW) {
    list.add(dV_flow);
  }

  double volume_change = 0.0;
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((area == BOTH) or
        (area == GROUNDED and cell_type.grounded(i, j)) or
        (area == SHELF and cell_type.ocean(i, j))) {

      double dV = term == FLOW ? dV_flow(i, j) : 0.0;

      // m^3 = m^2 * m
      volume_change += cell_area * ((*thickness_change)(i, j) + dV);
    }
  }

  // (kg / m^3) * m^3 = kg
  return ice_density * GlobalSum(grid.com, volume_change);
}

//! \brief Reports the total bottom surface ice flux.
class IceMassFluxBasal : public TSDiag<TSFluxDiagnostic, IceModel>
{
public:
  IceMassFluxBasal(const IceModel *m)
    : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_basal_mass_flux") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "Gt year-1");
    m_ts.variable().set_string("long_name", "total over ice domain of bottom surface ice mass flux");
    m_ts.variable().set_string("comment", "positive means ice gain");
  }

  double compute() {
    return mass_change(model, BMB, BOTH);
  }
};

//! \brief Reports the total top surface ice flux.
class IceMassFluxSurface : public TSDiag<TSFluxDiagnostic, IceModel>
{
public:
  IceMassFluxSurface(const IceModel *m)
    : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_surface_mass_flux") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "Gt year-1");
    m_ts.variable().set_string("long_name", "total over ice domain of top surface ice mass flux");
    m_ts.variable().set_string("comment", "positive means ice gain");
  }

  double compute() {
    return mass_change(model, SMB, BOTH);
  }
};

//! \brief Reports the total basal ice flux over the grounded region.
class IceMassFluxBasalGrounded : public TSDiag<TSFluxDiagnostic, IceModel>
{
public:
  IceMassFluxBasalGrounded(const IceModel *m)
    : TSDiag<TSFluxDiagnostic, IceModel>(m, "basal_mass_flux_grounded") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "Gt year-1");
    m_ts.variable().set_string("long_name", "total over grounded ice domain of basal mass flux");
    m_ts.variable().set_string("comment", "positive means ice gain");
  }

  double compute() {
    return mass_change(model, BMB, GROUNDED);
  }
};

//! \brief Reports the total sub-shelf ice flux.
class IceMassFluxBasalFloating : public TSDiag<TSFluxDiagnostic, IceModel>
{
public:
  IceMassFluxBasalFloating(const IceModel *m)
    : TSDiag<TSFluxDiagnostic, IceModel>(m, "basal_mass_flux_floating") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "Gt year-1");
    m_ts.variable().set_string("long_name", "total sub-shelf ice flux");
    m_ts.variable().set_string("comment", "positive means ice gain");
  }

  double compute() {
    return mass_change(model, BMB, SHELF);
  }
};

//! \brief Reports the total numerical mass flux needed to preserve
//! non-negativity of ice thickness.
class IceMassFluxConservationError : public TSDiag<TSFluxDiagnostic, IceModel>
{
public:
  IceMassFluxConservationError(const IceModel *m)
    : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_conservation_error") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "Gt year-1");
    m_ts.variable().set_string("long_name", "total numerical flux needed to preserve non-negativity"
                               " of ice thickness");
    m_ts.variable().set_string("comment", "positive means ice gain");
  }

  double compute() {
    return mass_change(model, ERROR, BOTH);
  }
};

//! \brief Reports the total discharge flux.
class IceMassFluxDischarge : public TSDiag<TSFluxDiagnostic, IceModel>
{
public:
  IceMassFluxDischarge(const IceModel *m)
    : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_discharge") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "Gt year-1");
    m_ts.variable().set_string("long_name", "discharge (calving & icebergs) flux");
    m_ts.variable().set_string("comment", "positive means ice gain");
  }

  double compute() {
    const double ice_density = m_config->get_double("constants.ice.density");

    const IceModelVec2S &discharge = model->discharge();

    auto cell_area = m_grid->cell_area();

    double volume_change = 0.0;

    IceModelVec::AccessList list{&discharge};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      // m^2 * m = m^3
      volume_change += cell_area * discharge(i, j);
    }

    // (kg/m^3) * m^3 = kg
    return ice_density * GlobalSum(m_grid->com, volume_change);
  }
};

} // end of namespace scalar


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

    std::string
      internal_units = "m second-1",
      external_units = "m year-1";

    set_attrs("ice thickness rate of change",
              "tendency_of_land_ice_thickness",
              internal_units, external_units, 0);

    units::Converter c(m_sys, external_units, internal_units);

    const double valid_range = c(1e6);

    m_vars[0].set_doubles("valid_range",  {-valid_range, valid_range});
    m_vars[0].set_double("_FillValue", c(m_fill_value));
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
    m_vars[0].set_doubles("valid_range", {-180, 180});
  } else {
    set_attrs("latitude bounds", "", "degree_north", "degree_north", 0);
    m_vars[0].set_doubles("valid_range", {-90, 90});
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

  const double
    ice_density   = m_config->get_double("constants.ice.density"),
    ocean_density = m_config->get_double("constants.sea_water.density");

  auto
    &ice_thickness  = model->geometry().ice_thickness,
    &sea_level      = model->geometry().sea_level_elevation,
    &bed_topography = model->geometry().bed_elevation;

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  compute_grounded_cell_fraction(ice_density, ocean_density, sea_level,
                                 ice_thickness, bed_topography,
                                 *result);

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

  set_attrs("ice thickness in excess of the maximum floating ice thickness",
            "", "m", "m", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
  m_vars[0].set_string("comment",
                       "shows how close to floatation the ice is at a given location");
}

IceModelVec::Ptr HeightAboveFloatation::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "height_above_flotation", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  const IceModelVec2CellType &cell_type = model->geometry().cell_type;

  const double
    ice_density   = m_config->get_double("constants.ice.density"),
    ocean_density = m_config->get_double("constants.sea_water.density");

  auto
    &sea_level      = model->geometry().sea_level_elevation,
    &ice_thickness  = model->geometry().ice_thickness,
    &bed_topography = model->geometry().bed_elevation;

  IceModelVec::AccessList list{&cell_type, result.get(), &ice_thickness, &bed_topography, &sea_level};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double
        thickness   = ice_thickness(i, j),
        bed         = bed_topography(i, j),
        ocean_depth = sea_level(i, j) - bed;

      if (cell_type.icy(i, j) and ocean_depth > 0.0) {
        const double max_floating_thickness = ocean_depth * (ocean_density / ice_density);
        (*result)(i, j) = thickness - max_floating_thickness;
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
    &ice_thickness = model->geometry().ice_thickness;

  auto cell_area = m_grid->cell_area();

  IceModelVec::AccessList list{&cell_type, result.get(), &ice_thickness};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // count all ice, including cells which have so little they
      // are considered "ice-free"
      if (ice_thickness(i, j) > 0.0) {
        (*result)(i,j) = ice_density * ice_thickness(i, j) * cell_area;
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
        (*result)(i, j) = ice_density * Href(i, j) * cell_area;
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

  auto
    &bed       = model->geometry().bed_elevation,
    &sea_level = model->geometry().sea_level_elevation;

  IceModelVec::AccessList list{&bed, &sea_level, result.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i, j) = bed(i, j) - sea_level(i, j);
  }

  return result;
}

IceHardness::IceHardness(const IceModel *m)
  : Diag<IceModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "hardness", m_grid->z())};

  const double power = 1.0 / m_config->get_double("stress_balance.sia.Glen_exponent");
  auto unitstr = pism::printf("Pa s%f", power);

  set_attrs("ice hardness computed using the SIA flow law", "",
            unitstr, unitstr, 0);
}

IceModelVec::Ptr IceHardness::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "hardness", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  EnthalpyConverter::Ptr EC = m_grid->ctx()->enthalpy_converter();

  const IceModelVec3  &ice_enthalpy  = model->energy_balance_model()->enthalpy();
  const IceModelVec2S &ice_thickness = model->geometry().ice_thickness;

  const rheology::FlowLaw &flow_law = *model->stress_balance()->modifier()->flow_law();

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

        hardness[k] = flow_law.hardness(E[k], pressure);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

IceViscosity::IceViscosity(IceModel *m)
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

IceModelVec::Ptr IceViscosity::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "effective_viscosity", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  IceModelVec3 W;
  W.create(m_grid, "wvel", WITH_GHOSTS);

  using mask::ice_free;

  EnthalpyConverter::Ptr EC = m_grid->ctx()->enthalpy_converter();

  const rheology::FlowLaw &flow_law = *model->stress_balance()->modifier()->flow_law();

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

        const double hardness = flow_law.hardness(E[k], pressure);

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
        flow_law.effective_viscosity(hardness, 0.5 * gamma, &nu, NULL);

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

  typedef Diagnostic      d;
  typedef Diagnostic::Ptr f;    // "f" for "field"
  m_diagnostics = {
    // geometry
    {"cell_grounded_fraction",              d::wrap(m_geometry.cell_grounded_fraction)},
    {"height_above_flotation",              f(new HeightAboveFloatation(this))},
    {"ice_area_specific_volume",            d::wrap(m_geometry.ice_area_specific_volume)},
    {"ice_mass",                            f(new IceMass(this))},
    {"lat",                                 d::wrap(m_geometry.latitude)},
    {"lon",                                 d::wrap(m_geometry.longitude)},
    {"mask",                                d::wrap(m_geometry.cell_type)},
    {"thk",                                 d::wrap(m_geometry.ice_thickness)},
    {"topg_sl_adjusted",                    f(new BedTopographySeaLevelAdjusted(this))},
    {"usurf",                               d::wrap(m_geometry.ice_surface_elevation)},
    {floating_ice_sheet_area_fraction_name, f(new IceAreaFractionFloating(this))},
    {grounded_ice_sheet_area_fraction_name, f(new IceAreaFractionGrounded(this))},
    {land_ice_area_fraction_name,           f(new IceAreaFraction(this))},

    // temperature, enthalpy, and liquid water fraction
    {"enthalpybase", f(new IceEnthalpyBasal(this))},
    {"enthalpysurf", f(new IceEnthalpySurface(this))},
    {"bedtoptemp",   d::wrap(m_bedtoptemp)},
    {"cts",          f(new CTS(this))},
    {"liqfrac",      f(new LiquidFraction(this))},
    {"temp",         f(new Temperature(this))},
    {"temp_pa",      f(new TemperaturePA(this))},
    {"tempbase",     f(new TemperatureBasal(this))},
    {"temppabase",   f(new TemperaturePABasal(this))},
    {"tempsurf",     f(new TemperatureSurface(this))},

    // rheology-related stuff
    {"tempicethk",          f(new TemperateIceThickness(this))},
    {"tempicethk_basal",    f(new TemperateIceThicknessBasal(this))},
    {"hardav",              f(new HardnessAverage(this))},
    {"hardness",            f(new IceHardness(this))},
    {"effective_viscosity", f(new IceViscosity(this))},

    // boundary conditions
    {"ssa_bc_mask",               d::wrap(m_ssa_dirichlet_bc_mask)},
    {"ssa_bc_vel",                d::wrap(m_ssa_dirichlet_bc_values)},
    {"ocean_pressure_difference", f(new CalvingFrontPressureDifference(this))},

    // balancing the books
    // tendency_of_ice_amount = (tendency_of_ice_amount_due_to_flow +
    //                           tendency_of_ice_amount_due_to_conservation_error +
    //                           tendency_of_ice_amount_due_to_surface_mass_balance +
    //                           tendency_of_ice_amount_due_to_basal_mass_balance +
    //                           tendency_of_ice_amount_due_to_discharge)
    {"tendency_of_ice_amount",                           f(new TendencyOfIceAmount(this,          AMOUNT))},
    {"tendency_of_ice_amount_due_to_flow",               f(new TendencyOfIceAmountDueToFlow(this, AMOUNT))},
    {"tendency_of_ice_amount_due_to_conservation_error", f(new ConservationErrorFlux(this,        AMOUNT))},
    {"tendency_of_ice_amount_due_to_surface_mass_flux",  f(new SurfaceFlux(this,                  AMOUNT))},
    {"tendency_of_ice_amount_due_to_basal_mass_flux",    f(new BasalFlux(this,                    AMOUNT))},
    {"tendency_of_ice_amount_due_to_discharge",          f(new DischargeFlux(this,                AMOUNT))},

    // same, in terms of mass
    // tendency_of_ice_mass = (tendency_of_ice_mass_due_to_flow +
    //                         tendency_of_ice_mass_due_to_conservation_error +
    //                         tendency_of_ice_mass_due_to_surface_mass_flux +
    //                         tendency_of_ice_mass_due_to_basal_mass_balance +
    //                         tendency_of_ice_mass_due_to_discharge)
    {"tendency_of_ice_mass",                           f(new TendencyOfIceAmount(this,          MASS))},
    {"tendency_of_ice_mass_due_to_flow",               f(new TendencyOfIceAmountDueToFlow(this, MASS))},
    {"tendency_of_ice_mass_due_to_conservation_error", f(new ConservationErrorFlux(this,        MASS))},
    {"tendency_of_ice_mass_due_to_surface_mass_flux",  f(new SurfaceFlux(this,                  MASS))},
    {"tendency_of_ice_mass_due_to_basal_mass_flux",    f(new BasalFlux(this,                    MASS))},
    {"tendency_of_ice_mass_due_to_discharge",          f(new DischargeFlux(this,                MASS))},

    // other rates and fluxes
    {"basal_mass_flux_grounded", f(new BMBSplit(this, GROUNDED))},
    {"basal_mass_flux_floating",    f(new BMBSplit(this, SHELF))},
    {"dHdt",                     f(new ThicknessRateOfChange(this))},
    {"bmelt",                    d::wrap(m_basal_melt_rate)},

    // misc
    {"rank", f(new Rank(this))},
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
    // area
    {"ice_area_glacierized",                s(new scalar::IceAreaGlacierized(this))},
    {"ice_area_glacierized_cold_base",      s(new scalar::IceAreaGlacierizedColdBase(this))},
    {"ice_area_glacierized_grounded",       s(new scalar::IceAreaGlacierizedGrounded(this))},
    {"ice_area_glacierized_floating",       s(new scalar::IceAreaGlacierizedShelf(this))},
    {"ice_area_glacierized_temperate_base", s(new scalar::IceAreaGlacierizedTemperateBase(this))},
    // mass
    {"ice_mass_glacierized",             s(new scalar::IceMassGlacierized(this))},
    {"ice_mass",                         s(new scalar::IceMass(this))},
    {"tendency_of_ice_mass_glacierized", s(new scalar::IceMassRateOfChangeGlacierized(this))},
    {"limnsw",                           s(new scalar::IceMassNotDisplacingSeaWater(this))},
    // volume
    {"ice_volume_glacierized",             s(new scalar::IceVolumeGlacierized(this))},
    {"ice_volume_glacierized_cold",        s(new scalar::IceVolumeGlacierizedCold(this))},
    {"ice_volume_glacierized_grounded",    s(new scalar::IceVolumeGlacierizedGrounded(this))},
    {"ice_volume_glacierized_floating",    s(new scalar::IceVolumeGlacierizedShelf(this))},
    {"ice_volume_glacierized_temperate",   s(new scalar::IceVolumeGlacierizedTemperate(this))},
    {"ice_volume",                         s(new scalar::IceVolume(this))},
    {"ice_volume_cold",                    s(new scalar::IceVolumeCold(this))},
    {"ice_volume_temperate",               s(new scalar::IceVolumeTemperate(this))},
    {"tendency_of_ice_volume_glacierized", s(new scalar::IceVolumeRateOfChangeGlacierized(this))},
    {"tendency_of_ice_volume",             s(new scalar::IceVolumeRateOfChange(this))},
    {"sea_level_rise_potential",           s(new scalar::SeaLevelRisePotential(this))},
    // energy
    {"ice_enthalpy_glacierized", s(new scalar::IceEnthalpyGlacierized(this))},
    {"ice_enthalpy",         s(new scalar::IceEnthalpy(this))},
    // time-stepping
    {"max_diffusivity", s(new scalar::MaxDiffusivity(this))},
    {"max_hor_vel",     s(new scalar::MaxHorizontalVelocity(this))},
    {"dt",              s(new scalar::TimeStepLength(this))},
    // balancing the books
    {"tendency_of_ice_mass",                           s(new scalar::IceMassRateOfChange(this))},
    {"tendency_of_ice_mass_due_to_flow",               s(new scalar::IceMassRateOfChangeDueToFlow(this))},
    {"tendency_of_ice_mass_due_to_conservation_error", s(new scalar::IceMassFluxConservationError(this))},
    {"tendency_of_ice_mass_due_to_basal_mass_flux",    s(new scalar::IceMassFluxBasal(this))},
    {"tendency_of_ice_mass_due_to_surface_mass_flux",  s(new scalar::IceMassFluxSurface(this))},
    {"tendency_of_ice_mass_due_to_discharge",          s(new scalar::IceMassFluxDischarge(this))},
    // other fluxes
    {"basal_mass_flux_grounded", s(new scalar::IceMassFluxBasalGrounded(this))},
    {"basal_mass_flux_floating", s(new scalar::IceMassFluxBasalFloating(this))},
  };

  // get diagnostics from submodels
  for (auto m : m_submodels) {
    m_diagnostics = pism::combine(m_diagnostics, m.second->diagnostics());
    m_ts_diagnostics = pism::combine(m_ts_diagnostics, m.second->ts_diagnostics());
  }
}

typedef std::map<std::string, std::vector<VariableMetadata>> Metadata;

static void print_diagnostics(const Logger &log, const Metadata &list) {
  for (const auto &d : list) {
    const std::string &name = d.first;
    log.message(1, " Name: %s\n", name.c_str());

    for (const auto &v : d.second) {

      std::string
        var_name            = v.get_name(),
        units               = v.get_string("units"),
        glaciological_units = v.get_string("glaciological_units"),
        long_name           = v.get_string("long_name"),
        comment             = v.get_string("comment");

      if (not glaciological_units.empty()) {
        units = glaciological_units;
      }

      log.message(1, "   %s [%s]\n", var_name.c_str(), units.c_str());
      log.message(1, "    %s\n", long_name.c_str());
      if (not comment.empty()) {
        log.message(1, "    %s\n", comment.c_str());
      }
    }
    log.message(1, "\n");
  }
}

static void print_diagnostics_json(const Logger &log, const Metadata &list) {
  log.message(1, "{\n");
  bool first_diagnostic = true;
  for (const auto &d : list) {

    if (not first_diagnostic) {
      log.message(1, ",\n");
    } else {
      first_diagnostic = false;
    }

    log.message(1, "\"%s\" : [\n", d.first.c_str());

    bool first_variable = true;
    for (const auto &variable : d.second) {

      std::string
        var_name            = variable.get_name(),
        units               = variable.get_string("units"),
        glaciological_units = variable.get_string("glaciological_units"),
        long_name           = variable.get_string("long_name"),
        standard_name       = variable.get_string("standard_name"),
        comment             = variable.get_string("comment");

      if (not glaciological_units.empty()) {
        units = glaciological_units;
      }

      if (not first_variable) {
        log.message(1, ",\n");
      } else {
        first_variable = false;
      }

      log.message(1, "[\"%s\", \"%s\", \"%s\", \"%s\", \"%s\"]",
                  var_name.c_str(), units.c_str(), long_name.c_str(), standard_name.c_str(), comment.c_str());
    }
    log.message(1, "]");
  }
  log.message(1, "}\n");
}

/*!
 * Return metadata of 2D and 3D diagnostics.
 */
static Metadata diag_metadata(const std::map<std::string,Diagnostic::Ptr> &diags) {
  Metadata result;

  for (auto f : diags) {
    Diagnostic::Ptr diag = f.second;

    for (unsigned int k = 0; k < diag->n_variables(); ++k) {
      result[f.first].push_back(diag->metadata(k));
    }
  }

  return result;
}

/*!
 * Return metadata of scalar diagnostics.
 */
static Metadata ts_diag_metadata(const std::map<std::string,TSDiagnostic::Ptr> &ts_diags) {
  Metadata result;

  for (auto d : ts_diags) {
    // always one variable per diagnostic
    result[d.first] = {d.second->metadata()};
  }

  return result;
}

void IceModel::list_diagnostics_json() const {

  m_log->message(1, "{\n");

  m_log->message(1, "\"spatial\" :\n");
  print_diagnostics_json(*m_log, diag_metadata(m_diagnostics));

  m_log->message(1, ",\n");        // separator

  m_log->message(1, "\"scalar\" :\n");
  print_diagnostics_json(*m_log, ts_diag_metadata(m_ts_diagnostics));

  m_log->message(1, "}\n");
}

void IceModel::list_diagnostics() const {

  m_log->message(1, "\n");
  m_log->message(1, "======== Available 2D and 3D diagnostics ========\n");

  print_diagnostics(*m_log, diag_metadata(m_diagnostics));

  // scalar time-series
  m_log->message(1, "======== Available time-series ========\n");

  print_diagnostics(*m_log, ts_diag_metadata(m_ts_diagnostics));
}

/*!
  Computes fraction of the base which is melted.

  Communication occurs here.
 */
double IceModel::compute_temperate_base_fraction(double total_ice_area) {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  auto cell_area = m_grid->cell_area();

  double result = 0.0, meltarea = 0.0;

  const IceModelVec3 &enthalpy = m_energy_model->enthalpy();

  IceModelVec::AccessList list{&enthalpy, &m_geometry.cell_type, &m_geometry.ice_thickness};
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_geometry.cell_type.icy(i, j)) {
        const double
          E_basal  = enthalpy.get_column(i, j)[0],
          pressure = EC->pressure(m_geometry.ice_thickness(i,j)); // FIXME issue #15
        // accumulate area of base which is at melt point
        if (EC->is_temperate_relaxed(E_basal, pressure)) {
          meltarea += cell_area;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // convert from m2 to km2
  meltarea = units::convert(m_sys, meltarea, "m2", "km2");
  // communication
  result = GlobalSum(m_grid->com, meltarea);

  // normalize fraction correctly
  if (total_ice_area > 0.0) {
    result = result / total_ice_area;
  } else {
    result = 0.0;
  }
  return result;
}


/*!
  Computes fraction of the ice which is as old as the start of the run (original).
  Communication occurs here.
 */
double IceModel::compute_original_ice_fraction(double total_ice_volume) {

  double result = -1.0;  // result value if not age.enabled

  if (m_age_model == NULL) {
    return result;  // leave now
  }

  const double a = m_grid->cell_area() * 1e-3 * 1e-3, // area unit (km^2)
    currtime = m_time->current(); // seconds

  const IceModelVec3 &ice_age = m_age_model->age();

  IceModelVec::AccessList list{&m_geometry.cell_type, &m_geometry.ice_thickness, &ice_age};

  const double one_year = units::convert(m_sys, 1.0, "year", "seconds");
  double original_ice_volume = 0.0;

  // compute local original volume
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_geometry.cell_type.icy(i, j)) {
        // accumulate volume of ice which is original
        const double *age = ice_age.get_column(i, j);
        const int  ks = m_grid->kBelowHeight(m_geometry.ice_thickness(i,j));
        for (int k = 1; k <= ks; k++) {
          // ice in segment is original if it is as old as one year less than current time
          if (0.5 * (age[k - 1] + age[k]) > currtime - one_year) {
            original_ice_volume += a * 1.0e-3 * (m_grid->z(k) - m_grid->z(k - 1));
          }
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  // communicate to turn into global original fraction
  result = GlobalSum(m_grid->com, original_ice_volume);

  // normalize fraction correctly
  if (total_ice_volume > 0.0) {
    result = result / total_ice_volume;
  } else {
    result = 0.0;
  }
  return result;
}

//! Computes the ice volume, in m^3.
double IceModel::ice_volume(double thickness_threshold) const {
  IceModelVec::AccessList list{&m_geometry.ice_thickness};

  double volume = 0.0;

  auto cell_area = m_grid->cell_area();

  {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_geometry.ice_thickness(i,j) >= thickness_threshold) {
        volume += m_geometry.ice_thickness(i,j) * cell_area;
      }
    }
  }

  // Add the volume of the ice in Href:
  if (m_config->get_boolean("geometry.part_grid.enabled")) {
    list.add(m_geometry.ice_area_specific_volume);
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      volume += m_geometry.ice_area_specific_volume(i,j) * cell_area;
    }
  }

  return GlobalSum(m_grid->com, volume);
}

double IceModel::ice_volume_not_displacing_seawater(double thickness_threshold) const {
  const double
    sea_water_density = m_config->get_double("constants.sea_water.density"),
    ice_density       = m_config->get_double("constants.ice.density"),
    cell_area         = m_grid->cell_area();

  IceModelVec::AccessList list{&m_geometry.cell_type, &m_geometry.ice_thickness,
      &m_geometry.bed_elevation, &m_geometry.sea_level_elevation};

  double volume = 0.0;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      bed       = m_geometry.bed_elevation(i, j),
      thickness = m_geometry.ice_thickness(i, j),
      sea_level = m_geometry.sea_level_elevation(i, j);

    if (m_geometry.cell_type.grounded(i, j) and thickness > thickness_threshold) {
      const double cell_ice_volume = thickness * cell_area;
      if (bed > sea_level) {
        volume += cell_ice_volume;
      } else {
        const double max_floating_volume = (sea_level - bed) * cell_area * (sea_water_density / ice_density);
        volume += cell_ice_volume - max_floating_volume;
      }
    }
  } // end of the loop over grid points

  return GlobalSum(m_grid->com, volume);
}

//! Computes the sea level rise that would result if all the ice were melted.
double IceModel::sea_level_rise_potential(double thickness_threshold) const {
  const double
    water_density = m_config->get_double("constants.fresh_water.density"),
    ice_density   = m_config->get_double("constants.ice.density"),
    ocean_area    = m_config->get_double("constants.global_ocean_area");

  const double
    volume                  = ice_volume_not_displacing_seawater(thickness_threshold),
    additional_water_volume = (ice_density / water_density) * volume,
    sea_level_change        = additional_water_volume / ocean_area;

  return sea_level_change;
}

//! Computes the temperate ice volume, in m^3.
double IceModel::ice_volume_temperate(double thickness_threshold) const {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  const IceModelVec3 &ice_enthalpy = m_energy_model->enthalpy();

  auto cell_area = m_grid->cell_area();

  double volume = 0.0;

  IceModelVec::AccessList list{&m_geometry.ice_thickness, &ice_enthalpy};
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_geometry.ice_thickness(i,j) >= thickness_threshold) {
        const int ks = m_grid->kBelowHeight(m_geometry.ice_thickness(i,j));
        const double *Enth = ice_enthalpy.get_column(i,j);

        for (int k = 0; k < ks; ++k) {
          if (EC->is_temperate_relaxed(Enth[k],EC->pressure(m_geometry.ice_thickness(i,j)))) { // FIXME issue #15
            volume += (m_grid->z(k + 1) - m_grid->z(k)) * cell_area;
          }
        }

        if (EC->is_temperate_relaxed(Enth[ks],EC->pressure(m_geometry.ice_thickness(i,j)))) { // FIXME issue #15
          volume += (m_geometry.ice_thickness(i,j) - m_grid->z(ks)) * cell_area;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return GlobalSum(m_grid->com, volume);
}

//! Computes the cold ice volume, in m^3.
double IceModel::ice_volume_cold(double thickness_threshold) const {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  const IceModelVec3 &ice_enthalpy = m_energy_model->enthalpy();

  double volume = 0.0;

  auto cell_area = m_grid->cell_area();

  IceModelVec::AccessList list{&m_geometry.ice_thickness, &ice_enthalpy};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double thickness = m_geometry.ice_thickness(i, j);

      // count all ice, including cells which have so little they
      // are considered "ice-free"
      if (thickness >= thickness_threshold) {
        const int ks = m_grid->kBelowHeight(thickness);
        const double *Enth = ice_enthalpy.get_column(i, j);

        for (int k=0; k<ks; ++k) {
          if (not EC->is_temperate_relaxed(Enth[k], EC->pressure(thickness))) { // FIXME issue #15
            volume += (m_grid->z(k+1) - m_grid->z(k)) * cell_area;
          }
        }

        if (not EC->is_temperate_relaxed(Enth[ks], EC->pressure(thickness))) { // FIXME issue #15
          volume += (thickness - m_grid->z(ks)) * cell_area;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return GlobalSum(m_grid->com, volume);
}

//! Computes ice area, in m^2.
double IceModel::ice_area(double thickness_threshold) const {
  double area = 0.0;

  auto cell_area = m_grid->cell_area();

  IceModelVec::AccessList list{&m_geometry.ice_thickness};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_geometry.ice_thickness(i, j) >= thickness_threshold) {
      area += cell_area;
    }
  }

  return GlobalSum(m_grid->com, area);
}

//! Computes area of basal ice which is temperate, in m^2.
double IceModel::ice_area_temperate(double thickness_threshold) const {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  const IceModelVec3 &ice_enthalpy = m_energy_model->enthalpy();

  double area = 0.0;

  auto cell_area = m_grid->cell_area();

  IceModelVec::AccessList list{&m_geometry.ice_thickness, &ice_enthalpy};
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double
        thickness      = m_geometry.ice_thickness(i, j),
        basal_enthalpy = ice_enthalpy.get_column(i, j)[0];

      if (thickness >= thickness_threshold and
          EC->is_temperate_relaxed(basal_enthalpy, EC->pressure(thickness))) { // FIXME issue #15
        area += cell_area;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return GlobalSum(m_grid->com, area);
}

//! Computes area of basal ice which is cold, in m^2.
double IceModel::ice_area_cold(double thickness_threshold) const {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  const IceModelVec3 &ice_enthalpy = m_energy_model->enthalpy();

  double area = 0.0;

  auto cell_area = m_grid->cell_area();

  IceModelVec::AccessList list{&ice_enthalpy, &m_geometry.ice_thickness};
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double
        thickness = m_geometry.ice_thickness(i, j),
        basal_enthalpy = ice_enthalpy.get_column(i, j)[0];

      if (thickness >= thickness_threshold and
          not EC->is_temperate_relaxed(basal_enthalpy, EC->pressure(thickness))) { // FIXME issue #15
        area += cell_area;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return GlobalSum(m_grid->com, area);
}

//! Computes grounded ice area, in m^2.
double IceModel::ice_area_grounded(double thickness_threshold) const {
  double area = 0.0;

  auto cell_area = m_grid->cell_area();

  IceModelVec::AccessList list{&m_geometry.cell_type, &m_geometry.ice_thickness};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_geometry.cell_type.grounded(i, j) and
        m_geometry.ice_thickness(i, j) >= thickness_threshold) {
      area += cell_area;
    }
  }

  return GlobalSum(m_grid->com, area);
}

//! Computes floating ice area, in m^2.
double IceModel::ice_area_floating(double thickness_threshold) const {
  double area = 0.0;

  auto cell_area = m_grid->cell_area();

  IceModelVec::AccessList list{&m_geometry.cell_type, &m_geometry.ice_thickness};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_geometry.cell_type.ocean(i, j) and
        m_geometry.ice_thickness(i, j) >= thickness_threshold) {
      area += cell_area;
    }
  }

  return GlobalSum(m_grid->com, area);
}


} // end of namespace pism
