/* Copyright (C) 2017, 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 * This source code is included in diagnostics.cc and so it does not need any #include
 * directives here.
 */

namespace pism {
namespace diagnostics {

enum AmountKind {AMOUNT, MASS};

//! @brief Computes tendency_of_ice_amount, the ice amount rate of change.
class TendencyOfIceAmount : public Diag<IceModel>
{
public:
  TendencyOfIceAmount(const IceModel *m, AmountKind kind)
    : Diag<IceModel>(m),
    m_kind(kind),
    m_last_amount(m_grid, "last_ice_amount", WITHOUT_GHOSTS),
    m_interval_length(0.0) {

    std::string
      name           = "tendency_of_ice_amount",
      long_name      = "rate of change of the ice amount",
      standard_name  = "",
      internal_units = "kg m-2 second-1",
      external_units = "kg m-2 year-1";
    if (kind == MASS) {
      name           = "tendency_of_ice_mass";
      long_name      = "rate of change of the ice mass";
      internal_units = "kg second-1";
      external_units = "Gt year-1" ;
    }

    // set metadata:
    m_vars = {SpatialVariableMetadata(m_sys, name)};

    set_attrs(long_name, standard_name, internal_units, external_units, 0);

    units::Converter c(m_sys, external_units, internal_units);

    m_fill_value = c(m_fill_value);

    const double valid_range = c(1e6);

    m_vars[0].set_doubles("valid_range",  {-valid_range, valid_range});
    m_vars[0].set_double("_FillValue", m_fill_value);
    m_vars[0].set_string("cell_methods", "time: mean");

    m_last_amount.set_attrs("internal",
                            "ice amount at the time of the last report of " + name,
                            internal_units + " second", "");
  }
protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "", WITHOUT_GHOSTS));
    result->metadata() = m_vars[0];

    if (m_interval_length > 0.0) {
      double ice_density = m_config->get_double("constants.ice.density");

      auto cell_area = m_grid->cell_area();

      const IceModelVec2S& thickness = model->geometry().ice_thickness;
      const IceModelVec2S& area_specific_volume = model->geometry().ice_area_specific_volume;

      IceModelVec::AccessList list{result.get(),
          &thickness, &area_specific_volume, &m_last_amount};

      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        // m * (kg / m^3) = kg / m^2
        double amount = (thickness(i, j) + area_specific_volume(i, j)) * ice_density;

        (*result)(i, j) = (amount - m_last_amount(i, j)) / m_interval_length;

        if (m_kind == MASS) {
          // kg / m^2 * m^2 = kg
          (*result)(i, j) *= cell_area;
        }
      }
    } else {
      result->set(m_fill_value);
    }

    return result;
  }

  void reset_impl() {
    m_interval_length = 0.0;

    const IceModelVec2S& thickness = model->geometry().ice_thickness;
    const IceModelVec2S& area_specific_volume = model->geometry().ice_area_specific_volume;

    double ice_density = m_config->get_double("constants.ice.density");

    IceModelVec::AccessList list{&m_last_amount, &thickness, &area_specific_volume};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // m * (kg / m^3) = kg / m^2
      m_last_amount(i, j) = (thickness(i, j) + area_specific_volume(i, j)) * ice_density;
    }
  }

  void update_impl(double dt) {
    m_interval_length += dt;
  }

protected:
  AmountKind m_kind;
  IceModelVec2S m_last_amount;
  double m_interval_length;
};

//! @brief Computes tendency_of_ice_amount_due_to_flow, the rate of change of ice amount due to
//! flow.
/*! @brief Report rate of change of ice amount due to flow. */
class TendencyOfIceAmountDueToFlow : public DiagAverageRate<IceModel>
{
public:
  TendencyOfIceAmountDueToFlow(const IceModel *m, AmountKind kind)
    : DiagAverageRate<IceModel>(m,
                                kind == AMOUNT
                                ? "tendency_of_ice_amount_due_to_flow"
                                : "tendency_of_ice_mass_due_to_flow", TOTAL_CHANGE),
    m_kind(kind) {

    std::string
      name              = "tendency_of_ice_amount_due_to_flow",
      long_name         = "rate of change of ice amount due to flow",
      standard_name     = "",
      accumulator_units = "kg m-2",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";

    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_flow";
      long_name         = "rate of change of ice mass due to flow";
      standard_name     = "";
      accumulator_units = "kg";
      internal_units    = "kg second-1";
      external_units    = "Gt year-1";
    }

    m_factor = m_config->get_double("constants.ice.density");

    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value,
                                       external_units, internal_units);
    m_vars[0].set_double("_FillValue", fill_value);
    m_vars[0].set_string("comment", "positive flux corresponds to ice gain");
  }

protected:
  void update_impl(double dt) {
    const IceModelVec2S
      &dH = model->geometry_evolution().thickness_change_due_to_flow(),
      &dV = model->geometry_evolution().area_specific_volume_change_due_to_flow();

    auto cell_area = m_grid->cell_area();

    IceModelVec::AccessList list{&m_accumulator, &dH, &dV};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double C = m_factor * (m_kind == AMOUNT ? 1.0 : cell_area);

      m_accumulator(i, j) += C * (dH(i, j) + dV(i, j));
    }

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

/*! @brief Report surface mass balance flux, averaged over the reporting interval */
class SurfaceFlux : public DiagAverageRate<IceModel>
{
public:
  SurfaceFlux(const IceModel *m, AmountKind kind)
    : DiagAverageRate<IceModel>(m,
                                kind == AMOUNT
                                ? "tendency_of_ice_amount_due_to_surface_mass_flux"
                                : "tendency_of_ice_mass_due_to_surface_mass_flux",
                                TOTAL_CHANGE),
    m_kind(kind) {
    m_factor = m_config->get_double("constants.ice.density");

    std::string
      name              = "tendency_of_ice_amount_due_to_surface_mass_flux",
      accumulator_units = "kg m-2",
      long_name         = "average surface mass flux over reporting interval",
      standard_name     = "",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_surface_mass_flux",
      accumulator_units = "kg",
      long_name         = "average surface mass flux over reporting interval",
      standard_name     = "",
      internal_units    = "kg second-1",
      external_units    = "Gt year-1";
    }

    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);

    double fill_value = units::convert(m_sys, m_fill_value,
                                       external_units, internal_units);
    m_vars[0].set_double("_FillValue", fill_value);
    m_vars[0].set_string("cell_methods", "time: mean");
    m_vars[0].set_string("comment", "positive flux corresponds to ice gain");
  }

protected:
  void update_impl(double dt) {
    const IceModelVec2S
      &SMB = model->geometry_evolution().top_surface_mass_balance();

    auto cell_area = m_grid->cell_area();

    IceModelVec::AccessList list{&m_accumulator, &SMB};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double C = m_factor * (m_kind == AMOUNT ? 1.0 : cell_area);

      m_accumulator(i, j) += C * SMB(i, j);
    }

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

/*! @brief Report basal mass balance flux, averaged over the reporting interval */
class BasalFlux : public DiagAverageRate<IceModel>
{
public:
  BasalFlux(const IceModel *m, AmountKind kind)
    : DiagAverageRate<IceModel>(m,
                                kind == AMOUNT
                                ? "tendency_of_ice_amount_due_to_basal_mass_flux"
                                : "tendency_of_ice_mass_due_to_basal_mass_flux",
                                TOTAL_CHANGE),
    m_kind(kind) {
    m_factor = m_config->get_double("constants.ice.density");

    std::string
      name              = "tendency_of_ice_amount_due_to_basal_mass_flux",
      accumulator_units = "kg m-2",
      long_name         = "average basal mass flux over reporting interval",
      standard_name     = "",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_basal_mass_flux",
      accumulator_units = "kg",
      long_name         = "average basal mass flux over reporting interval",
      standard_name     = "tendency_of_land_ice_mass_due_to_basal_mass_balance",
      internal_units    = "kg second-1",
      external_units    = "Gt year-1";
    }
    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_double("_FillValue", fill_value);
    m_vars[0].set_string("cell_methods", "time: mean");
    m_vars[0].set_string("comment", "positive flux corresponds to ice gain");
  }

protected:
  void update_impl(double dt) {
    const IceModelVec2S
      &BMB = model->geometry_evolution().bottom_surface_mass_balance();

    auto cell_area = m_grid->cell_area();

    IceModelVec::AccessList list{&m_accumulator, &BMB};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double C = m_factor * (m_kind == AMOUNT ? 1.0 : cell_area);

      m_accumulator(i, j) += C * BMB(i, j);
    }

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

class ConservationErrorFlux : public DiagAverageRate<IceModel>
{
public:
  ConservationErrorFlux(const IceModel *m, AmountKind kind)
    : DiagAverageRate<IceModel>(m,
                                kind == AMOUNT
                                ? "tendency_of_ice_amount_due_to_conservation_error"
                                : "tendency_of_ice_mass_due_to_conservation_error" ,
                                TOTAL_CHANGE),
    m_kind(kind) {
    m_factor = m_config->get_double("constants.ice.density");

    std::string
      name              = "tendency_of_ice_amount_due_to_conservation_error",
      accumulator_units = "kg m-2",
      long_name         = "average mass conservation error flux over reporting interval",
      standard_name     = "",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_conservation_error",
      accumulator_units = "kg",
      long_name         = "average mass conservation error flux over reporting interval",
      standard_name     = "",
      internal_units    = "kg second-1",
      external_units    = "Gt year-1";
    }

    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_double("_FillValue", fill_value);
    m_vars[0].set_string("cell_methods", "time: mean");
    m_vars[0].set_string("comment", "positive flux corresponds to ice gain");
  }

protected:
  void update_impl(double dt) {
    const IceModelVec2S
      &error = model->geometry_evolution().conservation_error();

    IceModelVec::AccessList list{&m_accumulator, &error};

    auto cell_area = m_grid->cell_area();

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double C = m_factor * (m_kind == AMOUNT ? 1.0 : cell_area);

      m_accumulator(i, j) += C * error(i, j);
    }

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

/*! @brief Report discharge (calving and frontal melt) flux. */
class DischargeFlux : public DiagAverageRate<IceModel>
{
public:
  DischargeFlux(const IceModel *m, AmountKind kind)
    : DiagAverageRate<IceModel>(m,
                                kind == AMOUNT
                                ? "tendency_of_ice_amount_due_to_discharge"
                                : "tendency_of_ice_mass_due_to_discharge",
                                TOTAL_CHANGE),
    m_kind(kind) {

    m_factor = m_config->get_double("constants.ice.density");

    std::string
      name              = "tendency_of_ice_amount_due_to_discharge",
      long_name         = "discharge (calving and frontal melt) flux",
      accumulator_units = "kg m-2",
      standard_name     = "land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_discharge";
      long_name         = "discharge (calving and frontal melt) flux";
      accumulator_units = "kg";
      standard_name     = "";
      internal_units    = "kg second-1";
      external_units    = "Gt year-1";
    }

    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_double("_FillValue", fill_value);
    m_vars[0].set_string("comment", "positive flux corresponds to ice gain");
  }

protected:
  void update_impl(double dt) {
    const IceModelVec2S
      &discharge = model->discharge();

    IceModelVec::AccessList list{&m_accumulator, &discharge};

    auto cell_area = m_grid->cell_area();

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double C = m_factor * (m_kind == AMOUNT ? 1.0 : cell_area);

      m_accumulator(i, j) += C * discharge(i, j);
    }

    m_interval_length += dt;
  }
  AmountKind m_kind;
};


} // end of namespace diagnostics
} // end of namespace pism
