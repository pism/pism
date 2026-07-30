// Copyright (C) 2010--2026 Constantine Khroulev
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

#include <algorithm>
#include <cassert>
#include <memory>

#include "pism/age/AgeModel.hh"
#include "pism/energy/EnergyModel.hh"
#include "pism/energy/utilities.hh"
#include "pism/geometry/grounded_cell_fraction.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/icemodel/IceModel.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/stressbalance/SSB_Modifier.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/basalstrength/YieldStress.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/earth/BedDef.hh"
#include "pism/hydrology/Routing.hh"
#include "pism/energy/BedThermalUnit.hh"
#include "pism/stressbalance/sia/SIAFD.hh"
#include "pism/stressbalance/sia/BedSmoother.hh"

#if (Pism_USE_PROJ == 1)
#include "pism/util/Proj.hh"
#endif

// Flux balance code
namespace pism {
namespace diagnostics {

enum AmountKind { AMOUNT, MASS };

//! @brief Computes tendency_of_ice_amount, the ice amount rate of change.
class TendencyOfIceAmount : public Diag<IceModel> {
public:
  TendencyOfIceAmount(const IceModel *Model, AmountKind kind)
      : Diag<IceModel>(Model),
        m_kind(kind),
        m_last_amount(m_grid, "last_ice_amount"),
        m_interval_length(0.0) {

    std::string name = "tendency_of_ice_amount", long_name = "rate of change of the ice amount",
                internal_units = "kg m^-2 second^-1", external_units = "kg m^-2 year^-1";
    if (kind == MASS) {
      name           = "tendency_of_ice_mass";
      long_name      = "rate of change of the ice mass";
      internal_units = "kg second^-1";
      external_units = "Gt year^-1";
    }

    // set metadata:
    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0].long_name(long_name).units(internal_units).output_units(external_units);

    m_vars[0]["_FillValue"]   = { fill_value() };
    m_vars[0]["cell_methods"] = "time: mean";

    m_last_amount.metadata()
        .long_name("ice amount at the time of the last report of " + name)
        .units(internal_units + " second");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &/*geometry*/) const {

    auto result        = std::make_shared<array::Scalar>(m_grid, "");
    result->metadata() = m_vars[0];

    if (m_interval_length > 0.0) {
      double ice_density = m_config->get_number("constants.ice.density");

      auto cell_area = m_grid->cell_area();

      const auto &thickness            = model->geometry().ice_thickness;
      const auto &area_specific_volume = model->geometry().ice_area_specific_volume;

      array::AccessScope list{ result.get(), &thickness, &area_specific_volume, &m_last_amount };

      for (auto p : m_grid->points()) {
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
      result->set(fill_value());
    }

    return result;
  }

  void reset_impl() {
    m_interval_length = 0.0;

    const array::Scalar &thickness            = model->geometry().ice_thickness;
    const array::Scalar &area_specific_volume = model->geometry().ice_area_specific_volume;

    double ice_density = m_config->get_number("constants.ice.density");

    array::AccessScope list{ &m_last_amount, &thickness, &area_specific_volume };

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      // m * (kg / m^3) = kg / m^2
      m_last_amount(i, j) = (thickness(i, j) + area_specific_volume(i, j)) * ice_density;
    }
  }

  void update_impl(double dt) {
    m_interval_length += dt;
  }

  AmountKind m_kind;
  array::Scalar m_last_amount;
  double m_interval_length;
};

//! @brief Computes tendency_of_ice_amount_due_to_flow, the rate of change of ice amount due to
//! flow.
/*! @brief Report rate of change of ice amount due to flow. */
class TendencyOfIceAmountDueToFlow : public DiagAverageRate<IceModel> {
public:
  TendencyOfIceAmountDueToFlow(const IceModel *Model, AmountKind kind)
      : DiagAverageRate<IceModel>(Model,
                                  kind == AMOUNT ? "tendency_of_ice_amount_due_to_flow" :
                                                   "tendency_of_ice_mass_due_to_flow",
                                  TOTAL_CHANGE),
        m_kind(kind) {

    std::string name              = "tendency_of_ice_amount_due_to_flow",
                long_name         = "rate of change of ice amount due to flow",
                accumulator_units = "kg m^-2", internal_units = "kg m^-2 second^-1",
                external_units = "kg m^-2 year^-1";

    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_flow";
      long_name         = "rate of change of ice mass due to flow";
      accumulator_units = "kg";
      internal_units    = "kg second^-1";
      external_units    = "Gt year^-1";
    }

    m_factor = m_config->get_number("constants.ice.density");

    m_accumulator.metadata().units(accumulator_units);

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0].long_name(long_name).units(internal_units).output_units(external_units);
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["_FillValue"]   = { fill_value() };
    m_vars[0]["comment"]      = "positive flux corresponds to ice gain";
  }

protected:
  void update_impl(double dt) {
    const array::Scalar &dH = model->geometry_evolution().thickness_change_due_to_flow(),
                        &dV = model->geometry_evolution().area_specific_volume_change_due_to_flow();

    auto cell_area = m_grid->cell_area();

    array::AccessScope list{ &m_accumulator, &dH, &dV };

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      double C = m_factor * (m_kind == AMOUNT ? 1.0 : cell_area);

      m_accumulator(i, j) += C * (dH(i, j) + dV(i, j));
    }

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

/*! @brief Report surface mass balance flux, averaged over the reporting interval */
class SurfaceFlux : public DiagAverageRate<IceModel> {
public:
  SurfaceFlux(const IceModel *m, AmountKind kind)
      : DiagAverageRate<IceModel>(m,
                                  kind == AMOUNT ?
                                      "tendency_of_ice_amount_due_to_surface_mass_flux" :
                                      "tendency_of_ice_mass_due_to_surface_mass_flux",
                                  TOTAL_CHANGE),
        m_kind(kind) {
    m_factor = m_config->get_number("constants.ice.density");

    auto ismip = m_config->get_flag("output.ISMIP");

    std::string name = ismip ? "acabf" : "tendency_of_ice_amount_due_to_surface_mass_flux",
                accumulator_units = "kg m^-2",
                long_name         = "average surface mass flux over reporting interval",
                standard_name     = "land_ice_surface_specific_mass_balance_flux",
                internal_units = "kg m^-2 s^-1", external_units = "kg m^-2 year^-1";
    if (kind == MASS) {
      name = "tendency_of_ice_mass_due_to_surface_mass_flux", accumulator_units = "kg",
      long_name = "average surface mass flux over reporting interval", standard_name = "",
      internal_units = "kg second^-1", external_units = "Gt year^-1";
    }

    m_accumulator.metadata()["units"] = accumulator_units;

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0]
        .long_name(long_name)
        .standard_name(standard_name)
        .units(internal_units)
        .output_units(external_units);
    m_vars[0]["_FillValue"]   = { fill_value() };
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["comment"]      = "positive flux corresponds to ice gain";
  }

protected:
  void update_impl(double dt) {
    const array::Scalar &SMB = model->geometry_evolution().top_surface_mass_balance();

    auto cell_area = m_grid->cell_area();

    array::AccessScope list{ &m_accumulator, &SMB };

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      double C = m_factor * (m_kind == AMOUNT ? 1.0 : cell_area);

      m_accumulator(i, j) += C * SMB(i, j);
    }

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

/*! @brief Report basal mass balance flux, averaged over the reporting interval */
class BasalFlux : public DiagAverageRate<IceModel> {
public:
  BasalFlux(const IceModel *m, AmountKind kind)
      : DiagAverageRate<IceModel>(m,
                                  kind == AMOUNT ? "tendency_of_ice_amount_due_to_basal_mass_flux" :
                                                   "tendency_of_ice_mass_due_to_basal_mass_flux",
                                  TOTAL_CHANGE),
        m_kind(kind) {
    m_factor = m_config->get_number("constants.ice.density");

    std::string name              = "tendency_of_ice_amount_due_to_basal_mass_flux",
                accumulator_units = "kg m^-2",
                long_name      = "average basal mass flux over reporting interval", standard_name,
                internal_units = "kg m^-2 second^-1", external_units = "kg m^-2 year^-1";
    if (kind == MASS) {
      name = "tendency_of_ice_mass_due_to_basal_mass_flux", accumulator_units = "kg",
      long_name      = "average basal mass flux over reporting interval",
      standard_name  = "tendency_of_land_ice_mass_due_to_basal_mass_balance",
      internal_units = "kg second^-1", external_units = "Gt year^-1";
    }
    m_accumulator.metadata()["units"] = accumulator_units;

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0]
        .long_name(long_name)
        .standard_name(standard_name)
        .units(internal_units)
        .output_units(external_units);
    m_vars[0]["_FillValue"]   = { fill_value() };
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["comment"]      = "positive flux corresponds to ice gain";
  }

protected:
  void update_impl(double dt) {
    const array::Scalar &BMB = model->geometry_evolution().bottom_surface_mass_balance();

    auto cell_area = m_grid->cell_area();

    array::AccessScope list{ &m_accumulator, &BMB };

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      double C = m_factor * (m_kind == AMOUNT ? 1.0 : cell_area);

      m_accumulator(i, j) += C * BMB(i, j);
    }

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

class ConservationErrorFlux : public DiagAverageRate<IceModel> {
public:
  ConservationErrorFlux(const IceModel *m, AmountKind kind)
      : DiagAverageRate<IceModel>(m,
                                  kind == AMOUNT ?
                                      "tendency_of_ice_amount_due_to_conservation_error" :
                                      "tendency_of_ice_mass_due_to_conservation_error",
                                  TOTAL_CHANGE),
        m_kind(kind) {
    m_factor = m_config->get_number("constants.ice.density");

    std::string name              = "tendency_of_ice_amount_due_to_conservation_error",
                accumulator_units = "kg m^-2",
                long_name         = "average mass conservation error flux over reporting interval",
                internal_units = "kg m^-2 second^-1", external_units = "kg m^-2 year^-1";
    if (kind == MASS) {
      name = "tendency_of_ice_mass_due_to_conservation_error", accumulator_units = "kg",
      long_name     = "average mass conservation error flux over reporting interval",
      internal_units = "kg second^-1", external_units = "Gt year^-1";
    }

    m_accumulator.metadata()["units"] = accumulator_units;

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0]
        .long_name(long_name)
        .units(internal_units)
        .output_units(external_units);
    m_vars[0]["_FillValue"] = { fill_value() };
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["comment"] = "positive flux corresponds to ice gain";
  }

protected:
  void update_impl(double dt) {
    const array::Scalar
      &error = model->geometry_evolution().conservation_error();

    array::AccessScope list{&m_accumulator, &error};

    auto cell_area = m_grid->cell_area();

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      double C = m_factor * (m_kind == AMOUNT ? 1.0 : cell_area);

      m_accumulator(i, j) += C * error(i, j);
    }

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

enum ChangeKind {CALVING, FRONTAL_MELT, FORCED_RETREAT, TOTAL_DISCHARGE};

static void accumulate_changes(const IceModel *model, double factor, ChangeKind kind,
                               array::Scalar &accumulator) {

  const auto &calving        = model->calving();
  const auto &frontal_melt   = model->frontal_melt();
  const auto &forced_retreat = model->forced_retreat();

  auto grid = accumulator.grid();

  bool add_calving        = (kind == CALVING or kind == TOTAL_DISCHARGE);
  bool add_frontal_melt   = (kind == FRONTAL_MELT or kind == TOTAL_DISCHARGE);
  bool add_forced_retreat = (kind == FORCED_RETREAT or kind == TOTAL_DISCHARGE);

  array::AccessScope scope{ &accumulator };
  if (add_calving) {
    scope.add(calving);
  }
  if (add_frontal_melt) {
    scope.add(frontal_melt);
  }
  if (add_forced_retreat) {
    scope.add(forced_retreat);
  }

  for (auto p : grid->points()) {
    const int i = p.i(), j = p.j();

    if (add_calving) {
      accumulator(i, j) += factor * calving(i, j);
    }
    if (add_frontal_melt) {
      accumulator(i, j) += factor * frontal_melt(i, j);
    }
    if (add_forced_retreat) {
      accumulator(i, j) += factor * forced_retreat(i, j);
    }
  }
}


/*! @brief Report discharge (calving and frontal melt) flux. */
class DischargeFlux : public DiagAverageRate<IceModel> {
public:
  DischargeFlux(const IceModel *m, AmountKind kind)
      : DiagAverageRate<IceModel>(m,
                                  kind == AMOUNT ? "tendency_of_ice_amount_due_to_discharge" :
                                                   "tendency_of_ice_mass_due_to_discharge",
                                  TOTAL_CHANGE),
        m_kind(kind) {

    m_factor = m_config->get_number("constants.ice.density");

    auto ismip = m_config->get_flag("output.ISMIP");

    std::string name      = ismip ? "lifmassbf" : "tendency_of_ice_amount_due_to_discharge",
                long_name = "discharge flux (calving, frontal melt, forced retreat)",
                accumulator_units = "kg m^-2",
                standard_name  = "land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting",
                internal_units = "kg m^-2 s^-1", external_units = "kg m^-2 year^-1";
    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_discharge";
      long_name         = "discharge flux (calving, frontal melt, forced retreat)";
      accumulator_units = "kg";
      standard_name     = "";
      internal_units    = "kg second^-1";
      external_units    = "Gt year^-1";
    }

    m_accumulator.metadata()["units"] = accumulator_units;

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0]
        .long_name(long_name)
        .standard_name(standard_name)
        .units(internal_units)
        .output_units(external_units);
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["_FillValue"]   = { fill_value() };
    m_vars[0]["comment"]      = "positive flux corresponds to ice gain";
  }

protected:
  void update_impl(double dt) {
    accumulate_changes(model,
                       m_factor * (m_kind == AMOUNT ? 1.0 : m_grid->cell_area()),
                       TOTAL_DISCHARGE,
                       m_accumulator);

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

/*! @brief Report the calving flux. */
class CalvingFlux : public DiagAverageRate<IceModel>
{
public:
  CalvingFlux(const IceModel *m, AmountKind kind)
    : DiagAverageRate<IceModel>(m,
                                kind == AMOUNT
                                ? "tendency_of_ice_amount_due_to_calving"
                                : "tendency_of_ice_mass_due_to_calving",
                                TOTAL_CHANGE),
    m_kind(kind) {

    m_factor = m_config->get_number("constants.ice.density");

    auto ismip = m_config->get_flag("output.ISMIP");

    std::string
      name              = ismip ? "licalvf" : "tendency_of_ice_amount_due_to_calving",
      long_name         = "calving flux",
      accumulator_units = "kg m^-2",
      standard_name     = "land_ice_specific_mass_flux_due_to_calving",
      internal_units    = "kg m^-2 s^-1",
      external_units    = "kg m^-2 year^-1";
    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_calving";
      long_name         = "calving flux";
      accumulator_units = "kg";
      standard_name     = "";
      internal_units    = "kg second^-1";
      external_units    = "Gt year^-1";
    }

    m_accumulator.metadata().units(accumulator_units);

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0]
        .long_name(long_name)
        .standard_name(standard_name)
        .units(internal_units)
        .output_units(external_units);
    m_vars[0]["cell_methods"] = "time: mean";

    m_vars[0]["_FillValue"] = { fill_value() };
    m_vars[0]["comment"] = "positive flux corresponds to ice gain";
  }

protected:
  void update_impl(double dt) {

    accumulate_changes(model,
                       m_factor * (m_kind == AMOUNT ? 1.0 : m_grid->cell_area()),
                       CALVING,
                       m_accumulator);

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

/*! @brief Report the frontal melt flux. */
class FrontalMeltFlux : public DiagAverageRate<IceModel>
{
public:
  FrontalMeltFlux(const IceModel *m, AmountKind kind)
    : DiagAverageRate<IceModel>(m,
                                kind == AMOUNT
                                ? "tendency_of_ice_amount_due_to_frontal_melt"
                                : "tendency_of_ice_mass_due_to_frontal_melt",
                                TOTAL_CHANGE),
    m_kind(kind) {

    m_factor = m_config->get_number("constants.ice.density");

    std::string name = "tendency_of_ice_amount_due_to_frontal_melt", accumulator_units = "kg m^-2",
                internal_units = "kg m^-2 s^-1", external_units = "kg m^-2 year^-1";
    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_frontal_melt";
      accumulator_units = "kg";
      internal_units    = "kg second^-1";
      external_units    = "Gt year^-1";
    }

    m_accumulator.metadata().units(accumulator_units);

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0].long_name("frontal melt flux").units(internal_units).output_units(external_units);
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["_FillValue"] = { fill_value() };
    m_vars[0]["comment"] = "positive flux corresponds to ice gain";
  }

protected:
  void update_impl(double dt) {

    accumulate_changes(model,
                       m_factor * (m_kind == AMOUNT ? 1.0 : m_grid->cell_area()),
                       FRONTAL_MELT,
                       m_accumulator);

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

/*! @brief Report the frontal melt flux. */
class ForcedRetreatFlux : public DiagAverageRate<IceModel>
{
public:
  ForcedRetreatFlux(const IceModel *m, AmountKind kind)
      : DiagAverageRate<IceModel>(m,
                                  kind == AMOUNT ? "tendency_of_ice_amount_due_to_forced_retreat" :
                                                   "tendency_of_ice_mass_due_to_forced_retreat",
                                  TOTAL_CHANGE),
        m_kind(kind) {

    m_factor = m_config->get_number("constants.ice.density");

    std::string name = "tendency_of_ice_amount_due_to_forced_retreat", accumulator_units = "kg m^-2",
                internal_units = "kg m^-2 s^-1", external_units = "kg m^-2 year^-1";
    if (kind == MASS) {
      name              = "tendency_of_ice_mass_due_to_forced_retreat";
      accumulator_units = "kg";
      internal_units    = "kg second^-1";
      external_units    = "Gt year^-1";
    }

    m_accumulator.metadata().units(accumulator_units);

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0]
        .long_name("forced (prescribed) retreat flux")
        .units(internal_units)
        .output_units(external_units);
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["_FillValue"] = { fill_value() };
    m_vars[0]["comment"] = "positive flux corresponds to ice gain";
  }

protected:
  void update_impl(double dt) {

    accumulate_changes(model,
                       m_factor * (m_kind == AMOUNT ? 1.0 : m_grid->cell_area()),
                       FORCED_RETREAT,
                       m_accumulator);

    m_interval_length += dt;
  }
  AmountKind m_kind;
};

//! \brief Computes basal frictional heating.
class StressBalanceBfrict : public Diag<IceModel> {
public:
  StressBalanceBfrict(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

StressBalanceBfrict::StressBalanceBfrict(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "bfrict", *m_grid } };
  m_vars[0].long_name("basal frictional heating").units("W m^-2");
}

std::shared_ptr<array::Array>
StressBalanceBfrict::compute_impl(const Geometry & /*geometry*/) const {

  auto result = allocate<array::Scalar>("bfrict");

  result->copy_from(model->stress_balance()->basal_frictional_heating());

  return result;
}

//! \brief Reports the volumetric strain heating (3D).
class StressBalanceStrainheat : public Diag<IceModel> {
public:
  StressBalanceStrainheat(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

StressBalanceStrainheat::StressBalanceStrainheat(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "strainheat", *m_grid, m_grid->z() } };
  m_vars[0]
      .long_name("rate of strain heating in ice (dissipation heating)")
      .units("W m^-3")
      .output_units("mW m^-3");
}

std::shared_ptr<array::Array>
StressBalanceStrainheat::compute_impl(const Geometry & /*geometry*/) const {
  auto result =
      std::make_shared<array::Array3D>(m_grid, "strainheat", array::WITHOUT_GHOSTS, m_grid->z());

  result->metadata() = m_vars[0];

  result->copy_from(model->stress_balance()->volumetric_strain_heating());

  return result;
}

} // end of namespace diagnostics
} // end of namespace pism

namespace pism {

namespace details {
enum IceKind {ICE_COLD, ICE_TEMPERATE};

static double ice_volume(const array::Scalar &ice_thickness,
                         const array::Array3D &ice_enthalpy,
                         IceKind kind,
                         double thickness_threshold) {

  auto grid = ice_thickness.grid();
  auto ctx = grid->ctx();
  auto EC = ctx->enthalpy_converter();

  auto cell_area = grid->cell_area();
  const auto& z = grid->z();

  double volume = 0.0;

  // count the volume of a 3D grid cell if
  //
  // - it is temperate and we're asked for the temperate ice volume
  // - it is cold and we're asked for the cold ice volume
  //
  // return zero otherwise
  //
  // uses the depth at the *bottom* of a cell to compute pressure
  auto volume_counter = [EC, kind, cell_area](double z_min, double z_max, double H, double E) {
    double depth = H - z_min;
    double P = EC->pressure(depth);
    double V = cell_area * (z_max - z_min);
    bool temperate = EC->is_temperate_relaxed(E, P); // FIXME issue #15

    switch (kind) {
    case ICE_TEMPERATE:
      return temperate ? V : 0.0;
    default:
    case ICE_COLD:
      return (not temperate) ? V : 0.0;
    }
  };

  array::AccessScope list{&ice_thickness, &ice_enthalpy};
  ParallelSection loop(grid->com);
  try {
    for (auto p : grid->points()) {
      const int i = p.i(), j = p.j();

      double H = ice_thickness(i, j);

      if (H >= thickness_threshold) {
        auto ks = grid->kBelowHeight(H);
        const double *E = ice_enthalpy.get_column(i, j);

        for (unsigned int k = 0; k < ks; ++k) {
          volume += volume_counter(z[k], z[k + 1], H, E[k]);
        }

        volume += volume_counter(z[ks], H, H, E[ks]);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return GlobalSum(grid->com, volume);
}

static double base_area(const array::Scalar &ice_thickness,
                        const array::Array3D &ice_enthalpy,
                        IceKind kind,
                        double thickness_threshold) {

  auto grid = ice_thickness.grid();
  auto ctx = grid->ctx();
  auto EC = ctx->enthalpy_converter();

  auto cell_area = grid->cell_area();

  double area = 0.0;

  array::AccessScope list{&ice_thickness, &ice_enthalpy};
  ParallelSection loop(grid->com);
  try {
    for (auto p : grid->points()) {
      const int i = p.i(), j = p.j();

      double thickness = ice_thickness(i, j);

      if (thickness >= thickness_threshold) {
        double basal_enthalpy = ice_enthalpy.get_column(i, j)[0];

        bool temperate = EC->is_temperate_relaxed(basal_enthalpy,
                                                  EC->pressure(thickness)); // FIXME issue #15

        switch (kind) {
        case ICE_TEMPERATE:
          area += temperate ? cell_area : 0.0;
          break;
        default:
        case ICE_COLD:
          area += (not temperate) ? cell_area : 0.0;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return GlobalSum(grid->com, area);
}

} // end of namespace details

// Horrendous names used by InitMIP (and ISMIP, and CMIP5). Ugh.
static const char* land_ice_area_fraction_name           = "sftgif";
static const char* grounded_ice_sheet_area_fraction_name = "sftgrf";
static const char* floating_ice_sheet_area_fraction_name = "sftflf";

namespace diagnostics {

enum AreaType {GROUNDED, SHELF, BOTH};

enum TermType {SMB, BMB, FLOW, ERROR};

/*! @brief Ocean pressure difference at calving fronts. Used to debug CF boundary conditins. */
class IceMarginPressureDifference : public Diag<IceModel>
{
public:
  IceMarginPressureDifference(IceModel *m);

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

IceMarginPressureDifference::IceMarginPressureDifference(IceModel *m) : Diag<IceModel>(m) {

  m_vars = { { m_sys, "ice_margin_pressure_difference", *m_grid } };
  m_vars[0]
      .long_name(
          "vertically-averaged pressure difference at ice margins (including calving fronts)")
      .units("Pa");
  m_vars[0]["_FillValue"] = { fill_value() };
}

std::shared_ptr<array::Array> IceMarginPressureDifference::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Scalar>("ice_margin_pressure_difference");

  array::CellType1 mask(m_grid, "mask");

  const auto &H         = geometry.ice_thickness;
  const auto &bed       = geometry.bed_elevation;
  const auto &sea_level = geometry.sea_level_elevation;

  {
    const double H_threshold = m_config->get_number("stress_balance.ice_free_thickness_standard");
    GeometryCalculator gc(*m_config);
    gc.set_icefree_thickness(H_threshold);

    gc.compute_mask(sea_level, bed, H, mask);
  }

  const double rho_ice   = m_config->get_number("constants.ice.density"),
               rho_ocean = m_config->get_number("constants.sea_water.density"),
               g         = m_config->get_number("constants.standard_gravity");
  double fill = fill_value();

  array::AccessScope list{ &H, &bed, &mask, &sea_level, result.get() };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      double delta_p = 0.0;
      if (mask.grounded_ice(i, j) and grid::domain_edge(*m_grid, i, j)) {
        delta_p = 0.0;
      } else if (mask.icy(i, j) and mask.next_to_ice_free_ocean(i, j)) {
        double P_ice   = 0.5 * rho_ice * g * H(i, j),
               P_water = average_water_column_pressure(H(i, j), bed(i, j), sea_level(i, j), rho_ice,
                                                       rho_ocean, g);

        delta_p = P_ice - P_water;
      } else {
        delta_p = fill;
      }

      (*result)(i, j) = delta_p;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return result;
}

/*! @brief Report average basal mass balance flux over the reporting interval (grounded or floating
  areas) */
class BMBSplit : public DiagAverageRate<IceModel> {
public:
  BMBSplit(const IceModel *m, AreaType flag)
      : DiagAverageRate<IceModel>(
            m, flag == GROUNDED ? "basal_mass_flux_grounded" : "basal_mass_flux_floating",
            TOTAL_CHANGE),
        m_kind(flag) {
    assert(flag != BOTH);

    auto ismip = m_config->get_flag("output.ISMIP");

    std::string name, description, standard_name;
    if (m_kind == GROUNDED) {
      name          = ismip ? "libmassbfgr" : "basal_mass_flux_grounded";
      description   = "average basal mass flux over the reporting interval (grounded areas)";
      standard_name = ismip ? "land_ice_basal_specific_mass_balance_flux" : "";
    } else {
      name          = ismip ? "libmassbffl" : "basal_mass_flux_floating";
      description   = "average basal mass flux over the reporting interval (floating areas)";
      standard_name = ismip ? "land_ice_basal_specific_mass_balance_flux" : "";
    }

    m_accumulator.metadata()["units"] = "kg m^-2";

    m_vars = { { m_sys, name, *m_grid } };
    m_vars[0]
        .long_name(description)
        .standard_name(standard_name)
        .units("kg m^-2 s^-1")
        .output_units("kg m^-2 year^-1");
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["_FillValue"]   = { fill_value() };
    m_vars[0]["comment"]      = "positive flux corresponds to ice gain";
  }

protected:
  AreaType m_kind;
  void update_impl(double dt) {
    const array::Scalar &input = model->geometry_evolution().bottom_surface_mass_balance();
    const auto &cell_type      = model->geometry().cell_type;

    double ice_density = m_config->get_number("constants.ice.density");

    // the accumulator has the units of kg/m^2, computed as
    //
    // accumulator += BMB (m) * ice_density (kg / m^3)

    array::AccessScope list{ &input, &cell_type, &m_accumulator };

    for (auto p : m_grid->points()) {
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

//! \brief Computes vertically-averaged ice hardness.
class HardnessAverage : public Diag<IceModel> {
public:
  HardnessAverage(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

HardnessAverage::HardnessAverage(const IceModel *m) : Diag<IceModel>(m) {

  // set metadata:
  m_vars = { { m_sys, "hardav", *m_grid } };
  m_vars[0]
      .long_name("vertical average of ice hardness")
      .set_units_without_validation(
          "Pa s^(1/n)"); // n is the Glen exponent used by the SSA (shallow stress balance) flow law
  m_vars[0]["valid_min"]  = { 0.0 };
  m_vars[0]["_FillValue"] = { fill_value() };
  m_vars[0]["comment"]    = "units depend on the Glen exponent used by the flow law";
}

//! \brief Computes vertically-averaged ice hardness.
std::shared_ptr<array::Array> HardnessAverage::compute_impl(const Geometry &geometry) const {

  const rheology::FlowLaw *flow_law = model->stress_balance()->shallow()->flow_law().get();
  if (flow_law == NULL) {
    flow_law = model->stress_balance()->modifier()->flow_law().get();
    if (flow_law == NULL) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "Can't compute vertically-averaged hardness: no flow law is used.");
    }
  }

  auto result        = std::make_shared<array::Scalar>(m_grid, "hardav");
  result->metadata() = m_vars[0];

  const auto &cell_type = geometry.cell_type;

  const array::Array3D &ice_enthalpy = model->energy_balance_model()->enthalpy();
  const array::Scalar &ice_thickness = geometry.ice_thickness;
  double fill = fill_value();

  array::AccessScope list{ &cell_type, &ice_enthalpy, &ice_thickness, result.get() };
  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      const double *Eij = ice_enthalpy.get_column(i, j);
      const double H    = ice_thickness(i, j);
      if (cell_type.icy(i, j)) {
        (*result)(i, j) = rheology::averaged_hardness(*flow_law, H, m_grid->kBelowHeight(H),
                                                      m_grid->z().data(), Eij);
      } else { // put negative value below valid range
        (*result)(i, j) = fill;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

//! \brief Computes a diagnostic field filled with processor rank values.
class Rank : public Diag<IceModel> {
public:
  Rank(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

Rank::Rank(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "rank", *m_grid } };
  m_vars[0]
      .long_name("processor rank")
      .units("1")
      .set_time_dependent(false)
      .set_output_type(io::PISM_INT);
}

std::shared_ptr<array::Array> Rank::compute_impl(const Geometry &/*geometry*/) const {

  auto result        = std::make_shared<array::Scalar>(m_grid, "rank");
  result->metadata() = m_vars[0];

  array::AccessScope list{ result.get() };

  for (auto p : m_grid->points()) {
    (*result)(p.i(), p.j()) = m_grid->rank();
  }

  return result;
}

//! \brief Computes CTS, CTS = E/E_s(p).
class CTS : public Diag<IceModel> {
public:
  CTS(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

CTS::CTS(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "cts", *m_grid, m_grid->z() } };
  m_vars[0]
      .long_name("cts = E/E_s(p), so cold-temperate transition surface is at cts = 1")
      .units("1");
}

std::shared_ptr<array::Array> CTS::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "cts", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata() = m_vars[0];

  energy::compute_cts(model->energy_balance_model()->enthalpy(), geometry.ice_thickness, *result);

  return result;
}

//! \brief Computes ice temperature from enthalpy.
class Temperature : public Diag<IceModel> {
public:
  Temperature(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

Temperature::Temperature(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "temp", *m_grid, m_grid->z() } };
  m_vars[0]
      .long_name("ice temperature")
      .standard_name("land_ice_temperature")
      .units("kelvin");
  m_vars[0]["valid_min"] = { 0.0 };
}

std::shared_ptr<array::Array> Temperature::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "temp", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata() = m_vars[0];

  const auto &thickness = geometry.ice_thickness;
  const auto &enthalpy  = model->energy_balance_model()->enthalpy();

  auto EC = model->ctx()->enthalpy_converter();

  double *Tij;
  const double *Enthij; // columns of these values

  array::AccessScope list{result.get(), &enthalpy, &thickness};

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
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

//! \brief Compute the pressure-adjusted temperature in degrees C corresponding
//! to ice temperature.
class TemperaturePA : public Diag<IceModel>
{
public:
  TemperaturePA(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};


TemperaturePA::TemperaturePA(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars = { { m_sys, "temp_pa", *m_grid, m_grid->z() } };
  m_vars[0]
      .long_name("pressure-adjusted ice temperature (degrees above pressure-melting point)")
      .units("deg_C");
  m_vars[0]["valid_max"] = {0};
}

std::shared_ptr<array::Array> TemperaturePA::compute_impl(const Geometry &geometry) const {
  bool cold_mode = set_member(m_config->get_string("energy.model"), {"cold", "none"});
  double melting_point_temp = m_config->get_number("constants.fresh_water.melting_point_temperature");

  auto result = std::make_shared<array::Array3D>(m_grid, "temp_pa", array::WITHOUT_GHOSTS, m_grid->z());
  result->metadata() = m_vars[0];

  const auto &thickness = geometry.ice_thickness;
  const auto &enthalpy  = model->energy_balance_model()->enthalpy();

  auto EC = model->ctx()->enthalpy_converter();

  double *Tij;
  const double *Enthij; // columns of these values

  array::AccessScope list{result.get(), &enthalpy, &thickness};

  ParallelSection loop(m_grid->com);
  try {
    for (auto pt : m_grid->points()) {
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

//! \brief Computes basal values of the pressure-adjusted temperature.
class TemperaturePABasal : public Diag<IceModel>
{
public:
  TemperaturePABasal(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

TemperaturePABasal::TemperaturePABasal(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars = { { m_sys, "temppabase", *m_grid } };
  m_vars[0].long_name("pressure-adjusted ice temperature at the base of ice").units("degree_Celsius");
}

std::shared_ptr<array::Array> TemperaturePABasal::compute_impl(const Geometry &geometry) const {

  bool cold_mode = set_member(m_config->get_string("energy.model"), {"cold", "none"});
  double melting_point_temp = m_config->get_number("constants.fresh_water.melting_point_temperature");

  auto result = std::make_shared<array::Scalar>(m_grid, "temp_pa_base");
  result->metadata() = m_vars[0];

  const auto &thickness = geometry.ice_thickness;
  const auto &enthalpy = model->energy_balance_model()->enthalpy();

  auto EC = model->ctx()->enthalpy_converter();

  array::AccessScope list{result.get(), &enthalpy, &thickness};

  ParallelSection loop(m_grid->com);
  try {
    for (auto pt : m_grid->points()) {
      const int i = pt.i(), j = pt.j();

      const auto *Enthij = enthalpy.get_column(i,j);

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

//! \brief Computes surface values of ice enthalpy.
class IceEnthalpySurface : public Diag<IceModel>
{
public:
  IceEnthalpySurface(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

IceEnthalpySurface::IceEnthalpySurface(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars = { { m_sys, "enthalpysurf", *m_grid } };
  m_vars[0].long_name("ice enthalpy at 1m below the ice surface").units("J kg^-1");
  m_vars[0]["_FillValue"] = {fill_value()};
}

std::shared_ptr<array::Array> IceEnthalpySurface::compute_impl(const Geometry &geometry) const {

  auto result = std::make_shared<array::Scalar>(m_grid, "enthalpysurf");
  result->metadata() = m_vars[0];

  // compute levels corresponding to 1 m below the ice surface:

  const array::Array3D& ice_enthalpy = model->energy_balance_model()->enthalpy();
  const array::Scalar& ice_thickness = geometry.ice_thickness;

  array::AccessScope list{&ice_thickness, result.get()};

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) = std::max(ice_thickness(i,j) - 1.0, 0.0);
  }

  extract_surface(ice_enthalpy, *result, *result);  // slice at 1 m below the surface

  double fill = fill_value();

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (ice_thickness(i,j) <= 1.0) {
      (*result)(i,j) = fill;
    }
  }

  return result;
}

//! \brief Computes enthalpy at the base of the ice.
class IceEnthalpyBasal : public Diag<IceModel>
{
public:
  IceEnthalpyBasal(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

IceEnthalpyBasal::IceEnthalpyBasal(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars = { { m_sys, "enthalpybase", *m_grid } };
  m_vars[0].long_name("ice enthalpy at the base of ice").units("J kg^-1");
  m_vars[0]["_FillValue"] = {fill_value()};
}

std::shared_ptr<array::Array> IceEnthalpyBasal::compute_impl(const Geometry &geometry) const {

  auto result = std::make_shared<array::Scalar>(m_grid, "enthalpybase");
  result->metadata() = m_vars[0];

  extract_surface(model->energy_balance_model()->enthalpy(), 0.0, *result);  // z=0 slice

  apply_mask(geometry.ice_thickness, fill_value(), *result);

  return result;
}

//! \brief Computes ice temperature at the base of the ice.
class TemperatureBasal : public Diag<IceModel>
{
public:
  TemperatureBasal(const IceModel *m, AreaType area_type);
private:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;

  AreaType m_area_type;
};

TemperatureBasal::TemperatureBasal(const IceModel *m, AreaType area_type)
  : Diag<IceModel>(m), m_area_type(area_type) {

  std::string name, long_name, standard_name;
  switch (area_type) {
  case GROUNDED:
    name          = "litempbotgr";
    long_name     = "ice temperature at the bottom surface of grounded ice";
    standard_name = "temperature_at_base_of_ice_sheet_model";
    break;
  case SHELF:
    name          = "litempbotfl";
    long_name     = "ice temperature at the bottom surface of floating ice";
    standard_name = "temperature_at_base_of_ice_sheet_model";
    break;
  case BOTH:
    name          = "tempbase";
    long_name     = "ice temperature at the base of ice";
    standard_name = "land_ice_basal_temperature";
    break;
  }

  m_vars = { { m_sys, name, *m_grid } };
  m_vars[0].long_name(long_name).standard_name(standard_name).units("kelvin");
  m_vars[0]["_FillValue"] = { fill_value() };
}

std::shared_ptr<array::Array> TemperatureBasal::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Scalar>("basal_temperature");

  const auto &thickness = geometry.ice_thickness;

  auto EC = model->ctx()->enthalpy_converter();

  extract_surface(model->energy_balance_model()->enthalpy(), 0.0, *result); // z=0 (basal) slice
  // Now result contains basal enthalpy.

  const auto &cell_type = geometry.cell_type;
  double fill = fill_value();

  array::AccessScope list{ &cell_type, result.get(), &thickness };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      double depth = thickness(i, j), pressure = EC->pressure(depth),
             T = EC->temperature((*result)(i, j), pressure);

      if ((m_area_type == BOTH and cell_type.icy(i, j)) or
          (m_area_type == GROUNDED and cell_type.grounded_ice(i, j)) or
          (m_area_type == SHELF and cell_type.floating_ice(i, j))) {
        (*result)(i, j) = T;
      } else {
        (*result)(i, j) = fill;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

//! \brief Computes ice temperature at the surface of the ice.
class TemperatureSurface : public Diag<IceModel> {
public:
  TemperatureSurface(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

TemperatureSurface::TemperatureSurface(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "tempsurf", *m_grid } };
  m_vars[0]
      .long_name("ice temperature at 1m below the ice surface")
      .standard_name("temperature_at_ground_level_in_snow_or_firn") // InitMIP "standard" name
      .units("kelvin");
  m_vars[0]["_FillValue"] = { fill_value() };
}

std::shared_ptr<array::Array> TemperatureSurface::compute_impl(const Geometry &geometry) const {

  const array::Scalar &thickness = geometry.ice_thickness;

  auto enth   = IceEnthalpySurface(model).compute(geometry);
  auto result = array::cast<array::Scalar>(enth);

  auto EC = model->ctx()->enthalpy_converter();
  double fill = fill_value();

  // result contains surface enthalpy; note that it is allocated by
  // IceEnthalpySurface::compute().

  array::AccessScope list{ result.get(), &thickness };

  double depth = 1.0, pressure = EC->pressure(depth);
  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      if (thickness(i, j) > 1) {
        (*result)(i, j) = EC->temperature((*result)(i, j), pressure);
      } else {
        (*result)(i, j) = fill;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->metadata(0) = m_vars[0];
  return result;
}

//! \brief Computes the liquid water fraction.
class LiquidFraction : public Diag<IceModel> {
public:
  LiquidFraction(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

LiquidFraction::LiquidFraction(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "liqfrac", *m_grid,  m_grid->z() } };
  m_vars[0].long_name("liquid water fraction in ice (between 0 and 1)").units("1");
  m_vars[0]["valid_range"] = { 0.0, 1.0 };
}

std::shared_ptr<array::Array> LiquidFraction::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "liqfrac", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata(0) = m_vars[0];

  bool cold_mode = set_member(m_config->get_string("energy.model"), {"cold", "none"});

  if (cold_mode) {
    result->set(0.0);
  } else {
    energy::compute_liquid_water_fraction(model->energy_balance_model()->enthalpy(),
                                          geometry.ice_thickness, *result);
  }

  return result;
}

//! \brief Computes the total thickness of temperate ice in a column.
class TemperateIceThickness : public Diag<IceModel> {
public:
  TemperateIceThickness(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

TemperateIceThickness::TemperateIceThickness(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "tempicethk", *m_grid } };
  m_vars[0].long_name("temperate ice thickness (total column content)").units("m");
  m_vars[0]["_FillValue"] = { fill_value() };
}

std::shared_ptr<array::Array> TemperateIceThickness::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Scalar>("tempicethk");

  const auto &cell_type     = geometry.cell_type;
  const auto &ice_enthalpy  = model->energy_balance_model()->enthalpy();
  const auto &ice_thickness = geometry.ice_thickness;

  array::AccessScope list{ &cell_type, result.get(), &ice_enthalpy, &ice_thickness };

  auto EC = model->ctx()->enthalpy_converter();
  double fill = fill_value();

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      if (cell_type.icy(i, j)) {
        const double *Enth    = ice_enthalpy.get_column(i, j);
        double H_temperate    = 0.0;
        const double H        = ice_thickness(i, j);
        const unsigned int ks = m_grid->kBelowHeight(H);

        for (unsigned int k = 0; k < ks; ++k) { // FIXME issue #15
          double pressure = EC->pressure(H - m_grid->z(k));

          if (EC->is_temperate_relaxed(Enth[k], pressure)) {
            H_temperate += m_grid->z(k + 1) - m_grid->z(k);
          }
        }

        double pressure = EC->pressure(H - m_grid->z(ks));
        if (EC->is_temperate_relaxed(Enth[ks], pressure)) {
          H_temperate += H - m_grid->z(ks);
        }

        (*result)(i, j) = H_temperate;
      } else {
        // ice-free
        (*result)(i, j) = fill;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

//! \brief Computes the thickness of the basal layer of temperate ice.
class TemperateIceThicknessBasal : public Diag<IceModel> {
public:
  TemperateIceThicknessBasal(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

TemperateIceThicknessBasal::TemperateIceThicknessBasal(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "tempicethk_basal", *m_grid } };
  m_vars[0].long_name("thickness of the basal layer of temperate ice").units("m");
  m_vars[0]["_FillValue"] = { fill_value() };
}

/*!
 * Uses linear interpolation to go beyond vertical grid resolution.
 */
std::shared_ptr<array::Array> TemperateIceThicknessBasal::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Scalar>("tempicethk_basal");

  auto EC = model->ctx()->enthalpy_converter();

  const auto &cell_type     = geometry.cell_type;
  const auto &ice_enthalpy  = model->energy_balance_model()->enthalpy();
  const auto &ice_thickness = geometry.ice_thickness;
  double fill = fill_value();

  array::AccessScope list{ &cell_type, result.get(), &ice_thickness, &ice_enthalpy };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      double H = ice_thickness(i, j);

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (cell_type.ice_free(i, j)) {
        (*result)(i, j) = fill;
        continue;
      }

      const double *Enth = ice_enthalpy.get_column(i, j);

      unsigned int ks = m_grid->kBelowHeight(H);

      unsigned int k  = 0;
      double pressure = EC->pressure(H - m_grid->z(k));
      while (k <= ks) { // FIXME issue #15
        pressure = EC->pressure(H - m_grid->z(k));

        if (EC->is_temperate_relaxed(Enth[k], pressure)) {
          k++;
        } else {
          break;
        }
      }
      // after this loop 'pressure' is equal to the pressure at the first level
      // that is cold

      // no temperate ice at all; go to the next grid point
      if (k == 0) {
        (*result)(i, j) = 0.0;
        continue;
      }

      // the whole column is temperate (except, possibly, some ice between
      // z(ks) and the total thickness; we ignore it)
      if (k == ks + 1) {
        (*result)(i, j) = m_grid->z(ks);
        continue;
      }

      double pressure_0 = EC->pressure(H - m_grid->z(k - 1)), dz = m_grid->z(k) - m_grid->z(k - 1),
             slope1 = (Enth[k] - Enth[k - 1]) / dz,
             slope2 = (EC->enthalpy_cts(pressure) - EC->enthalpy_cts(pressure_0)) / dz;

      if (slope1 != slope2) {
        (*result)(i, j) =
            m_grid->z(k - 1) + (EC->enthalpy_cts(pressure_0) - Enth[k - 1]) / (slope1 - slope2);

        // check if the resulting thickness is valid:
        (*result)(i, j) = std::max((*result)(i, j), m_grid->z(k - 1));
        (*result)(i, j) = std::min((*result)(i, j), m_grid->z(k));
      } else {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Linear interpolation of the thickness of"
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
class IceVolumeGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceVolumeGlacierized(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized") {

    set_units("m^3", "m^3");
    m_variable["long_name"] = "volume of the ice in glacierized areas";
    m_variable["valid_min"] = { 0.0 };
  }
  double compute() {
    return ice_volume(model->geometry(),
                      m_config->get_number("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total ice volume.
class IceVolume : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceVolume(IceModel *m) : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume") {

    set_units("m^3", "m^3");
    m_variable["long_name"] = "volume of the ice, including seasonal cover";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    return ice_volume(model->geometry(), 0.0);
  }
};

//! \brief Computes the total ice volume which is relevant for sea-level
class SeaLevelRisePotential : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  SeaLevelRisePotential(const IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "sea_level_rise_potential") {

    set_units("m", "m");
    m_variable["long_name"] = "the sea level rise that would result if all the ice were melted";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    return sea_level_rise_potential(model->geometry(),
                                    m_config->get_number("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the rate of change of the total ice volume in glacierized areas.
class IceVolumeRateOfChangeGlacierized : public TSDiag<TSRateDiagnostic, IceModel> {
public:
  IceVolumeRateOfChangeGlacierized(IceModel *m)
      : TSDiag<TSRateDiagnostic, IceModel>(m, "tendency_of_ice_volume_glacierized") {

    set_units("m^3 s^-1", "m^3 year^-1");
    m_variable["long_name"] = "rate of change of the ice volume in glacierized areas";
  }

  double compute() {
    return ice_volume(model->geometry(),
                      m_config->get_number("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the rate of change of the total ice volume.
class IceVolumeRateOfChange : public TSDiag<TSRateDiagnostic, IceModel> {
public:
  IceVolumeRateOfChange(IceModel *m)
      : TSDiag<TSRateDiagnostic, IceModel>(m, "tendency_of_ice_volume") {

    set_units("m^3 s^-1", "m^3 year^-1");
    m_variable["long_name"] = "rate of change of the ice volume, including seasonal cover";
  }

  double compute() {
    return ice_volume(model->geometry(), 0.0);
  }
};

//! \brief Computes the total ice area.
class IceAreaGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceAreaGlacierized(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized") {

    set_units("m^2", "m^2");
    m_variable["long_name"] = "glacierized area";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    return ice_area(model->geometry(), m_config->get_number("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total mass of the ice not displacing sea water.
class IceMassNotDisplacingSeaWater : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceMassNotDisplacingSeaWater(const IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "limnsw") {

    set_units("kg", "kg");
    m_variable["long_name"]     = "mass of the ice not displacing sea water";
    m_variable["standard_name"] = "land_ice_mass_not_displacing_sea_water";
    m_variable["valid_min"]     = { 0.0 };
  }

  double compute() {

    const double thickness_standard = m_config->get_number("output.ice_free_thickness_standard"),
                 ice_density        = m_config->get_number("constants.ice.density"),
                 ice_volume =
                     ice_volume_not_displacing_seawater(model->geometry(), thickness_standard),
                 ice_mass = ice_volume * ice_density;

    return ice_mass;
  }
};

//! \brief Computes the total ice mass in glacierized areas.
class IceMassGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceMassGlacierized(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_mass_glacierized") {

    set_units("kg", "kg");
    m_variable["long_name"] = "mass of the ice in glacierized areas";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    double ice_density        = m_config->get_number("constants.ice.density"),
           thickness_standard = m_config->get_number("output.ice_free_thickness_standard");
    return ice_volume(model->geometry(), thickness_standard) * ice_density;
  }
};

//! \brief Computes the total ice mass.
class IceMass : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceMass(IceModel *m) : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_mass") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("lim");
    }

    set_units("kg", "kg");
    m_variable["long_name"]     = "mass of the ice, including seasonal cover";
    m_variable["standard_name"] = "land_ice_mass";
    m_variable["valid_min"]     = { 0.0 };
  }

  double compute() {
    return (ice_volume(model->geometry(), 0.0) * m_config->get_number("constants.ice.density"));
  }
};

//! \brief Computes the rate of change of the total ice mass in glacierized areas.
class IceMassRateOfChangeGlacierized : public TSDiag<TSRateDiagnostic, IceModel> {
public:
  IceMassRateOfChangeGlacierized(IceModel *m)
      : TSDiag<TSRateDiagnostic, IceModel>(m, "tendency_of_ice_mass_glacierized") {

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"] = "rate of change of the ice mass in glacierized areas";
  }

  double compute() {
    double ice_density         = m_config->get_number("constants.ice.density"),
           thickness_threshold = m_config->get_number("output.ice_free_thickness_standard");
    return ice_volume(model->geometry(), thickness_threshold) * ice_density;
  }
};

//! \brief Computes the rate of change of the total ice mass due to flow (influx due to
//! prescribed constant-in-time ice thickness).
/*!
 * This is the change in mass resulting from prescribing (fixing) ice thickness.
 */
class IceMassRateOfChangeDueToFlow : public TSDiag<TSFluxDiagnostic, IceModel> {
public:
  IceMassRateOfChangeDueToFlow(IceModel *m)
      : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_flow") {

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"] = "rate of change of the mass of ice due to flow"
                              " (i.e. prescribed ice thickness)";
  }

  double compute() {

    const double ice_density = m_config->get_number("constants.ice.density");

    const array::Scalar &dH = model->geometry_evolution().thickness_change_due_to_flow(),
                        &dV = model->geometry_evolution().area_specific_volume_change_due_to_flow();

    auto cell_area = m_grid->cell_area();

    array::AccessScope list{ &dH, &dV };

    double volume_change = 0.0;
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      // m * m^2 = m^3
      volume_change += (dH(i, j) + dV(i, j)) * cell_area;
    }

    // (kg/m^3) * m^3 = kg
    return ice_density * GlobalSum(m_grid->com, volume_change);
  }
};

//! \brief Computes the rate of change of the total ice mass.
class IceMassRateOfChange : public TSDiag<TSRateDiagnostic, IceModel> {
public:
  IceMassRateOfChange(IceModel *m) : TSDiag<TSRateDiagnostic, IceModel>(m, "tendency_of_ice_mass") {

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"] = "rate of change of the mass of ice, including seasonal cover";
  }

  double compute() {
    const double ice_density = m_config->get_number("constants.ice.density");
    return ice_volume(model->geometry(), 0.0) * ice_density;
  }
};


//! \brief Computes the total volume of the temperate ice in glacierized areas.
class IceVolumeGlacierizedTemperate : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceVolumeGlacierizedTemperate(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized_temperate") {

    set_units("m^3", "m^3");
    m_variable["long_name"] = "volume of temperate ice in glacierized areas";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    auto thickness_threshold = m_config->get_number("output.ice_free_thickness_standard");
    return details::ice_volume(model->geometry().ice_thickness,
                               model->energy_balance_model()->enthalpy(), details::ICE_TEMPERATE,
                               thickness_threshold);
  }
};

//! \brief Computes the total volume of the temperate ice.
class IceVolumeTemperate : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceVolumeTemperate(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_temperate") {

    set_units("m^3", "m^3");
    m_variable["long_name"] = "volume of temperate ice, including seasonal cover";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    return details::ice_volume(model->geometry().ice_thickness,
                               model->energy_balance_model()->enthalpy(), details::ICE_TEMPERATE,
                               0.0);
  }
};

//! \brief Computes the total volume of the cold ice in glacierized areas.
class IceVolumeGlacierizedCold : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceVolumeGlacierizedCold(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized_cold") {

    set_units("m^3", "m^3");
    m_variable["long_name"] = "volume of cold ice in glacierized areas";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    auto thickness_threshold = m_config->get_number("output.ice_free_thickness_standard");
    return details::ice_volume(model->geometry().ice_thickness,
                               model->energy_balance_model()->enthalpy(), details::ICE_COLD,
                               thickness_threshold);
  }
};

//! \brief Computes the total volume of the cold ice.
class IceVolumeCold : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceVolumeCold(IceModel *m) : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_cold") {

    set_units("m^3", "m^3");
    m_variable["long_name"] = "volume of cold ice, including seasonal cover";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    return details::ice_volume(model->geometry().ice_thickness,
                               model->energy_balance_model()->enthalpy(), details::ICE_COLD,
                               0.0);
  }
};

//! \brief Computes the total area of the temperate ice.
class IceAreaGlacierizedTemperateBase : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceAreaGlacierizedTemperateBase(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized_temperate_base") {

    set_units("m^2", "m^2");
    m_variable["long_name"] = "glacierized area where basal ice is temperate";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    auto thickness_threshold = m_config->get_number("output.ice_free_thickness_standard");
    return details::base_area(model->geometry().ice_thickness,
                              model->energy_balance_model()->enthalpy(), details::ICE_TEMPERATE,
                              thickness_threshold);
  }
};

//! \brief Computes the total area of the cold ice.
class IceAreaGlacierizedColdBase : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceAreaGlacierizedColdBase(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized_cold_base") {

    set_units("m^2", "m^2");
    m_variable["long_name"] = "glacierized area where basal ice is cold";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    auto thickness_threshold = m_config->get_number("output.ice_free_thickness_standard");
    return details::base_area(model->geometry().ice_thickness,
                              model->energy_balance_model()->enthalpy(), details::ICE_COLD,
                              thickness_threshold);
  }
};

//! \brief Computes the total ice enthalpy in glacierized areas.
class IceEnthalpyGlacierized : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceEnthalpyGlacierized(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_enthalpy_glacierized") {

    set_units("J", "J");
    m_variable["long_name"] = "enthalpy of the ice in glacierized areas";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    return energy::total_ice_enthalpy(m_config->get_number("output.ice_free_thickness_standard"),
                                      model->energy_balance_model()->enthalpy(),
                                      model->geometry().ice_thickness);
  }
};

//! \brief Computes the total ice enthalpy.
class IceEnthalpy : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceEnthalpy(IceModel *m) : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_enthalpy") {

    set_units("J", "J");
    m_variable["long_name"] = "enthalpy of the ice, including seasonal cover";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    return energy::total_ice_enthalpy(0.0, model->energy_balance_model()->enthalpy(),
                                      model->geometry().ice_thickness);
  }
};

//! \brief Computes the total grounded ice area.
class IceAreaGlacierizedGrounded : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceAreaGlacierizedGrounded(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized_grounded") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("iareagr");
    }

    set_units("m^2", "m^2");
    m_variable["long_name"]     = "area of grounded ice in glacierized areas";
    m_variable["standard_name"] = "grounded_ice_sheet_area";
    m_variable["valid_min"]     = { 0.0 };
  }

  double compute() {
    return ice_area_grounded(model->geometry(),
                             m_config->get_number("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total floating ice area.
class IceAreaGlacierizedShelf : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceAreaGlacierizedShelf(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_area_glacierized_floating") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("iareafl");
    }

    set_units("m^2", "m^2");
    m_variable["long_name"]     = "area of ice shelves in glacierized areas";
    m_variable["standard_name"] = "floating_ice_shelf_area";
    m_variable["valid_min"]     = { 0.0 };
  }

  double compute() {
    return ice_area_floating(model->geometry(),
                             m_config->get_number("output.ice_free_thickness_standard"));
  }
};

//! \brief Computes the total grounded ice volume.
class IceVolumeGlacierizedGrounded : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceVolumeGlacierizedGrounded(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized_grounded") {

    set_units("m^3", "m^3");
    m_variable["long_name"] = "volume of grounded ice in glacierized areas";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    const auto &cell_type = model->geometry().cell_type;

    const array::Scalar &ice_thickness = model->geometry().ice_thickness;

    const double thickness_threshold = m_config->get_number("output.ice_free_thickness_standard"),
                 cell_area           = m_grid->cell_area();

    array::AccessScope list{ &ice_thickness, &cell_type };

    double volume = 0.0;
    for (auto p : m_grid->points()) {
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
class IceVolumeGlacierizedShelf : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  IceVolumeGlacierizedShelf(IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "ice_volume_glacierized_floating") {

    set_units("m^3", "m^3");
    m_variable["long_name"] = "volume of ice shelves in glacierized areas";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    const auto &cell_type = model->geometry().cell_type;

    const array::Scalar &ice_thickness = model->geometry().ice_thickness;

    const double thickness_threshold = m_config->get_number("output.ice_free_thickness_standard"),
                 cell_area           = m_grid->cell_area();

    array::AccessScope list{ &ice_thickness, &cell_type };

    double volume = 0.0;
    for (auto p : m_grid->points()) {
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
class TimeStepLength : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  TimeStepLength(const IceModel *m) : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "dt") {

    set_units("second", "year");
    m_variable["long_name"] = "mass continuity time step";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    return model->dt();
  }
};

//! \brief Reports the mass continuity time step.
class TimeStepRatio : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  TimeStepRatio(const IceModel *m) : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "dt_ratio") {

    set_units("1", "1");
    m_variable["long_name"] = "ratio of max. allowed time steps "
        "according to CFL and SIA diffusivity criteria";
    m_variable["valid_min"] = { 0.0 };
    m_variable["_FillValue"] = { -1.0 };
  }

  double compute() {

    const auto *stress_balance = model->stress_balance();

    auto cfl_2d = model->max_timestep_cfl_2d();
    auto cfl_3d = model->max_timestep_cfl_3d();

    auto dt_diff = max_timestep_diffusivity(stress_balance->max_diffusivity(),
                                            model->grid()->dx(),
                                            model->grid()->dy(),
                                            m_config->get_number("time_stepping.adaptive_ratio"));

    auto dt_cfl = std::min(cfl_2d.dt_max, cfl_3d.dt_max);

    if (dt_cfl.finite() and dt_diff.finite() and dt_diff.value() > 0.0) {
      return dt_cfl.value() / dt_diff.value();
    }

    return -1.0;
  }
};

//! \brief Reports maximum diffusivity.
class MaxDiffusivity : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  MaxDiffusivity(const IceModel *m) : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "max_diffusivity") {

    set_units("m^2 s^-1", "m^2 s^-1");
    m_variable["long_name"] = "maximum diffusivity of the flow (usually from the SIA model)";
    m_variable["valid_min"] = { 0.0 };
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
class MaxHorizontalVelocity : public TSDiag<TSSnapshotDiagnostic, IceModel> {
public:
  MaxHorizontalVelocity(const IceModel *m)
      : TSDiag<TSSnapshotDiagnostic, IceModel>(m, "max_horizontal_vel") {

    set_units("m second^-1", "m year^-1");
    m_variable["long_name"] = "max(max(abs(u)), max(abs(v))) for the horizontal velocity of ice"
                              " over volume of the ice in last time step during time-series reporting interval";
    m_variable["valid_min"] = { 0.0 };
  }

  double compute() {
    CFLData cfl = model->max_timestep_cfl_3d();

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
  const Grid &grid     = *model->grid();
  const Config &config = *grid.ctx()->config();

  const double ice_density = config.get_number("constants.ice.density"),
               cell_area   = grid.cell_area();

  const auto &cell_type = model->geometry().cell_type;

  const array::Scalar *thickness_change = nullptr;

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

  const array::Scalar &dV_flow =
      model->geometry_evolution().area_specific_volume_change_due_to_flow();

  array::AccessScope list{ &cell_type, thickness_change };

  if (term == FLOW) {
    list.add(dV_flow);
  }

  double volume_change = 0.0;
  for (auto p : grid.points()) {
    const int i = p.i(), j = p.j();

    if ((area == BOTH) or (area == GROUNDED and cell_type.grounded(i, j)) or
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
class IceMassFluxBasal : public TSDiag<TSFluxDiagnostic, IceModel> {
public:
  IceMassFluxBasal(const IceModel *m)
      : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_basal_mass_flux") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("tendlibmassbf");
    }

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"]     = "total over ice domain of bottom surface ice mass flux";
    m_variable["standard_name"] = "tendency_of_land_ice_mass_due_to_basal_mass_balance";
    m_variable["comment"]       = "positive means ice gain";
  }

  double compute() {
    return mass_change(model, BMB, BOTH);
  }
};

//! \brief Reports the total top surface ice flux.
class IceMassFluxSurface : public TSDiag<TSFluxDiagnostic, IceModel> {
public:
  IceMassFluxSurface(const IceModel *m)
      : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_surface_mass_flux") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("tendacabf");
    }

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"]     = "total over ice domain of top surface ice mass flux";
    m_variable["standard_name"] = "tendency_of_land_ice_mass_due_to_surface_mass_balance";
    m_variable["comment"]       = "positive means ice gain";
  }

  double compute() {
    return mass_change(model, SMB, BOTH);
  }
};

//! \brief Reports the total basal ice flux over the grounded region.
class IceMassFluxBasalGrounded : public TSDiag<TSFluxDiagnostic, IceModel> {
public:
  IceMassFluxBasalGrounded(const IceModel *m)
      : TSDiag<TSFluxDiagnostic, IceModel>(m, "basal_mass_flux_grounded") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("tendlibmassbfgr");
    }

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"]     = "total over grounded ice domain of basal mass flux";
    m_variable["standard_name"] = "tendency_of_land_ice_mass_due_to_basal_mass_balance";
    m_variable["comment"]       = "positive means ice gain";
  }

  double compute() {
    return mass_change(model, BMB, GROUNDED);
  }
};

//! \brief Reports the total sub-shelf ice flux.
class IceMassFluxBasalFloating : public TSDiag<TSFluxDiagnostic, IceModel> {
public:
  IceMassFluxBasalFloating(const IceModel *m)
      : TSDiag<TSFluxDiagnostic, IceModel>(m, "basal_mass_flux_floating") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("tendlibmassbffl");
    }

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"]     = "total sub-shelf ice flux";
    m_variable["standard_name"] = "tendency_of_land_ice_mass_due_to_basal_mass_balance";
    m_variable["comment"]       = "positive means ice gain";
  }

  double compute() {
    return mass_change(model, BMB, SHELF);
  }
};

//! \brief Reports the total numerical mass flux needed to preserve
//! non-negativity of ice thickness.
class IceMassFluxConservationError : public TSDiag<TSFluxDiagnostic, IceModel> {
public:
  IceMassFluxConservationError(const IceModel *m)
      : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_conservation_error") {

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"] = "total numerical flux needed to preserve non-negativity"
                              " of ice thickness";
    m_variable["comment"]   = "positive means ice gain";
  }

  double compute() {
    return mass_change(model, ERROR, BOTH);
  }
};

//! \brief Reports the total discharge flux.
class IceMassFluxDischarge : public TSDiag<TSFluxDiagnostic, IceModel> {
public:
  IceMassFluxDischarge(const IceModel *m)
      : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_discharge") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("tendlifmassbf");
    }

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"]     = "discharge flux (frontal melt, calving, forced retreat)";
    m_variable["standard_name"] = "tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting";
    m_variable["comment"]       = "positive means ice gain";
  }

  double compute() {
    const double ice_density = m_config->get_number("constants.ice.density");

    const array::Scalar &calving        = model->calving();
    const array::Scalar &frontal_melt   = model->frontal_melt();
    const array::Scalar &forced_retreat = model->forced_retreat();

    auto cell_area = m_grid->cell_area();

    double volume_change = 0.0;

    array::AccessScope list{ &calving, &frontal_melt, &forced_retreat };

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      // m^2 * m = m^3
      volume_change += cell_area * (calving(i, j) + frontal_melt(i, j) + forced_retreat(i, j));
    }

    // (kg/m^3) * m^3 = kg
    return ice_density * GlobalSum(m_grid->com, volume_change);
  }
};

//! \brief Reports the total calving flux.
class IceMassFluxCalving : public TSDiag<TSFluxDiagnostic, IceModel> {
public:
  IceMassFluxCalving(const IceModel *m)
      : TSDiag<TSFluxDiagnostic, IceModel>(m, "tendency_of_ice_mass_due_to_calving") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("tendlicalvf");
    }

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"]     = "calving flux";
    m_variable["standard_name"] = "tendency_of_land_ice_mass_due_to_calving";
    m_variable["comment"]       = "positive means ice gain";
  }

  double compute() {
    const double ice_density = m_config->get_number("constants.ice.density");

    const array::Scalar &calving = model->calving();

    auto cell_area = m_grid->cell_area();

    double volume_change = 0.0;

    array::AccessScope list{ &calving };

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      // m^2 * m = m^3
      volume_change += cell_area * calving(i, j);
    }

    // (kg/m^3) * m^3 = kg
    return ice_density * GlobalSum(m_grid->com, volume_change);
  }
};

//! @brief Reports the total flux across the grounding line.
class IceMassFluxAtGroundingLine : public TSDiag<TSFluxDiagnostic, IceModel> {
public:
  IceMassFluxAtGroundingLine(const IceModel *m)
      : TSDiag<TSFluxDiagnostic, IceModel>(m, "grounding_line_flux") {

    if (m_config->get_flag("output.ISMIP")) {
      m_variable.set_name("tendligroundf");
      m_variable["standard_name"] = "tendency_of_grounded_ice_mass";
    }

    set_units("kg s^-1", "Gt year^-1");
    m_variable["long_name"] = "total ice flux across the grounding line";
    m_variable["comment"]   = "negative flux corresponds to ice loss into the ocean";
  }

  double compute() {
    return total_grounding_line_flux(model->geometry().cell_type,
                                     model->geometry_evolution().flux_staggered(), model->dt());
  }
};

} // end of namespace scalar


//! \brief Computes dHdt, the ice thickness rate of change.
class ThicknessRateOfChange : public Diag<IceModel> {
public:
  ThicknessRateOfChange(const IceModel *m)
      : Diag<IceModel>(m), m_last_thickness(m_grid, "last_ice_thickness"), m_interval_length(0.0) {

    auto ismip = m_config->get_flag("output.ISMIP");

    // set metadata:
    m_vars = { { m_sys, ismip ? "dlithkdt" : "dHdt", *m_grid } };
    m_vars[0]
        .long_name("ice thickness rate of change")
        .standard_name("tendency_of_land_ice_thickness")
        .units("m s^-1")
        .output_units("m year^-1");

    m_vars[0]["_FillValue"]   = { fill_value() };
    m_vars[0]["cell_methods"] = "time: mean";

    m_last_thickness.metadata(0)
        .long_name(
            "ice thickness at the time of the last report of the rate of change of ice thickness")
        .units("m")
        .standard_name("land_ice_thickness");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const {

    auto result        = std::make_shared<array::Scalar>(m_grid, "dHdt");
    result->metadata() = m_vars[0];

    if (m_interval_length > 0.0) {
      geometry.ice_thickness.add(-1.0, m_last_thickness, *result);
      result->scale(1.0 / m_interval_length);
    } else {
      result->set(fill_value());
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

  array::Scalar m_last_thickness;
  double m_interval_length;
};

//! \brief Computes latitude and longitude bounds.
class LatLonBounds : public Diag<IceModel> {
public:
  LatLonBounds(const IceModel *m, const std::string &var_name, const std::string &proj_string);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;

  std::string m_var_name, m_proj_string;
};

LatLonBounds::LatLonBounds(const IceModel *m, const std::string &var_name,
                           const std::string &proj_string)
    : Diag<IceModel>(m) {
  assert(var_name == "lat" || var_name == "lon");
  m_var_name = var_name;

  // set metadata:
  m_vars = { { m_sys, m_var_name + "_bnds", *m_grid, { 0.0, 1.0, 2.0, 3.0 } } };
  m_vars[0].dimension("z").clear().set_name("nv4");

  m_vars[0].set_time_dependent(false);
  if (m_var_name == "lon") {
    m_vars[0].long_name("longitude bounds").units("degree_east");
    m_vars[0]["valid_range"] = { -180, 180 };
  } else {
    m_vars[0].long_name("latitude bounds").units("degree_north");
    m_vars[0]["valid_range"] = { -90, 90 };
  }

  m_proj_string = proj_string;

#if (Pism_USE_PROJ == 1)
  // create the transformation from the provided projection to lat,lon to check if
  // proj_string is valid.
  Proj crs(m_proj_string, "EPSG:4326");
#endif
  // If PISM_USE_PROJ is not 1 we don't need to check validity of m_proj_string: this diagnostic
  // will not be available and so this code will not run.
}

std::shared_ptr<array::Array> LatLonBounds::compute_impl(const Geometry &/*geometry*/) const {
  std::shared_ptr<array::Array3D> result(new array::Array3D(
      m_grid, m_var_name + "_bnds", array::WITHOUT_GHOSTS, { 0.0, 1.0, 2.0, 3.0 }));
  result->metadata(0) = m_vars[0];

  if (m_var_name == "lat") {
    compute_lat_bounds(m_proj_string, *result);
  } else {
    compute_lon_bounds(m_proj_string, *result);
  }

  return result;
}

class IceAreaFraction : public Diag<IceModel> {
public:
  IceAreaFraction(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

IceAreaFraction::IceAreaFraction(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, land_ice_area_fraction_name, *m_grid } };
  m_vars[0]
      .long_name("fraction of a grid cell covered by ice (grounded or floating)")
      .standard_name("land_ice_area_fraction") // InitMIP "standard" name
      .units("1");
}

std::shared_ptr<array::Array> IceAreaFraction::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Scalar>(land_ice_area_fraction_name);

  const array::Scalar1 &thickness         = geometry.ice_thickness,
                       &surface_elevation = geometry.ice_surface_elevation,
                       &bed_topography    = geometry.bed_elevation;

  const array::CellType1 &cell_type = geometry.cell_type;

  array::AccessScope list{ &thickness, &surface_elevation, &bed_topography, &cell_type,
                           result.get() };

  const bool do_part_grid   = m_config->get_flag("geometry.part_grid.enabled");
  const array::Scalar &Href = geometry.ice_area_specific_volume;
  ;
  if (do_part_grid) {
    list.add(Href);
  }

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      if (cell_type.icy(i, j)) {
        // an "icy" cell: the area fraction is one
        (*result)(i, j) = 1.0;
      } else if (cell_type.ice_free_ocean(i, j)) {
        // an ice-free ocean cell may be "partially-filled", in which case we need to compute its
        // ice area fraction by dividing Href by the threshold thickness.

        double H_reference = do_part_grid ? Href(i, j) : 0.0;

        if (H_reference > 0.0) {
          const double H_threshold =
              part_grid_threshold_thickness(cell_type.star_int(i, j), thickness.star(i, j),
                                            surface_elevation.star(i, j), bed_topography(i, j));
          // protect from a division by zero
          if (H_threshold > 0.0) {
            (*result)(i, j) = std::min(H_reference / H_threshold, 1.0);
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

class IceAreaFractionGrounded : public Diag<IceModel> {
public:
  IceAreaFractionGrounded(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

IceAreaFractionGrounded::IceAreaFractionGrounded(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, grounded_ice_sheet_area_fraction_name, *m_grid } };
  m_vars[0]
      .long_name("fraction of a grid cell covered by grounded ice")
      .standard_name("grounded_ice_sheet_area_fraction") // InitMIP "standard" name
      .units("1");
}

std::shared_ptr<array::Array> IceAreaFractionGrounded::compute_impl(const Geometry &geometry) const {
  auto result = std::make_shared<array::Scalar>(m_grid, grounded_ice_sheet_area_fraction_name);
  result->metadata() = m_vars[0];

  const double ice_density   = m_config->get_number("constants.ice.density"),
               ocean_density = m_config->get_number("constants.sea_water.density");

  const auto &ice_thickness  = geometry.ice_thickness;
  const auto &sea_level      = geometry.sea_level_elevation;
  const auto &bed_topography = geometry.bed_elevation;

  const auto &cell_type = geometry.cell_type;

  compute_grounded_cell_fraction(ice_density, ocean_density, sea_level, ice_thickness,
                                 bed_topography, *result);

  // All grounded areas have the grounded cell fraction of one, so now we make sure that ice-free
  // areas get the value of 0 (they are grounded but not covered by a grounded ice sheet).

  array::AccessScope list{ &cell_type, result.get() };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
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

class IceAreaFractionFloating : public Diag<IceModel> {
public:
  IceAreaFractionFloating(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

IceAreaFractionFloating::IceAreaFractionFloating(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, floating_ice_sheet_area_fraction_name, *m_grid } };
  m_vars[0]
      .long_name("fraction of a grid cell covered by floating ice")
      .standard_name("floating_ice_shelf_area_fraction")
      .units("1");
}

std::shared_ptr<array::Array> IceAreaFractionFloating::compute_impl(const Geometry &geometry) const {

  auto ice_area_fraction      = IceAreaFraction(model).compute(geometry);
  auto grounded_area_fraction = IceAreaFractionGrounded(model).compute(geometry);

  auto result        = ice_area_fraction;
  result->metadata() = m_vars[0];

  // Floating area fraction is total area fraction minus grounded area fraction.
  result->add(-1.0, *grounded_area_fraction);

  return result;
}

//! \brief Computes the 2D height above flotation.
class HeightAboveFloatation : public Diag<IceModel> {
public:
  HeightAboveFloatation(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

HeightAboveFloatation::HeightAboveFloatation(const IceModel *m) : Diag<IceModel>(m) {

  // set metadata:
  m_vars = { { m_sys, "height_above_flotation", *m_grid } };
  m_vars[0].long_name("ice thickness in excess of the maximum floating ice thickness").units("m");
  m_vars[0]["_FillValue"] = { fill_value() };
  m_vars[0]["comment"]    = "shows how close to floatation the ice is at a given location";
}

std::shared_ptr<array::Array> HeightAboveFloatation::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Scalar>("height_above_flotation");

  const auto &cell_type = geometry.cell_type;

  const double ice_density   = m_config->get_number("constants.ice.density"),
               ocean_density = m_config->get_number("constants.sea_water.density");
  double fill = fill_value();

  const auto &sea_level      = geometry.sea_level_elevation;
  const auto &ice_thickness  = geometry.ice_thickness;
  const auto &bed_topography = geometry.bed_elevation;

  array::AccessScope list{ &cell_type, result.get(), &ice_thickness, &bed_topography, &sea_level };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      const double thickness = ice_thickness(i, j), bed = bed_topography(i, j),
                   ocean_depth = sea_level(i, j) - bed;

      if (cell_type.icy(i, j) and ocean_depth > 0.0) {
        const double max_floating_thickness = ocean_depth * (ocean_density / ice_density);
        (*result)(i, j)                     = thickness - max_floating_thickness;
      } else {
        (*result)(i, j) = fill;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

//! \brief Computes the mass per cell.
class IceMass : public Diag<IceModel> {
public:
  IceMass(const IceModel *m);

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

IceMass::IceMass(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "ice_mass", *m_grid } };
  m_vars[0].long_name("ice mass per cell").units("kg");
  m_vars[0]["_FillValue"] = { fill_value() };
}

std::shared_ptr<array::Array> IceMass::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Scalar>("ice_mass");

  const auto &cell_type = geometry.cell_type;

  const double ice_density = m_config->get_number("constants.ice.density");

  const array::Scalar &ice_thickness = geometry.ice_thickness;

  auto cell_area = m_grid->cell_area();

  array::AccessScope list{ &cell_type, result.get(), &ice_thickness };

  double fill = fill_value();

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      // count all ice, including cells which have so little they
      // are considered "ice-free"
      if (ice_thickness(i, j) > 0.0) {
        (*result)(i, j) = ice_density * ice_thickness(i, j) * cell_area;
      } else {
        (*result)(i, j) = fill;
      }
    } // end of loop over grid points

  } catch (...) {
    loop.failed();
  }
  loop.check();

  // Add the mass of ice in Href:
  if (m_config->get_flag("geometry.part_grid.enabled")) {
    const array::Scalar &Href = geometry.ice_area_specific_volume;
    list.add(Href);
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      if (ice_thickness(i, j) <= 0.0 and Href(i, j) > 0.0) {
        (*result)(i, j) = ice_density * Href(i, j) * cell_area;
      }
    }
  }

  return result;
}

/*! @brief Sea-level adjusted bed topography (zero at sea level). */
class BedTopographySeaLevelAdjusted : public Diag<IceModel> {
public:
  BedTopographySeaLevelAdjusted(const IceModel *m);

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

BedTopographySeaLevelAdjusted::BedTopographySeaLevelAdjusted(const IceModel *m)
    : Diag<IceModel>(m) {
  m_vars = { { m_sys, "topg_sl_adjusted", *m_grid } };
  m_vars[0].long_name("sea-level adjusted bed topography (zero is at sea level)").units("meters");
}

std::shared_ptr<array::Array> BedTopographySeaLevelAdjusted::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Scalar>("topg_sl_adjusted");

  const auto &bed       = geometry.bed_elevation;
  const auto &sea_level = geometry.sea_level_elevation;

  array::AccessScope list{ &bed, &sea_level, result.get() };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    (*result)(i, j) = bed(i, j) - sea_level(i, j);
  }

  return result;
}

/*! @brief Ice hardness computed using the SIA flow law. */
class IceHardness : public Diag<IceModel> {
public:
  IceHardness(const IceModel *m);

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

IceHardness::IceHardness(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "hardness", *m_grid, m_grid->z() } };
  m_vars[0]
      .long_name("ice hardness computed using the SIA flow law")
      .set_units_without_validation(
          "Pa s^(1/n)"); // n is the Glen exponent used by the SIA (modifier) flow law
  m_vars[0]["comment"] = "units depend on the Glen exponent used by the flow law";
}

std::shared_ptr<array::Array> IceHardness::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "hardness", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata(0) = m_vars[0];

  auto EC = m_grid->ctx()->enthalpy_converter();

  const array::Array3D &ice_enthalpy = model->energy_balance_model()->enthalpy();
  const array::Scalar &ice_thickness = geometry.ice_thickness;

  const rheology::FlowLaw &flow_law = *model->stress_balance()->modifier()->flow_law();

  array::AccessScope list{ &ice_enthalpy, &ice_thickness, result.get() };

  const unsigned int Mz = m_grid->Mz();

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      const double *E = ice_enthalpy.get_column(i, j);
      const double H  = ice_thickness(i, j);

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

/*! @brief Effective viscosity of ice (3D). */
class IceViscosity : public Diag<IceModel> {
public:
  IceViscosity(IceModel *m);

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

IceViscosity::IceViscosity(IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "effective_viscosity", *m_grid, m_grid->z() } };
  m_vars[0]
      .long_name("effective viscosity of ice")
      .units("Pascal second")
      .output_units("kPascal second");
  m_vars[0]["valid_min"]  = { 0.0 };
  m_vars[0]["_FillValue"] = { fill_value() };
}

static inline double square(double x) {
  return x * x;
}

std::shared_ptr<array::Array> IceViscosity::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "effective_viscosity", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata(0) = m_vars[0];

  array::Array3D W(m_grid, "wvel", array::WITH_GHOSTS, m_grid->z());

  using mask::ice_free;

  auto EC = m_grid->ctx()->enthalpy_converter();

  const rheology::FlowLaw &flow_law = *model->stress_balance()->modifier()->flow_law();

  const array::Scalar &ice_thickness = geometry.ice_thickness;

  const array::Array3D &ice_enthalpy     = model->energy_balance_model()->enthalpy(),
                       &U                = model->stress_balance()->velocity_u(),
                       &V                = model->stress_balance()->velocity_v(),
                       &W_without_ghosts = model->vertical_velocity();

  W.copy_from(W_without_ghosts);

  const unsigned int Mz = m_grid->Mz();
  const double dx = m_grid->dx(), dy = m_grid->dy();
  const std::vector<double> &z = m_grid->z();

  double fill = fill_value();

  const array::CellType1 &mask = geometry.cell_type;

  array::AccessScope list{ &U, &V, &W, &ice_enthalpy, &ice_thickness, &mask, result.get() };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      const double *E = ice_enthalpy.get_column(i, j);
      const double H  = ice_thickness(i, j);

      const double *u = U.get_column(i, j), *u_n = U.get_column(i, j + 1),
                   *u_e = U.get_column(i + 1, j), *u_s = U.get_column(i, j - 1),
                   *u_w = U.get_column(i - 1, j);

      const double *v = V.get_column(i, j), *v_n = V.get_column(i, j + 1),
                   *v_e = V.get_column(i + 1, j), *v_s = V.get_column(i, j - 1),
                   *v_w = V.get_column(i - 1, j);

      const double *w = W.get_column(i, j), *w_n = W.get_column(i, j + 1),
                   *w_e = W.get_column(i + 1, j), *w_s = W.get_column(i, j - 1),
                   *w_w = W.get_column(i - 1, j);

      auto m                  = mask.star_int(i, j);
      const unsigned int east = ice_free(m.e) ? 0 : 1, west = ice_free(m.w) ? 0 : 1,
                         south = ice_free(m.s) ? 0 : 1, north = ice_free(m.n) ? 0 : 1;

      double *viscosity = result->get_column(i, j);

      if (ice_free(m.c)) {
        result->set_column(i, j, fill);
        continue;
      }

      for (unsigned int k = 0; k < Mz; ++k) {
        const double depth = H - z[k];

        if (depth < 0.0) {
          viscosity[k] = fill;
          continue;
        }

        // EC->pressure() handles negative depths correctly
        const double pressure = EC->pressure(depth);

        const double hardness = flow_law.hardness(E[k], pressure);

        double u_x = 0.0, v_x = 0.0, w_x = 0.0;
        if (west + east > 0) {
          const double D = 1.0 / (dx * (west + east));
          u_x            = D * (west * (u[k] - u_w[k]) + east * (u_e[k] - u[k]));
          v_x            = D * (west * (v[k] - v_w[k]) + east * (v_e[k] - v[k]));
          w_x            = D * (west * (w[k] - w_w[k]) + east * (w_e[k] - w[k]));
        }

        double u_y = 0.0, v_y = 0.0, w_y = 0.0;
        if (south + north > 0) {
          const double D = 1.0 / (dy * (south + north));
          u_y            = D * (south * (u[k] - u_s[k]) + north * (u_n[k] - u[k]));
          v_y            = D * (south * (v[k] - v_s[k]) + north * (v_n[k] - v[k]));
          w_y            = D * (south * (w[k] - w_s[k]) + north * (w_n[k] - w[k]));
        }

        double u_z = 0.0, v_z = 0.0, w_z = 0.0;

        if (k == 0) {
          const double dz = z[1] - z[0];
          u_z             = (u[1] - u[0]) / dz;
          v_z             = (v[1] - v[0]) / dz;
          w_z             = (w[1] - w[0]) / dz;
        } else if (k == Mz - 1) {
          const double dz = z[Mz - 1] - z[Mz - 2];
          u_z             = (u[Mz - 1] - u[Mz - 2]) / dz;
          v_z             = (v[Mz - 1] - v[Mz - 2]) / dz;
          w_z             = (w[Mz - 1] - w[Mz - 2]) / dz;
        } else {
          const double dz_p = z[k + 1] - z[k], dz_m = z[k] - z[k - 1];
          u_z = 0.5 * ((u[k + 1] - u[k]) / dz_p + (u[k] - u[k - 1]) / dz_m);
          v_z = 0.5 * ((v[k + 1] - v[k]) / dz_p + (v[k] - v[k - 1]) / dz_m);
          w_z = 0.5 * ((w[k + 1] - w[k]) / dz_p + (w[k] - w[k - 1]) / dz_m);
        }

        // These should be "epsilon dot", but that's just too long.
        const double eps_xx = u_x, eps_yy = v_y, eps_zz = w_z, eps_xy = 0.5 * (u_y + v_x),
                     eps_xz = 0.5 * (u_z + w_x), eps_yz = 0.5 * (v_z + w_y);

        // The second invariant of the 3D strain rate tensor; see equation 4.8 in [@ref
        // GreveBlatter2009]. Unlike secondInvariant_2D(), this code does not make assumptions about
        // the input velocity field: we do not ignore w_x and w_y and do not assume that u_z and v_z
        // are zero.
        const double gamma = (square(eps_xx) + square(eps_yy) + square(eps_zz) +
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

/*! @brief Report ice thickness */
class IceThickness : public Diag<IceModel> {
public:
  IceThickness(const IceModel *m) : Diag<IceModel>(m) {

    auto ismip = m_config->get_flag("output.ISMIP");

    m_vars = { { m_sys, ismip ? "lithk" : "thk", *m_grid } };

    m_vars[0].long_name("land ice thickness").standard_name("land_ice_thickness").units("m");
    m_vars[0]["valid_min"] = { 0.0 };
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const {

    auto result = allocate<array::Scalar>("thk");

    result->copy_from(geometry.ice_thickness);

    return result;
  }
};

/*! @brief Report ice top surface elevation */
class IceBottomSurfaceElevation : public Diag<IceModel> {
public:
  IceBottomSurfaceElevation(const IceModel *m) : Diag<IceModel>(m) {

    auto ismip = m_config->get_flag("output.ISMIP");

    m_vars = { { m_sys, ismip ? "base" : "ice_base_elevation", *m_grid } };
    m_vars[0].long_name("ice bottom surface elevation").units("m");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const {

    auto result = allocate<array::Scalar>("ice_base_elevation");

    ice_bottom_surface(geometry, *result);

    return result;
  }
};

/*! @brief Report ice top surface elevation */
class IceSurfaceElevation : public Diag<IceModel> {
public:
  IceSurfaceElevation(const IceModel *m) : Diag<IceModel>(m) {

    auto ismip = m_config->get_flag("output.ISMIP");

    m_vars = { { m_sys, ismip ? "orog" : "usurf", *m_grid } };
    m_vars[0].long_name("ice top surface elevation").standard_name("surface_altitude").units("m");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const {
    auto result = allocate<array::Scalar>("usurf");

    result->copy_from(geometry.ice_surface_elevation);

    return result;
  }
};

/*! @brief Report grounding line flux. */
class GroundingLineFlux : public DiagAverageRate<IceModel> {
public:
  GroundingLineFlux(const IceModel *m) : DiagAverageRate<IceModel>(m, "grounding_line_flux", RATE) {

    m_accumulator.metadata()["units"] = "kg m^-2";

    auto ismip = m_config->get_flag("output.ISMIP");

    m_vars = { { m_sys, ismip ? "ligroundf" : "grounding_line_flux", *m_grid } };

    m_vars[0]
        .long_name("grounding line flux")
        .units("kg m^-2 second^-1")
        .output_units("kg m^-2 year^-1");
    m_vars[0]["cell_methods"] = "time: mean";

    m_vars[0]["_FillValue"] = { fill_value() };
    m_vars[0]["comment"] =
      "Positive flux corresponds to mass moving from the ocean to"
      " an icy grounded area. This convention makes it easier to compare"
      " grounding line flux to the total discharge into the ocean";
  }

protected:
  void update_impl(double dt) {
    auto grid      = m_accumulator.grid();
    auto cell_area = grid->cell_area(); // units: m^2
    auto ice_density =
        grid->ctx()->config()->get_number("constants.ice.density"); // units: kg / m^3

    // factor used to convert from m^3/s to kg/m^2
    double unit_conversion_factor = dt * (ice_density / cell_area); // units: kg * s / m^5

    ice_flow_rate_across_grounding_line(model->geometry().cell_type,
                                        model->geometry_evolution().flux_staggered(),
                                        unit_conversion_factor, m_accumulator);

    m_interval_length += dt;
  }
};

class MassTransportAcrossGroundingLine : public DiagAverageRate<IceModel> {
public:
  MassTransportAcrossGroundingLine(const IceModel *m)
      : DiagAverageRate<IceModel>(m, "ice_mass_transport_across_grounding_line", RATE) {

    m_accumulator.metadata()["units"] = "kg";

    m_vars = { { m_sys, "ice_mass_transport_across_grounding_line", *m_grid } };

    m_vars[0]
        .long_name("ice mass flow rate across the grounding line")
        .units("kg s^-1")
        .output_units("Gt year^-1");
    m_vars[0]["cell_methods"] = "time: mean";

    m_vars[0]["_FillValue"] = { fill_value() };
    m_vars[0]["comment"] =
        "Negative values correspond to mass moving from an icy grounded area into a lake or ocean."
        " This convention makes it easier to compare to calving, frontal melt, and discharge fluxes.";
  }

protected:
  void update_impl(double dt) {
    auto grid = model->grid();
    double ice_density =
        grid->ctx()->config()->get_number("constants.ice.density"); // units: kg / m^3

    // factor used to convert from m^3/s to kg
    double unit_conversion_factor = dt * ice_density; // units: kg * s / m^3

    ice_flow_rate_across_grounding_line(model->geometry().cell_type,
                                        model->geometry_evolution().flux_staggered(),
                                        unit_conversion_factor, m_accumulator);

    m_interval_length += dt;
  }
};

//! \brief Reports the pressure within the ice (3D).
class PressureInIce : public Diag<IceModel>
{
public:
  PressureInIce(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};


PressureInIce::PressureInIce(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "pressure", *m_grid, m_grid->z() } };
  m_vars[0].long_name("pressure in ice (hydrostatic)").units("Pa");
}

std::shared_ptr<array::Array> PressureInIce::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "pressure", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata(0) = m_vars[0];

  const array::Scalar &thickness = geometry.ice_thickness;

  array::AccessScope list{ &thickness, result.get() };

  const double rg = m_config->get_number("constants.ice.density") *
                    m_config->get_number("constants.standard_gravity");

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      unsigned int ks  = m_grid->kBelowHeight(thickness(i, j));
      double *P_out_ij = result->get_column(i, j);
      const double H   = thickness(i, j);
      // within the ice:
      for (unsigned int k = 0; k <= ks; ++k) {
        P_out_ij[k] = rg * (H - m_grid->z(k)); // FIXME: add atmospheric pressure?
      }
      // above the ice:
      for (unsigned int k = ks + 1; k < m_grid->Mz(); ++k) {
        P_out_ij[k] = 0.0; // FIXME: use atmospheric pressure?
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

//! \brief Computes the gravitational driving stress (diagnostically).
class DrivingShearStress : public Diag<IceModel>
{
public:
  DrivingShearStress(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

DrivingShearStress::DrivingShearStress(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = { { m_sys, "taud_x", *m_grid }, { m_sys, "taud_y", *m_grid } };
  m_vars[0].long_name("X-component of the driving shear stress at the base of ice");
  m_vars[1].long_name("Y-component of the driving shear stress at the base of ice");

  for (auto &v : m_vars) {
    v.units("Pa");
    v["comment"] = "this field is purely diagnostic (not used by the model)";
  }
}

/*!
 * The driving stress computed here is not used by the model, so this
 * implementation intentionally does not use the eta-transformation or special
 * cases at ice margins.
 */
std::shared_ptr<array::Array> DrivingShearStress::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Vector>("taud");

  const array::Scalar &thickness = geometry.ice_thickness;
  const array::Scalar &surface   = geometry.ice_surface_elevation;

  double standard_gravity = m_config->get_number("constants.standard_gravity"),
         ice_density      = m_config->get_number("constants.ice.density");

  array::AccessScope list{ &surface, &thickness, result.get() };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    double pressure = ice_density * standard_gravity * thickness(i, j);
    if (pressure <= 0.0) {
      (*result)(i, j).u = 0.0;
      (*result)(i, j).v = 0.0;
    } else {
      (*result)(i, j).u = -pressure * diff_x_p(surface, i, j);
      (*result)(i, j).v = -pressure * diff_y_p(surface, i, j);
    }
  }

  return result;
}

//! @brief Computes the basal shear stress @f$ \tau_b @f$.
class BasalShearStress : public Diag<IceModel>
{
public:
  BasalShearStress(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};


BasalShearStress::BasalShearStress(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "taub_x", *m_grid }, { m_sys, "taub_y", *m_grid } };

  m_vars[0].long_name("X-component of the shear stress at the base of ice");
  m_vars[1].long_name("Y-component of the shear stress at the base of ice");

  for (auto &v : m_vars) {
    v.units("Pa");
    v["comment"] = "this field is purely diagnostic (not used by the model)";
  }
}


std::shared_ptr<array::Array> BasalShearStress::compute_impl(const Geometry &geometry) const {

  auto result = allocate<array::Vector>("taub");

  const auto *yield_stress_model = model->basal_yield_stress_model();

  if (yield_stress_model == nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot compute 'taub': no yield stress model available");
  }

  const auto &velocity = model->stress_balance()->shallow()->velocity();
  const auto &tauc     = yield_stress_model->basal_material_yield_stress();
  const auto &mask     = geometry.cell_type;

  const auto *basal_sliding_law = model->stress_balance()->shallow()->sliding_law();

  array::AccessScope list{ &tauc, &velocity, &mask, result.get() };
  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (mask.grounded_ice(i, j)) {
      double beta     = basal_sliding_law->drag(tauc(i, j), velocity(i, j).u, velocity(i, j).v);
      (*result)(i, j) = -beta * velocity(i, j);
    } else {
      (*result)(i, j) = 0.0;
    }
  }

  return result;
}

//! \brief Computes the magnitude of the gravitational driving stress
//! (diagnostically).
class DrivingShearStressmagnitude : public Diag<IceModel>
{
public:
  DrivingShearStressmagnitude(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

DrivingShearStressmagnitude::DrivingShearStressmagnitude(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "taud_mag", *m_grid } };
  m_vars[0]
      .long_name("magnitude of the gravitational driving stress at the base of ice")
      .units("Pa");
  m_vars[0]["comment"] = "this field is purely diagnostic (not used by the model)";
}

std::shared_ptr<array::Array> DrivingShearStressmagnitude::compute_impl(const Geometry &geometry) const {
  auto result = allocate<array::Scalar>("taud_mag");
  auto taud = array::cast<array::Vector>(DrivingShearStress(model).compute(geometry));

  compute_magnitude(*taud, *result);

  return result;
}

//! \brief Computes the magnitude of the basal shear stress
//! (diagnostically).
class BasalShearStressMagnitude : public Diag<IceModel>
{
public:
  BasalShearStressMagnitude(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

BasalShearStressMagnitude::BasalShearStressMagnitude(const IceModel *m) : Diag<IceModel>(m) {

  auto ismip = m_config->get_flag("output.ISMIP");

  m_vars = { { m_sys, ismip ? "strbasemag" : "taub_mag", *m_grid } };
  m_vars[0]
      .long_name("magnitude of the basal shear stress at the base of ice")
      .standard_name("land_ice_basal_drag") // ISMIP "standard" name
      .units("Pa");
  m_vars[0]["comment"] = "this field is purely diagnostic (not used by the model)";
}

std::shared_ptr<array::Array> BasalShearStressMagnitude::compute_impl(const Geometry &geometry) const {
  auto result = allocate<array::Scalar>("taub_mag");

  auto taub = array::cast<array::Vector>(BasalShearStress(model).compute(geometry));

  compute_magnitude(*taub, *result);

  return result;
}

//! \brief Computes uflux and vflux, components of vertically-integrated horizontal
//! flux of ice.
class StressBalanceFlux : public Diag<IceModel>
{
public:
  StressBalanceFlux(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes flux_mag, the magnitude of vertically-integrated horizontal
//! flux of ice.
class StressBalanceFluxMag : public Diag<IceModel>
{
public:
  StressBalanceFluxMag(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes the vertically-averaged ice velocity.
class StressBalanceVelbar : public Diag<IceModel>
{
public:
  StressBalanceVelbar(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes velbar_mag, the magnitude of vertically-integrated horizontal
//! velocity of ice and masks out ice-free areas.
class StressBalanceVelbarMag : public Diag<IceModel>
{
public:
  StressBalanceVelbarMag(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Reports the vertically-integrated (2D) principal strain rates.
class StressBalanceStrainRates : public Diag<IceModel>
{
public:
  StressBalanceStrainRates(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Reports the vertically-integrated (2D) deviatoric stresses.
class StressBalanceDeviatoricStresses : public Diag<IceModel>
{
public:
  StressBalanceDeviatoricStresses(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

/*! @brief tensile von Mises stress */
class StressBalanceVonmisesStress : public Diag<IceModel>
{
public:
  StressBalanceVonmisesStress(const IceModel *m);
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes the x-component of the horizontal ice velocity.
class StressBalanceUvel : public Diag<IceModel>
{
public:
  StressBalanceUvel(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes the y-component of the horizontal ice velocity.
class StressBalanceVvel : public Diag<IceModel>
{
public:
  StressBalanceVvel(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes vertical velocity of ice, relative to the bed directly
//! below.
class StressBalanceWvelRel : public Diag<IceModel>
{
public:
  StressBalanceWvelRel(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Reports the xz component of the shear stress within the ice (3D), according to the SIA formula.
class StressBalanceTauxz : public Diag<IceModel>
{
public:
  StressBalanceTauxz(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Reports the yz component of the shear stress within the ice (3D), according to the SIA formula.
class StressBalanceTauyz : public Diag<IceModel>
{
public:
  StressBalanceTauyz(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes velbase_mag, the magnitude of horizontal velocity of ice at base
//! of ice and masks out ice-free areas.
class StressBalanceVelbaseMag : public Diag<IceModel>
{
public:
  StressBalanceVelbaseMag(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes velsurf_mag, the magnitude of horizontal ice velocity at the
//! surface.
class StressBalanceVelsurfMag : public Diag<IceModel>
{
public:
  StressBalanceVelsurfMag(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes velsurf, the horizontal velocity of ice at ice surface.
class StressBalanceVelsurf : public Diag<IceModel>
{
public:
  StressBalanceVelsurf(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! Computes vertical ice velocity (relative to the geoid).
/*!
  \f[
  w(s) = \tilde w(s) + \frac{\partial b}{\partial t} + U(s) \cdot \nabla b
  \f]
  in grounded areas. In floating shelves
  \f[
  w(s) = \tilde w(s) - \tilde  w(z_{\text{sea level}}).
  \f]

  This ensures that \f$\tilde w(z_{\text{sea level}}) = 0\f$.
*/
class StressBalanceWvel : public Diag<IceModel>
{
public:
  StressBalanceWvel(const IceModel *m);
  virtual std::shared_ptr<array::Array> compute(const Geometry &geometry, bool zero_above_ice) const;
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! Computes wvelsurf, the vertical velocity of ice at ice surface.
class StressBalanceWvelsurf : public Diag<IceModel>
{
public:
  StressBalanceWvelsurf(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! Computes wvelbase, the vertical velocity of ice at the base of ice.
class StressBalanceWvelbase : public Diag<IceModel>
{
public:
  StressBalanceWvelbase(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

//! \brief Computes horizontal ice velocity at the base of ice.
class StressBalanceVelbase : public Diag<IceModel>
{
public:
  StressBalanceVelbase(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

class ShallowStressBalanceBeta : public Diag<IceModel>
{
public:
  ShallowStressBalanceBeta(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

StressBalanceUvel::StressBalanceUvel(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "uvel", *m_grid } };
  m_vars[0]
      .long_name("horizontal velocity of ice in the X direction")
      .standard_name("land_ice_x_velocity")
      .units("m s^-1")
      .output_units("m year^-1");
}

/*!
 * Copy F to result and set it to zero above the surface of the ice.
 */
static void zero_above_ice(const array::Array3D &F, const array::Scalar &H,
                           array::Array3D &result) {

  array::AccessScope list{ &F, &H, &result };

  auto grid = result.grid();

  auto Mz = grid->Mz();

  ParallelSection loop(grid->com);
  try {
    for (auto p : grid->points()) {
      const int i = p.i(), j = p.j();

      int ks = grid->kBelowHeight(H(i, j));

      const double *F_ij = F.get_column(i, j);
      double *F_out_ij   = result.get_column(i, j);

      // in the ice:
      for (int k = 0; k <= ks; k++) {
        F_out_ij[k] = F_ij[k];
      }
      // above the ice:
      for (unsigned int k = ks + 1; k < Mz; k++) {
        F_out_ij[k] = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

std::shared_ptr<array::Array> StressBalanceUvel::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "uvel", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata() = m_vars[0];

  zero_above_ice(model->stress_balance()->velocity_u(), geometry.ice_thickness, *result);

  return result;
}

StressBalanceVvel::StressBalanceVvel(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "vvel", *m_grid } };
  m_vars[0]
      .long_name("horizontal velocity of ice in the Y direction")
      .standard_name("land_ice_y_velocity")
      .units("m s^-1")
      .output_units("m year^-1");
}

std::shared_ptr<array::Array> StressBalanceVvel::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "vvel", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata() = m_vars[0];

  zero_above_ice(model->stress_balance()->velocity_v(), geometry.ice_thickness, *result);

  return result;
}

StressBalanceWvelRel::StressBalanceWvelRel(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "wvel_rel", *m_grid } };
  m_vars[0]
      .long_name("vertical velocity of ice, relative to base of ice directly below")
      .units("m s^-1")
      .output_units("m year^-1");
}

std::shared_ptr<array::Array> StressBalanceWvelRel::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "wvel_rel", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata() = m_vars[0];

  zero_above_ice(model->vertical_velocity(), geometry.ice_thickness, *result);

  return result;
}

StressBalanceStrainRates::StressBalanceStrainRates(const IceModel *m) : Diag<IceModel>(m) {

  m_vars = { { m_sys, "eigen1", *m_grid }, { m_sys, "eigen2", *m_grid } };
  m_vars[0]
      .long_name("first eigenvalue of the horizontal, vertically-integrated strain rate tensor")
      .units("s^-1");
  m_vars[1]
      .long_name("second eigenvalue of the horizontal, vertically-integrated strain rate tensor")
      .units("s^-1");
}

std::shared_ptr<array::Array> StressBalanceStrainRates::compute_impl(const Geometry &geometry) const {

  auto result = std::make_shared<array::Array2D<stressbalance::PrincipalStrainRates> >(
      m_grid, "strain_rates", array::WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  array::Vector1 velbar_with_ghosts(m_grid, "velbar");
  {
    auto velbar = array::cast<array::Vector>(StressBalanceVelbar(model).compute(geometry));

    // copy_from communicates ghosts
    velbar_with_ghosts.copy_from(*velbar);
  }

  array::CellType1 cell_type(m_grid, "cell_type");
  cell_type.copy_from(geometry.cell_type);

  compute_2D_principal_strain_rates(velbar_with_ghosts, cell_type, *result);

  return result;
}

StressBalanceDeviatoricStresses::StressBalanceDeviatoricStresses(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "sigma_xx", *m_grid },
             { m_sys, "sigma_yy", *m_grid },
             { m_sys, "sigma_xy", *m_grid } };

  m_vars[0].long_name("deviatoric stress in X direction").units("Pa");
  m_vars[1].long_name("deviatoric stress in Y direction").units("Pa");
  m_vars[2].long_name("deviatoric shear stress").units("Pa");
}

std::shared_ptr<array::Array> StressBalanceDeviatoricStresses::compute_impl(const Geometry &geometry) const {

  auto result = std::make_shared<array::Array2D<stressbalance::DeviatoricStresses> >(
      m_grid, "deviatoric_stresses", array::WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];
  result->metadata(2) = m_vars[2];

  const array::Array3D &enthalpy = model->energy_balance_model()->enthalpy();
  const array::Scalar &thickness = geometry.ice_thickness;

  array::Scalar hardness(m_grid, "hardness");
  array::Vector1 velocity(m_grid, "velocity");

  averaged_hardness_vec(*model->stress_balance()->shallow()->flow_law(), thickness, enthalpy,
                        hardness);

  // copy_from updates ghosts
  velocity.copy_from(*array::cast<array::Vector>(StressBalanceVelbar(model).compute(geometry)));

  array::CellType1 cell_type(m_grid, "cell_type");
  cell_type.copy_from(geometry.cell_type);

  stressbalance::compute_2D_stresses(*model->stress_balance()->shallow()->flow_law(), velocity,
                                     hardness, cell_type, *result);

  return result;
}

StressBalanceFlux::StressBalanceFlux(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = { { m_sys, "uflux", *m_grid }, { m_sys, "vflux", *m_grid } };
  m_vars[0]
      .long_name("Vertically integrated horizontal flux of ice in the X direction")
      .units("m^2 s^-1")
      .output_units("m^2 year^-1");
  m_vars[1]
      .long_name("Vertically integrated horizontal flux of ice in the Y direction")
      .units("m^2 s^-1")
      .output_units("m^2 year^-1");
}

std::shared_ptr<array::Array> StressBalanceFlux::compute_impl(const Geometry &geometry) const {
  double H_threshold = m_config->get_number("geometry.ice_free_thickness_standard");

  auto result = allocate<array::Vector>("flux");

  // get the thickness
  const array::Scalar &thickness = geometry.ice_thickness;

  const array::Array3D
    &u3 = model->stress_balance()->velocity_u(),
    &v3 = model->stress_balance()->velocity_v();

  array::AccessScope list{&u3, &v3, &thickness, result.get()};

  const auto &z = m_grid->z();

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      double H = thickness(i,j);

      // an "ice-free" cell:
      if (H < H_threshold) {
        (*result)(i, j) = 0.0;
        continue;
      }

      // an icy cell:
      {
        const auto *u = u3.get_column(i, j);
        const auto *v = v3.get_column(i, j);

        Vector2d Q(0.0, 0.0);

        // ks is "k just below the surface"
        auto ks = m_grid->kBelowHeight(H);

        if (ks > 0) {
          Vector2d v0(u[0], v[0]);

          for (unsigned int k = 1; k <= ks; ++k) {
            Vector2d v1(u[k], v[k]);

            // trapezoid rule
            Q += (z[k] - z[k - 1]) * 0.5 * (v0 + v1);

            v0 = v1;
          }
        }

        // rectangle method to integrate over the last level
        Q += (H - z[ks]) * Vector2d(u[ks], v[ks]);

        (*result)(i, j) = Q;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


StressBalanceFluxMag::StressBalanceFluxMag(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = { { m_sys, "flux_mag", *m_grid } };
  m_vars[0]
      .long_name("magnitude of vertically-integrated horizontal flux of ice")
      .units("m^2 s^-1")
      .output_units("m^2 year^-1");
  m_vars[0]["_FillValue"] = { fill_value() };
  m_vars[0]["valid_min"]  = { 0.0 };
}

std::shared_ptr<array::Array> StressBalanceFluxMag::compute_impl(const Geometry &geometry) const {
  const array::Scalar &thickness = geometry.ice_thickness;

  // Compute the vertically-averaged horizontal ice velocity:
  auto result = array::cast<array::Scalar>(StressBalanceVelbarMag(model).compute(geometry));

  result->metadata() = m_vars[0];

  array::AccessScope list{&thickness, result.get()};

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) *= thickness(i,j);
  }

  apply_mask(thickness, fill_value(), *result);

  return result;
}

StressBalanceVelbar::StressBalanceVelbar(const IceModel *m)
  : Diag<IceModel>(m) {

  auto ismip = m_config->get_flag("output.ISMIP");

  // set metadata:
  m_vars = {{m_sys, ismip ? "xvelmean" : "ubar", *m_grid },
            {m_sys, ismip ? "yvelmean" : "vbar", *m_grid }};

  m_vars[0]
      .long_name("vertical mean of horizontal ice velocity in the X direction")
      .standard_name("land_ice_vertical_mean_x_velocity")
      .units("m s^-1")
      .output_units("m year^-1");
  m_vars[1]
      .long_name("vertical mean of horizontal ice velocity in the Y direction")
      .standard_name("land_ice_vertical_mean_y_velocity")
      .units("m s^-1")
      .output_units("m year^-1");
}

std::shared_ptr<array::Array> StressBalanceVelbar::compute_impl(const Geometry &geometry) const {
  // get the thickness
  const array::Scalar& thickness = geometry.ice_thickness;

  // Compute the vertically-integrated horizontal ice flux:
  auto result = array::cast<array::Vector>(StressBalanceFlux(model).compute(geometry));

  // Override metadata set by the flux computation
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  array::AccessScope list{&thickness, result.get()};

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();
    double thk = thickness(i,j);

    // Ice flux is masked already, but we need to check for division
    // by zero anyway.
    if (thk > 0.0) {
      (*result)(i,j) /= thk;
    } else {
      (*result)(i,j) = 0.0;
    }
  }

  return result;
}

StressBalanceVelbarMag::StressBalanceVelbarMag(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars = {{m_sys, "velbar_mag", *m_grid }};

  m_vars[0]
      .long_name("magnitude of vertically-integrated horizontal velocity of ice")
      .units("m second^-1")
      .output_units("m year^-1");

  m_vars[0]["_FillValue"] = { fill_value() };
  m_vars[0]["valid_min"]  = { 0.0 };
}

std::shared_ptr<array::Array> StressBalanceVelbarMag::compute_impl(const Geometry &geometry) const {
  auto result = allocate<array::Scalar>("velbar_mag");

  // compute vertically-averaged horizontal velocity:
  auto velbar = array::cast<array::Vector>(StressBalanceVelbar(model).compute(geometry));

  // compute its magnitude:
  compute_magnitude(*velbar, *result);

  // mask out ice-free areas:
  apply_mask(geometry.ice_thickness, fill_value(), *result);

  return result;
}

StressBalanceVonmisesStress::StressBalanceVonmisesStress(const IceModel *m)
    : Diag<IceModel>(m) {
  m_vars = { { m_sys, "vonmises_stress", *m_grid } };
  m_vars[0].long_name("tensile von Mises stress").units("Pascal");
}

std::shared_ptr<array::Array> StressBalanceVonmisesStress::compute_impl(const Geometry &geometry) const {

  using std::max;
  using std::sqrt;
  using std::pow;

  auto result = allocate<array::Scalar>("vonmises_stress");

  array::Scalar &vonmises_stress = *result;

  auto velbar = array::cast<array::Vector>(StressBalanceVelbar(model).compute(geometry));
  array::Vector &velocity = *velbar;

  using StrainRates = array::Array2D<stressbalance::PrincipalStrainRates>;
  auto eigen12 = array::cast<StrainRates>(StressBalanceStrainRates(model).compute(geometry));
  const auto &strain_rates = *eigen12;

  const array::Scalar &ice_thickness = geometry.ice_thickness;
  const array::Array3D &enthalpy = model->energy_balance_model()->enthalpy();
  const auto &mask = geometry.cell_type;

  std::shared_ptr<const rheology::FlowLaw> flow_law;

  if (m_config->get_flag("calving.vonmises_calving.use_custom_flow_law")) {
    rheology::FlowLawFactory factory(m_config, m_grid->ctx()->enthalpy_converter());
    flow_law = factory.create(m_config->get_string("calving.vonmises_calving.flow_law"),
                              m_config->get_number("calving.vonmises_calving.Glen_exponent"));
  } else {
    flow_law = model->stress_balance()->shallow()->flow_law();
  }

  double glen_exponent = flow_law->exponent();

  array::AccessScope list{&vonmises_stress, &velocity, &strain_rates, &ice_thickness,
      &enthalpy, &mask};

  for (auto pt : m_grid->points()) {
    const int i = pt.i(), j = pt.j();

    if (mask.icy(i, j)) {

      const double       H = ice_thickness(i, j);
      const unsigned int k = m_grid->kBelowHeight(H);

      const double
        *enthalpy_column   = enthalpy.get_column(i, j),
        hardness           = averaged_hardness(*flow_law, H, k, m_grid->z().data(), enthalpy_column),
        eigen1             = strain_rates(i, j).eigen1,
        eigen2             = strain_rates(i, j).eigen2;

      // [\ref Morlighem2016] equation 6
      const double effective_tensile_strain_rate = sqrt(0.5 * (pow(max(0.0, eigen1), 2) +
                                                               pow(max(0.0, eigen2), 2)));
      // [\ref Morlighem2016] equation 7
      vonmises_stress(i, j) = sqrt(3.0) * hardness * pow(effective_tensile_strain_rate,
                                                         1.0 / glen_exponent);

    } else { // end of "if (mask.icy(i, j))"
      vonmises_stress(i, j) = 0.0;
    }
  }   // end of loop over grid points

  return result;
}

StressBalanceTauxz::StressBalanceTauxz(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = {{m_sys, "tauxz", *m_grid, m_grid->z()}};
  m_vars[0].long_name("shear stress xz component (in shallow ice approximation SIA)").units("Pa");
}


/*!
 * The SIA-applicable shear stress component tauxz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH StressBalanceTauyz
 */
std::shared_ptr<array::Array> StressBalanceTauxz::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "tauxz", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata() = m_vars[0];

  const array::Scalar &thickness = geometry.ice_thickness,
                      &surface   = geometry.ice_surface_elevation;

  array::AccessScope list{ &surface, &thickness, result.get() };

  const double rg = m_config->get_number("constants.ice.density") *
                    m_config->get_number("constants.standard_gravity");

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      unsigned int ks      = m_grid->kBelowHeight(thickness(i, j));
      double *tauxz_out_ij = result->get_column(i, j);
      const double H = thickness(i, j), dhdx = diff_x_p(surface, i, j);

      // within the ice:
      for (unsigned int k = 0; k <= ks; ++k) {
        tauxz_out_ij[k] = -rg * (H - m_grid->z(k)) * dhdx;
      }
      // above the ice:
      for (unsigned int k = ks + 1; k < m_grid->Mz(); ++k) {
        tauxz_out_ij[k] = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


StressBalanceTauyz::StressBalanceTauyz(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "tauyz", *m_grid, m_grid->z() } };
  m_vars[0].long_name("shear stress yz component (in shallow ice approximation SIA)").units("Pa");
}


/*!
 * The SIA-applicable shear stress component tauyz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH StressBalanceTauxz
 */
std::shared_ptr<array::Array> StressBalanceTauyz::compute_impl(const Geometry &geometry) const {

  std::shared_ptr<array::Array3D> result(
      new array::Array3D(m_grid, "tauyz", array::WITHOUT_GHOSTS, m_grid->z()));
  result->metadata(0) = m_vars[0];

  const array::Scalar &thickness = geometry.ice_thickness,
                      &surface   = geometry.ice_surface_elevation;

  array::AccessScope list{ &surface, &thickness, result.get() };

  const double rg = m_config->get_number("constants.ice.density") *
                    m_config->get_number("constants.standard_gravity");

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      unsigned int ks      = m_grid->kBelowHeight(thickness(i, j));
      double *tauyz_out_ij = result->get_column(i, j);
      const double H = thickness(i, j), dhdy = diff_y_p(surface, i, j);

      // within the ice:
      for (unsigned int k = 0; k <= ks; ++k) {
        tauyz_out_ij[k] = -rg * (H - m_grid->z(k)) * dhdy;
      }
      // above the ice:
      for (unsigned int k = ks + 1; k < m_grid->Mz(); ++k) {
        tauyz_out_ij[k] = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

StressBalanceVelbaseMag::StressBalanceVelbaseMag(const IceModel *m)
  : Diag<IceModel>(m) {

  m_vars = { { m_sys, "velbase_mag", *m_grid } };
  m_vars[0]
      .long_name("magnitude of horizontal velocity of ice at base of ice")
      .units("m s^-1")
      .output_units("m year^-1");
  m_vars[0]["_FillValue"] = { fill_value() };
  m_vars[0]["valid_min"]  = { 0.0 };
}

std::shared_ptr<array::Array> StressBalanceVelbaseMag::compute_impl(const Geometry &geometry) const {
  auto result = allocate<array::Scalar>("velbase_mag");

  compute_magnitude(*array::cast<array::Vector>(StressBalanceVelbase(model).compute(geometry)),
                    *result);

  double fill = fill_value();

  const auto &mask = geometry.cell_type;

  array::AccessScope list{&mask, result.get()};

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j) = fill;
    }
  }

  return result;
}

StressBalanceVelsurfMag::StressBalanceVelsurfMag(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars = { { m_sys, "velsurf_mag", *m_grid } };
  m_vars[0]
      .long_name("magnitude of horizontal velocity of ice at ice surface")
      .units("m s^-1")
      .output_units("m year^-1");
  m_vars[0]["_FillValue"] = { fill_value() };
  m_vars[0]["valid_min"] = {0.0};
}

std::shared_ptr<array::Array> StressBalanceVelsurfMag::compute_impl(const Geometry &geometry) const {
  double fill = fill_value();

  auto result = allocate<array::Scalar>("velsurf_mag");

  compute_magnitude(*array::cast<array::Vector>(StressBalanceVelsurf(model).compute(geometry)),
                    *result);

  const auto &mask = geometry.cell_type;

  array::AccessScope list{&mask, result.get()};

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j) = fill;
    }
  }

  return result;
}


StressBalanceVelsurf::StressBalanceVelsurf(const IceModel *m)
  : Diag<IceModel>(m) {

  auto ismip = m_config->get_flag("output.ISMIP");
  m_vars = { { m_sys, ismip ? "xvelsurf" : "uvelsurf", *m_grid },
             { m_sys, ismip ? "yvelsurf" : "vvelsurf", *m_grid } };
  m_vars[0]
      .long_name("x-component of the horizontal velocity of ice at ice surface")
      .standard_name("land_ice_surface_x_velocity"); // InitMIP "standard" name
  m_vars[1]
      .long_name("y-component of the horizontal velocity of ice at ice surface")
      .standard_name("land_ice_surface_y_velocity"); // InitMIP "standard" name

  for (auto &v : m_vars) {
    v.units("m s^-1").output_units("m year^-1");
    v["_FillValue"] = { fill_value() };
  }
}

std::shared_ptr<array::Array> StressBalanceVelsurf::compute_impl(const Geometry &geometry) const {
  double fill = fill_value();

  auto result = allocate<array::Vector>("surf");

  array::Scalar u_surf(m_grid, "u_surf");
  array::Scalar v_surf(m_grid, "v_surf");

  const array::Array3D
    &u3 = model->stress_balance()->velocity_u(),
    &v3 = model->stress_balance()->velocity_v();

  const auto &thickness = geometry.ice_thickness;
  const auto &cell_type = geometry.cell_type;

  extract_surface(u3, thickness, u_surf);
  extract_surface(v3, thickness, v_surf);

  array::AccessScope list{ &cell_type, &u_surf, &v_surf, result.get() };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free(i, j)) {
      (*result)(i, j) = fill;
    } else {
      (*result)(i, j) = { u_surf(i, j), v_surf(i, j) };
    }
  }

  return result;
}

StressBalanceWvel::StressBalanceWvel(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "wvel", *m_grid } };
  m_vars[0]
      .long_name("vertical velocity of ice, relative to geoid")
      .units("m s^-1")
      .output_units("m year^-1");
}

std::shared_ptr<array::Array> StressBalanceWvel::compute(const Geometry &geometry,
                                                         bool zero_above_ice) const {
  std::shared_ptr<array::Array3D> result3(
      new array::Array3D(m_grid, "wvel", array::WITHOUT_GHOSTS, m_grid->z()));
  result3->metadata() = m_vars[0];

  const auto &bed = geometry.bed_elevation;
  const auto &uplift = model->bed_deformation_model()->uplift();

  const auto &thickness = geometry.ice_thickness;
  const auto &mask      = geometry.cell_type;

  const array::Array3D &u3 = model->stress_balance()->velocity_u(),
                       &v3 = model->stress_balance()->velocity_v(),
                       &w3 = model->vertical_velocity();

  array::AccessScope list{ &thickness, &mask, &bed, &u3, &v3, &w3, &uplift, result3.get() };

  const double ice_density       = m_config->get_number("constants.ice.density"),
               sea_water_density = m_config->get_number("constants.sea_water.density"),
               R                 = ice_density / sea_water_density;

  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      const double *u = u3.get_column(i, j), *v = v3.get_column(i, j), *w = w3.get_column(i, j);
      double *result = result3->get_column(i, j);

      int ks = m_grid->kBelowHeight(thickness(i, j));

      // in the ice:
      if (mask.grounded(i, j)) {
        const double bed_dx = diff_x_p(bed, i, j), bed_dy = diff_y_p(bed, i, j),
                     uplift_ij = uplift(i, j);
        for (int k = 0; k <= ks; k++) {
          result[k] = w[k] + uplift_ij + u[k] * bed_dx + v[k] * bed_dy;
        }

      } else { // floating
        const double z_sl = R * thickness(i, j), w_sl = w3.interpolate(i, j, z_sl);

        for (int k = 0; k <= ks; k++) {
          result[k] = w[k] - w_sl;
        }
      }

      // above the ice:
      if (zero_above_ice) {
        // set to zeros
        for (unsigned int k = ks + 1; k < m_grid->Mz(); k++) {
          result[k] = 0.0;
        }
      } else {
        // extrapolate using the topmost value
        for (unsigned int k = ks + 1; k < m_grid->Mz(); k++) {
          result[k] = result[ks];
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result3;
}

std::shared_ptr<array::Array> StressBalanceWvel::compute_impl(const Geometry &geometry) const {
  return this->compute(geometry, true); // fill wvel above the ice with zeros
}

StressBalanceWvelsurf::StressBalanceWvelsurf(const IceModel *m) : Diag<IceModel>(m) {

  auto ismip = m_config->get_flag("output.ISMIP");

  // set metadata:
  m_vars = { { m_sys, ismip ? "zvelsurf" : "wvelsurf", *m_grid } };

  m_vars[0]
      .long_name("vertical velocity of ice at ice surface, relative to the geoid")
      .standard_name("land_ice_surface_upward_velocity") // InitMIP "standard" name
      .units("m s^-1")
      .output_units("m year^-1");

  m_vars[0]["_FillValue"]  = { fill_value() };
}

std::shared_ptr<array::Array> StressBalanceWvelsurf::compute_impl(const Geometry &geometry) const {
  double fill = fill_value();

  auto result = allocate<array::Scalar>("wvelsurf");

  // here "false" means "don't fill w3 above the ice surface with zeros"
  auto w3 = array::cast<array::Array3D>(StressBalanceWvel(model).compute(geometry, false));

  const auto &thickness = geometry.ice_thickness;

  extract_surface(*w3, thickness, *result);

  const auto &mask = geometry.cell_type;

  array::AccessScope list{ &mask, result.get() };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j) = fill;
    }
  }

  return result;
}

StressBalanceWvelbase::StressBalanceWvelbase(const IceModel *m) : Diag<IceModel>(m) {

  auto ismip = m_config->get_flag("output.ISMIP");

  m_vars = { { m_sys, ismip ? "zvelbase" : "wvelbase", *m_grid } };
  m_vars[0]
      .long_name("vertical velocity of ice at the base of ice, relative to the geoid")
      .standard_name("land_ice_basal_upward_velocity") // InitMIP "standard" name
      .units("m s^-1")
      .output_units("m year^-1");

  m_vars[0]["_FillValue"]  = { fill_value() };
}

std::shared_ptr<array::Array> StressBalanceWvelbase::compute_impl(const Geometry &geometry) const {
  double fill = fill_value();

  auto result = allocate<array::Scalar>("wvelbase");

  // here "false" means "don't fill w3 above the ice surface with zeros"
  auto w3 = array::cast<array::Array3D>(StressBalanceWvel(model).compute(geometry, false));

  extract_surface(*w3, 0.0, *result);

  const auto &mask = geometry.cell_type;

  array::AccessScope list{ &mask, result.get() };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j) = fill;
    }
  }

  return result;
}

StressBalanceVelbase::StressBalanceVelbase(const IceModel *m) : Diag<IceModel>(m) {

  auto ismip = m_config->get_flag("output.ISMIP");

  // set metadata:
  m_vars = { { m_sys, ismip ? "xvelbase" : "uvelbase", *m_grid },
             { m_sys, ismip ? "yvelbase" : "vvelbase", *m_grid } };
  m_vars[0]
      .long_name("x-component of the horizontal velocity of ice at the base of ice")
      .standard_name("land_ice_basal_x_velocity"); // InitMIP "standard" name
  m_vars[1]
      .long_name("y-component of the horizontal velocity of ice at the base of ice")
      .standard_name("land_ice_basal_y_velocity"); // InitMIP "standard" name

  for (auto &v : m_vars) {
    v.units("m s^-1").output_units("m year^-1");
    v["_FillValue"]  = { fill_value() };
  }
}

std::shared_ptr<array::Array> StressBalanceVelbase::compute_impl(const Geometry &geometry) const {
  double fill = fill_value();

  auto result = allocate<array::Vector>("base");

  array::Scalar u_base(m_grid, "u_base");
  array::Scalar v_base(m_grid, "v_base");

  const array::Array3D &u3 = model->stress_balance()->velocity_u(),
                       &v3 = model->stress_balance()->velocity_v();

  extract_surface(u3, 0.0, u_base);
  extract_surface(v3, 0.0, v_base);

  const auto &mask = geometry.cell_type;

  array::AccessScope list{ &mask, &u_base, &v_base, result.get() };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j) = fill;
    } else {
      (*result)(i, j) = { u_base(i, j), v_base(i, j) };
    }
  }

  return result;
}

ShallowStressBalanceBeta::ShallowStressBalanceBeta(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "beta", *m_grid } };
  m_vars[0].long_name("basal drag coefficient").units("Pa s / m");
}

std::shared_ptr<array::Array> ShallowStressBalanceBeta::compute_impl(const Geometry &/*geometry*/) const {
  auto result = allocate<array::Scalar>("beta");

  const auto *yield_stress_model = model->basal_yield_stress_model();

  if (yield_stress_model == nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "basal yield stress model is not available");
  }

  const auto &tauc = yield_stress_model->basal_material_yield_stress();

  const auto *shallow_stress_balance = model->stress_balance()->shallow();

  const auto *basal_sliding_law = shallow_stress_balance->sliding_law();

  if (basal_sliding_law == nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "basal sliding law is not available");
  }

  const auto &velocity = shallow_stress_balance->velocity();

  array::AccessScope list{&tauc, &velocity, result.get()};
  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) =  basal_sliding_law->drag(tauc(i,j), velocity(i,j).u, velocity(i,j).v);
  }

  return result;
}

//! \brief Report the wall melt rate from dissipation of the potential energy of the
//! transportable water.
class HydrologyWallMelt : public Diag<IceModel>
{
public:
  HydrologyWallMelt(const IceModel *m)
    : Diag<IceModel>(m) {
    m_vars = { { m_sys, "wallmelt", *m_grid } };
    m_vars[0]
        .long_name("wall melt into subglacial hydrology layer from (turbulent)"
                   " dissipation of energy in transportable water")
        .units("m s^-1")
        .output_units("m year^-1");

    const auto *routing_hydrology =
        dynamic_cast<const hydrology::Routing *>(model->subglacial_hydrology_model());

    if (routing_hydrology == nullptr) {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION, "cannot compute 'wallmelt': routing hydrology is not available");
    }
  }

protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const {
    auto result = allocate<array::Scalar>("wallmelt");

    const array::Scalar &bed_elevation = geometry.bed_elevation;

    const auto *routing_hydrology =
        dynamic_cast<const hydrology::Routing *>(model->subglacial_hydrology_model());

    if (routing_hydrology == nullptr) {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION, "cannot compute 'wallmelt': routing hydrology is not available");
    }

    hydrology::wall_melt(*routing_hydrology, bed_elevation, *result);
    return result;
  }
};

/*! @brief Report hydraulic potential in the subglacial hydrology system */
class HydraulicPotential : public Diag<IceModel>
{
public:
  HydraulicPotential(const IceModel *m) : Diag<IceModel>(m) {
    m_vars = { { m_sys, "hydraulic_potential", *m_grid } };
    m_vars[0].long_name("hydraulic potential in the subglacial hydrology system").units("Pa");

    const auto *routing_hydrology =
        dynamic_cast<const hydrology::Routing *>(model->subglacial_hydrology_model());

    if (routing_hydrology == nullptr) {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "cannot compute 'hydraulic_potential': routing hydrology is not available");
    }
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const {

    auto result = allocate<array::Scalar>("hydraulic_potential");

    const auto &sea_level     = geometry.sea_level_elevation;
    const auto &bed_elevation = geometry.bed_elevation;
    const auto &ice_thickness = geometry.ice_thickness;

    const auto *routing_hydrology =
        dynamic_cast<const hydrology::Routing *>(model->subglacial_hydrology_model());

    if (routing_hydrology == nullptr) {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "cannot compute 'hydraulic_potential': routing hydrology is not available");
    }

    hydrology::hydraulic_potential(routing_hydrology->subglacial_water_thickness(),
                                   routing_hydrology->subglacial_water_pressure(), sea_level,
                                   bed_elevation, ice_thickness, *result);

    return result;
  }
};

class HeatFluxFromBedrock : public DiagAverageRate<IceModel>
{
public:
  HeatFluxFromBedrock(const IceModel *m)
      : DiagAverageRate<IceModel>(m, "heat_flux_from_bedrock", RATE) {

    m_accumulator.metadata()["units"] = "joule m^-2";

    auto ismip = m_config->get_flag("output.ISMIP");
    m_vars      = { { m_sys, ismip ? "hfgeoubed" : "heat_flux_from_bedrock", *m_grid } };
    m_vars[0]
        .long_name("upward heat flux from bedrock into ice at the ice base")
        .standard_name((ismip ? "upward_geothermal_heat_flux_in_land_ice" :
                                "upward_geothermal_heat_flux_at_ground_level_in_land_ice"))
        .units("W m^-2")
        .output_units("mW m^-2");
    m_vars[0]["cell_methods"] = "time: mean";
    m_vars[0]["comment"]      = "positive values correspond to an upward flux";
    m_vars[0]["_FillValue"]   = { fill_value() };
  }

protected:
  virtual void update_impl(double dt) {
    const auto &Q         = model->bedrock_thermal_model()->flux_through_top_surface();
    const auto &cell_type = model->geometry().cell_type;

    array::AccessScope list{&Q, &cell_type, &m_accumulator};

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      if (cell_type.grounded_ice(i, j)) {
        m_accumulator(i, j) += dt * Q(i, j);
      }
    }

    m_interval_length += dt;
  }
};

//! \brief Computes the multiplier \f$\theta\f$ in Schoof's (2003) theory of the
//! effect of bed roughness on the diffusivity of the SIA.
/*!
  See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
*/
class SIAFDSchoofs_Theta : public Diag<IceModel>
{
public:
  SIAFDSchoofs_Theta(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

SIAFDSchoofs_Theta::SIAFDSchoofs_Theta(const IceModel *m) : Diag<IceModel>(m) {
  m_vars = { { m_sys, "schoofs_theta", *m_grid } };

  m_vars[0]
      .long_name("multiplier 'theta' in Schoof's (2003) theory of bed roughness in SIA")
      .units("1");
  m_vars[0]["valid_range"] = { 0.0, 1.0 };

  const auto *siafd =
    dynamic_cast<const stressbalance::SIAFD*>(model->stress_balance()->modifier());

  if (siafd == nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot compute 'schoofs_theta': SIAFD is not available");
  }
}

std::shared_ptr<array::Array> SIAFDSchoofs_Theta::compute_impl(const Geometry &geometry) const {
  auto result = allocate<array::Scalar>("schoofs_theta");

  const auto *siafd =
    dynamic_cast<const stressbalance::SIAFD*>(model->stress_balance()->modifier());

  if (siafd == nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot compute 'schoofs_theta': SIAFD is not available");
  }

  siafd->bed_smoother().theta(geometry.ice_surface_elevation, *result);

  return result;
}

//! \brief Computes the thickness relative to the smoothed bed elevation in
//! Schoof's (2003) theory of the effect of bed roughness on the SIA.
/*!
  See page \ref bedrough and reference [\ref Schoofbasaltopg2003].
*/
class SIAFDThksmooth : public Diag<IceModel>
{
public:
  SIAFDThksmooth(const IceModel *m);
protected:
  virtual std::shared_ptr<array::Array> compute_impl(const Geometry &geometry) const;
};

SIAFDThksmooth::SIAFDThksmooth(const IceModel *m)
  : Diag<IceModel>(m) {

  m_vars = { { m_sys, "thksmooth", *m_grid } };
  m_vars[0]
      .long_name(
          "thickness relative to smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA")
      .units("m");
}

std::shared_ptr<array::Array> SIAFDThksmooth::compute_impl(const Geometry &geometry) const {

  array::CellType2 cell_type(m_grid, "cell_type");
  cell_type.copy_from(geometry.cell_type);

  const auto *siafd =
    dynamic_cast<const stressbalance::SIAFD*>(model->stress_balance()->modifier());

  if (siafd == nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot compute 'schoofs_theta': SIAFD is not available");
  }

  auto result = allocate<array::Scalar>("thksmooth");
  siafd->bed_smoother().smoothed_thk(geometry.ice_surface_elevation,
                                     geometry.ice_thickness, cell_type, *result);
  return result;
}

} // end of namespace diagnostics

void IceModel::init_outputs(InputOptions options, DiagnosticReport report_type) {
  m_available_spatial_diagnostics = allocate_spatial_diagnostics();
  m_available_scalar_diagnostics = allocate_scalar_diagnostics();

  list_diagnostics(report_type);

  // initialize diagnostics first: we need to know which spatial and scalar variables will
  // be saved to determine which *state* variables for these diagnostics need to be saved
  // to output files that can be used for re-starting.
  init_scalar_diagnostics();
  init_spatial_diagnostics();
  // initialize outputs that can be used for re-starting:
  init_snapshots();
  init_checkpoints();
  init_final_output();

  deallocate_unused_diagnostics();

  // reset: this gives diagnostics a chance to capture the current state of the model at the
  // beginning of the run
  for (auto &d : m_available_spatial_diagnostics) {
    d.second->reset();
  }

  // read in the state (accumulators) if we are re-starting a run
  if (options.type == INIT_RESTART) {
    File file(m_grid->com, options.filename, io::PISM_GUESS, io::PISM_READONLY);
    for (const auto &d : m_available_spatial_diagnostics) {
      d.second->init(file, options.record);
    }
  }

  // Tell the output writer about all the variables we may need to write:
  {
    std::set<VariableMetadata> all_variables;
    all_variables = pism::combine(all_variables, m_output_file_contents);
    all_variables = pism::combine(all_variables, m_snapshot_file_contents);
    all_variables = pism::combine(all_variables, m_spatial_file_contents);
    all_variables = pism::combine(all_variables, m_checkpoint_file_contents);

    m_output_writer->initialize(all_variables);

    if (m_snapshot_writer != m_output_writer) {
      m_snapshot_writer->initialize(all_variables);
    }

    // note: no need to call m_spatial_writer->initialize() because it is equal to
    // m_snapshot_writer.
  }
}

std::map<std::string, Diagnostic::Ptr> IceModel::allocate_spatial_diagnostics() {
  std::map<std::string, Diagnostic::Ptr> result;

  using namespace diagnostics;

  using d       = Diagnostic;
  using f       = Diagnostic::Ptr; // "f" for "field"
  result        = {
           // geometry
    { "cell_grounded_fraction", d::wrap(m_geometry.cell_grounded_fraction) },
    { "height_above_flotation", f(new HeightAboveFloatation(this)) },
    { "ice_area_specific_volume", d::wrap(m_geometry.ice_area_specific_volume) },
    { "ice_mass", f(new IceMass(this)) },
    { "mask", d::wrap(m_geometry.cell_type) },
    { "pressure", f(new PressureInIce(this)) },
    { "thk", f(new IceThickness(this)) },
    { "topg_sl_adjusted", f(new BedTopographySeaLevelAdjusted(this)) },
    { "usurf", f(new IceSurfaceElevation(this)) },
    { "ice_base_elevation", f(new IceBottomSurfaceElevation(this)) },
    { floating_ice_sheet_area_fraction_name, f(new IceAreaFractionFloating(this)) },
    { grounded_ice_sheet_area_fraction_name, f(new IceAreaFractionGrounded(this)) },
    { land_ice_area_fraction_name, f(new IceAreaFraction(this)) },

    // temperature, enthalpy, and liquid water fraction
    { "enthalpybase", f(new IceEnthalpyBasal(this)) },
    { "enthalpysurf", f(new IceEnthalpySurface(this)) },
    { "bedtoptemp", d::wrap(m_bedtoptemp) },
    { "cts", f(new CTS(this)) },
    { "liqfrac", f(new LiquidFraction(this)) },
    { "temp", f(new Temperature(this)) },
    { "temp_pa", f(new TemperaturePA(this)) },
    { "tempbase", f(new TemperatureBasal(this, BOTH)) },
    { "temppabase", f(new TemperaturePABasal(this)) },
    { "tempsurf", f(new TemperatureSurface(this)) },
    { "heat_flux_from_bedrock", f(new HeatFluxFromBedrock(this)) },

    // rheology-related stuff
    { "tempicethk", f(new TemperateIceThickness(this)) },
    { "tempicethk_basal", f(new TemperateIceThicknessBasal(this)) },
    { "hardav", f(new HardnessAverage(this)) },
    { "hardness", f(new IceHardness(this)) },
    { "effective_viscosity", f(new IceViscosity(this)) },

    // boundary conditions
    { "vel_bc_mask", d::wrap(m_velocity_bc_mask) },
    { "vel_bc_values", d::wrap(m_velocity_bc_values) },
    { "ice_margin_pressure_difference", f(new IceMarginPressureDifference(this)) },
    { "thk_bc_mask", d::wrap(m_ice_thickness_bc_mask) },

    // balancing the books
    // tendency_of_ice_amount = (tendency_of_ice_amount_due_to_flow +
    //                           tendency_of_ice_amount_due_to_conservation_error +
    //                           tendency_of_ice_amount_due_to_surface_mass_balance +
    //                           tendency_of_ice_amount_due_to_basal_mass_balance +
    //                           tendency_of_ice_amount_due_to_discharge)
    //
    // Also,
    // tendency_of_ice_amount_due_to_discharge = (tendency_of_ice_amount_due_to_calving +
    //                                            tendency_of_ice_amount_due_to_frontal_melt +
    //                                            tendency_of_ice_amount_due_to_forced_retreat)
    { "tendency_of_ice_amount", f(new TendencyOfIceAmount(this, AMOUNT)) },
    { "tendency_of_ice_amount_due_to_flow", f(new TendencyOfIceAmountDueToFlow(this, AMOUNT)) },
    { "tendency_of_ice_amount_due_to_conservation_error",
             f(new ConservationErrorFlux(this, AMOUNT)) },
    { "tendency_of_ice_amount_due_to_surface_mass_flux", f(new SurfaceFlux(this, AMOUNT)) },
    { "tendency_of_ice_amount_due_to_basal_mass_flux", f(new BasalFlux(this, AMOUNT)) },
    { "tendency_of_ice_amount_due_to_discharge", f(new DischargeFlux(this, AMOUNT)) },
    { "tendency_of_ice_amount_due_to_calving", f(new CalvingFlux(this, AMOUNT)) },
    { "tendency_of_ice_amount_due_to_frontal_melt", f(new FrontalMeltFlux(this, AMOUNT)) },
    { "tendency_of_ice_amount_due_to_forced_retreat", f(new ForcedRetreatFlux(this, AMOUNT)) },

    // same, in terms of mass
    // tendency_of_ice_mass = (tendency_of_ice_mass_due_to_flow +
    //                         tendency_of_ice_mass_due_to_conservation_error +
    //                         tendency_of_ice_mass_due_to_surface_mass_flux +
    //                         tendency_of_ice_mass_due_to_basal_mass_balance +
    //                         tendency_of_ice_mass_due_to_discharge)
    //
    // Also,
    // tendency_of_ice_mass_due_to_discharge = (tendency_of_ice_mass_due_to_calving +
    //                                          tendency_of_ice_mass_due_to_frontal_melt +
    //                                          tendency_of_ice_mass_due_to_forced_retreat)
    { "tendency_of_ice_mass", f(new TendencyOfIceAmount(this, MASS)) },
    { "tendency_of_ice_mass_due_to_flow", f(new TendencyOfIceAmountDueToFlow(this, MASS)) },
    { "tendency_of_ice_mass_due_to_conservation_error", f(new ConservationErrorFlux(this, MASS)) },
    { "tendency_of_ice_mass_due_to_surface_mass_flux", f(new SurfaceFlux(this, MASS)) },
    { "tendency_of_ice_mass_due_to_basal_mass_flux", f(new BasalFlux(this, MASS)) },
    { "tendency_of_ice_mass_due_to_discharge", f(new DischargeFlux(this, MASS)) },
    { "tendency_of_ice_mass_due_to_calving", f(new CalvingFlux(this, MASS)) },
    { "tendency_of_ice_mass_due_to_frontal_melt", f(new FrontalMeltFlux(this, MASS)) },
    { "tendency_of_ice_mass_due_to_forced_retreat", f(new ForcedRetreatFlux(this, MASS)) },

    // other rates and fluxes
    { "basal_mass_flux_grounded", f(new BMBSplit(this, GROUNDED)) },
    { "basal_mass_flux_floating", f(new BMBSplit(this, SHELF)) },
    { "dHdt", f(new ThicknessRateOfChange(this)) },
    { "bmelt", d::wrap(m_basal_melt_rate) },
    { "grounding_line_flux", f(new GroundingLineFlux(this)) },
    { "ice_mass_transport_across_grounding_line", f(new MassTransportAcrossGroundingLine(this)) },
    { "taud", f(new DrivingShearStress(this)) },
    { "taud_mag", f(new DrivingShearStressmagnitude(this)) },

    // velocities and fluxes that use ice geometry
    { "velbar_mag", Diagnostic::Ptr(new StressBalanceVelbarMag(this)) },
    { "flux", Diagnostic::Ptr(new StressBalanceFlux(this)) },
    { "flux_mag", Diagnostic::Ptr(new StressBalanceFluxMag(this)) },
    { "velbar", Diagnostic::Ptr(new StressBalanceVelbar(this)) },
    { "strain_rates", Diagnostic::Ptr(new StressBalanceStrainRates(this)) },
    { "deviatoric_stresses", Diagnostic::Ptr(new StressBalanceDeviatoricStresses(this)) },
    { "vonmises_stress", Diagnostic::Ptr(new StressBalanceVonmisesStress(this)) },
    { "uvel", Diagnostic::Ptr(new StressBalanceUvel(this)) },
    { "vvel", Diagnostic::Ptr(new StressBalanceVvel(this)) },
    { "wvel_rel", Diagnostic::Ptr(new StressBalanceWvelRel(this)) },
    { "tauxz", Diagnostic::Ptr(new StressBalanceTauxz(this)) },
    { "tauyz", Diagnostic::Ptr(new StressBalanceTauyz(this)) },
    { "velbase_mag", Diagnostic::Ptr(new StressBalanceVelbaseMag(this)) },
    { "velsurf_mag", Diagnostic::Ptr(new StressBalanceVelsurfMag(this)) },
    { "velbase", Diagnostic::Ptr(new StressBalanceVelbase(this)) },
    { "velsurf", Diagnostic::Ptr(new StressBalanceVelsurf(this)) },
    { "wvel", Diagnostic::Ptr(new StressBalanceWvel(this)) },
    { "wvelbase", Diagnostic::Ptr(new StressBalanceWvelbase(this)) },
    { "wvelsurf", Diagnostic::Ptr(new StressBalanceWvelsurf(this)) },
    { "bfrict", Diagnostic::Ptr(new StressBalanceBfrict(this)) },
    { "strainheat", Diagnostic::Ptr(new StressBalanceStrainheat(this)) },

    // misc
    { "rank", f(new Rank(this)) },
  };

  if (dynamic_cast<const stressbalance::SIAFD *>(stress_balance()->modifier()) != nullptr) {
    result["schoofs_theta"] = f(new SIAFDSchoofs_Theta(this));
    result["thksmooth"] = f(new SIAFDThksmooth(this));
  }

  if (stress_balance()->shallow()->sliding_law() != nullptr and
      basal_yield_stress_model() != nullptr) {
    result["beta"] = f(new ShallowStressBalanceBeta(this));
  }

  if (dynamic_cast<const hydrology::Routing *>(subglacial_hydrology_model()) != nullptr) {
    result["wallmelt"] = f(new HydrologyWallMelt(this));
    result["hydraulic_potential"] = f(new HydraulicPotential(this));
  }

  if (m_grid->has_longitude_latitude()) {
    result["lon"] = d::wrap(m_grid->longitude());
    result["lat"] = d::wrap(m_grid->latitude());
  }

  if (basal_yield_stress_model() != nullptr) {
    result["taub"]     = f(new BasalShearStress(this));
    result["taub_mag"] = f(new BasalShearStressMagnitude(this));
  }

#if (Pism_USE_PROJ == 1)
  std::string proj = m_grid->get_mapping_info()["proj_params"];
  if (not proj.empty()) {
    result["lat_bnds"] = f(new LatLonBounds(this, "lat", proj));
    result["lon_bnds"] = f(new LatLonBounds(this, "lon", proj));
  }
#endif

  // get diagnostics from submodels (may override some diagnostics allocated above)
  for (const auto& m : m_submodels) {
    result = pism::combine(result, m.second->spatial_diagnostics());
  }

  // add ISMIP variable names
  if (m_config->get_flag("output.ISMIP")) {
    result["acabf"]       = result["tendency_of_ice_amount_due_to_surface_mass_flux"];
    result["base"]        = result["ice_base_elevation"];
    result["dlithkdt"]    = result["dHdt"];
    result["hfgeoubed"]   = result["heat_flux_from_bedrock"];
    result["libmassbffl"] = result["basal_mass_flux_floating"];
    result["libmassbfgr"] = result["basal_mass_flux_grounded"];
    result["licalvf"]     = result["tendency_of_ice_amount_due_to_calving"];
    result["lifmassbf"]   = result["tendency_of_ice_amount_due_to_discharge"];
    result["ligroundf"]   = result["grounding_line_flux"];
    result["litempbotfl"] = f(new TemperatureBasal(this, SHELF));
    result["litempbotgr"] = f(new TemperatureBasal(this, GROUNDED));
    result["lithk"]       = result["thk"];
    result["orog"]        = result["usurf"];
    result["strbasemag"]  = result["taub_mag"];
    result["velmean"]     = result["velbar"];
    result["zvelbase"]    = result["wvelbase"];
    result["zvelsurf"]    = result["wvelsurf"];
  }

  return result;
}

std::map<std::string, TSDiagnostic::Ptr> IceModel::allocate_scalar_diagnostics() {
  std::map<std::string, TSDiagnostic::Ptr> result;

  using namespace diagnostics;

  using s = TSDiagnostic::Ptr; // "s" for "scalar"

  result = {
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
    {"max_horizontal_vel", s(new scalar::MaxHorizontalVelocity(this))},
    {"dt",              s(new scalar::TimeStepLength(this))},
    {"dt_ratio",        s(new scalar::TimeStepRatio(this))},
    // balancing the books
    {"tendency_of_ice_mass",                           s(new scalar::IceMassRateOfChange(this))},
    {"tendency_of_ice_mass_due_to_flow",               s(new scalar::IceMassRateOfChangeDueToFlow(this))},
    {"tendency_of_ice_mass_due_to_conservation_error", s(new scalar::IceMassFluxConservationError(this))},
    {"tendency_of_ice_mass_due_to_basal_mass_flux",    s(new scalar::IceMassFluxBasal(this))},
    {"tendency_of_ice_mass_due_to_surface_mass_flux",  s(new scalar::IceMassFluxSurface(this))},
    {"tendency_of_ice_mass_due_to_discharge",          s(new scalar::IceMassFluxDischarge(this))},
    {"tendency_of_ice_mass_due_to_calving",            s(new scalar::IceMassFluxCalving(this))},
    // other fluxes
    {"basal_mass_flux_grounded", s(new scalar::IceMassFluxBasalGrounded(this))},
    {"basal_mass_flux_floating", s(new scalar::IceMassFluxBasalFloating(this))},
    {"grounding_line_flux",      s(new scalar::IceMassFluxAtGroundingLine(this))},
  };

  // add ISMIP variable names
  if (m_config->get_flag("output.ISMIP")) {
    result["iareafl"]         = result["ice_area_glacierized_floating"];
    result["iareagr"]         = result["ice_area_glacierized_grounded"];
    result["lim"]             = result["ice_mass"];
    result["tendacabf"]       = result["tendency_of_ice_mass_due_to_surface_mass_flux"];
    result["tendlibmassbf"]   = result["tendency_of_ice_mass_due_to_basal_mass_flux"];
    result["tendlibmassbffl"] = result["basal_mass_flux_floating"];
    result["tendlibmassbfgr"] = result["basal_mass_flux_grounded"];
    result["tendlicalvf"]     = result["tendency_of_ice_mass_due_to_calving"];
    result["tendlifmassbf"]   = result["tendency_of_ice_mass_due_to_discharge"];
    result["tendligroundf"]   = result["grounding_line_flux"];
  }

  // get diagnostics from submodels
  for (const auto& m : m_submodels) {
    result = pism::combine(result, m.second->scalar_diagnostics());
  }

  return result;
}

using Metadata = std::map<std::string, std::vector<VariableMetadata>>;

static void print_diagnostics(const Logger &log, const Metadata &list) {
  for (const auto &d : list) {
    const std::string &name = d.first;
    log.message(1, " Name: %s\n", name.c_str());

    for (const auto &v : d.second) {

      std::string
        var_name     = v.get_name(),
        units        = v["units"],
        output_units = v["output_units"],
        long_name    = v["long_name"],
        comment      = v["comment"];

      if (not output_units.empty()) {
        units = output_units;
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
        var_name      = variable.get_name(),
        units         = variable["units"],
        output_units  = variable["output_units"],
        long_name     = variable["long_name"],
        standard_name = variable["standard_name"],
        comment       = variable["comment"];

      if (not output_units.empty()) {
        units = output_units;
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
static Metadata spatial_diag_metadata(const std::map<std::string,Diagnostic::Ptr> &diags) {
  Metadata result;

  for (const auto& f : diags) {
    auto diag = f.second;

    for (unsigned int k = 0; k < diag->n_variables(); ++k) {
      result[f.first].push_back(diag->metadata(k));
    }
  }

  return result;
}

/*!
 * Return metadata of scalar diagnostics.
 */
static Metadata scalar_diag_metadata(const std::map<std::string,TSDiagnostic::Ptr> &diags) {
  Metadata result;

  for (const auto& d : diags) {
    // always one variable per diagnostic
    result[d.first] = {d.second->metadata()};
  }

  return result;
}

void IceModel::list_diagnostics(DiagnosticReport type) const {

  if (type == DIAG_JSON) {
    m_log->message(1, "{\n");

    m_log->message(1, "\"spatial\" :\n");
    print_diagnostics_json(*m_log, spatial_diag_metadata(m_available_spatial_diagnostics));

    m_log->message(1, ",\n");        // separator

    m_log->message(1, "\"scalar\" :\n");
    print_diagnostics_json(*m_log, scalar_diag_metadata(m_available_scalar_diagnostics));

    m_log->message(1, "}\n");

    return;
  }

  if (type == DIAG_ALL or type == DIAG_SPATIAL) {
    m_log->message(1, "\n");
    m_log->message(1, "======== Available 2D and 3D diagnostics ========\n");

    print_diagnostics(*m_log, spatial_diag_metadata(m_available_spatial_diagnostics));
  }

  if (type == DIAG_ALL or type == DIAG_SCALAR) {
    // scalar time-series
    m_log->message(1, "======== Available time-series ========\n");

    print_diagnostics(*m_log, scalar_diag_metadata(m_available_scalar_diagnostics));
  }
}

/*!
  Computes fraction of the base which is melted.

  Communication occurs here.
 */
double IceModel::compute_temperate_base_fraction(double total_ice_area) {

  auto EC = m_ctx->enthalpy_converter();

  auto cell_area = m_grid->cell_area();

  double result = 0.0, meltarea = 0.0;

  const array::Array3D &enthalpy = m_energy_model->enthalpy();

  array::AccessScope list{&enthalpy, &m_geometry.cell_type, &m_geometry.ice_thickness};
  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
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
  meltarea = units::convert(m_sys, meltarea, "m^2", "km^2");
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

  if (not m_age_model) {
    return result;  // leave now
  }

  const double a = m_grid->cell_area() * 1e-3 * 1e-3, // area unit (km^2)
    currtime = m_time->current(); // seconds

  const array::Array3D &ice_age = m_age_model->age();

  array::AccessScope list{&m_geometry.cell_type, &m_geometry.ice_thickness, &ice_age};

  const double one_year = units::convert(m_sys, 1.0, "year", "seconds");
  double original_ice_volume = 0.0;

  // compute local original volume
  ParallelSection loop(m_grid->com);
  try {
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      if (m_geometry.cell_type.icy(i, j)) {
        // accumulate volume of ice which is original
        const double *age = ice_age.get_column(i, j);
        auto ks = m_grid->kBelowHeight(m_geometry.ice_thickness(i,j));
        for (unsigned int k = 1; k <= ks; k++) {
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

static void warn_about_missing(const Logger &log,
                        const std::set<std::string> &vars,
                        const std::string &type,
                        const std::set<std::string> &available,
                        bool stop) {
  std::vector<std::string> missing;
  for (const auto &v : vars) {
    if (available.find(v) == available.end()) {
      missing.push_back(v);
    }
  }

  if (not missing.empty()) {
    size_t N = missing.size();
    const char *ending = N > 1 ? "s" : "";
    const char *verb   = N > 1 ? "are" : "is";
    if (stop) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "%s variable%s %s %s not available!\n"
                                    "Available variables:\n- %s",
                                    type.c_str(),
                                    ending,
                                    join(missing, ",").c_str(),
                                    verb,
                                    set_join(available, ",\n- ").c_str());
    }

    log.message(2,
                "\nWARNING: %s variable%s %s %s not available!\n\n",
                type.c_str(),
                ending,
                join(missing, ",").c_str(),
                verb);
  }
}

/*!
 * De-allocate diagnostics that were not requested.
 *
 * Checks viewers, -spatial_vars, -checkpoint, -save_vars, and regular output.
 *
 * FIXME: I need to make sure that these reporting mechanisms are active. It is possible that
 * variables are on a list, but that list is not actually used.
 */
void IceModel::deallocate_unused_diagnostics() {

  // get the list of available diagnostics
  std::set<std::string> available;
  for (const auto &d : m_available_spatial_diagnostics) {
    available.insert(d.first);
  }

  auto stop = m_config->get_flag("output.spatial.stop_missing");
  warn_about_missing(*m_log, m_spatial_vars, "diagnostic", available, stop);

  // get the list of requested diagnostics
  auto requested = set_split(m_config->get_string("output.runtime.viewer.variables"), ',');
  requested = combine(requested, m_output_vars);
  requested = combine(requested, m_snapshot_vars);
  requested = combine(requested, m_spatial_vars);
  requested = combine(requested, m_checkpoint_vars);

  // de-allocate diagnostics that were not requested
  for (const auto &v : available) {
    if (requested.find(v) == requested.end()) {
      m_available_spatial_diagnostics.erase(v);
    }
  }
}

/*!
 * Update diagnostics.
 *
 * This usually involves accumulating data needed to computed time-averaged quantities.
 *
 * Call this after deallocate_unused_diagnostics() to avoid unnecessary work.
 */
void IceModel::update_diagnostics(double t, double dt) {
  for (const auto &d : m_available_spatial_diagnostics) {
    d.second->update(dt);
  }

  for (const auto &d : m_available_scalar_diagnostics) {
    d.second->update(t - dt, t);
  }
}

//! Writes variables listed in variable_names to file.
void IceModel::write_diagnostics(const OutputFile &file,
                                 const std::set<std::string> &variable_names) const {
  for (const auto &variable : variable_names) {
    auto diag = m_available_spatial_diagnostics.find(variable);

    if (diag != m_available_spatial_diagnostics.end()) {
      diag->second->compute(m_geometry)->write(file);
    }
  }
}

std::set<VariableMetadata>
IceModel::diagnostic_variables(const std::set<std::string> &variable_names) const {
  std::set<VariableMetadata> result{};
  {
    for (const auto &var : variable_names) {
      auto diag = m_available_spatial_diagnostics.find(var);

      if (diag != m_available_spatial_diagnostics.end()) {
        const auto &D = diag->second;
        for (unsigned int k = 0; k < D->n_variables(); ++k) {
          result.insert(D->metadata(k));
        }
      }
    }
  }
  return result;
}

std::set<VariableMetadata>
IceModel::state_variables_for_diagnostics(const std::set<std::string> &variable_names) const {
  std::set<VariableMetadata> result{};
  {
    for (const auto &var : variable_names) {
      auto diag = m_available_spatial_diagnostics.find(var);

      if (diag != m_available_spatial_diagnostics.end()) {
        const auto &D = diag->second;
        result = pism::combine(result, D->state());
      }
    }
  }
  return result;
}

void IceModel::write_state_for_diagnostics(const OutputFile &file,
                                       const std::set<std::string> &variable_names) const {

  for (const auto &var : variable_names) {
    auto diag = m_available_spatial_diagnostics.find(var);

    if (diag != m_available_spatial_diagnostics.end()) {
      const auto &D = diag->second;
      D->write_state(file);
    }
  }
}

} // end of namespace pism
