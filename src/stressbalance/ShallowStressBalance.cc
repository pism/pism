// Copyright (C) 2010--2023 Constantine Khroulev and Ed Bueler
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

#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/stressbalance/SSB_diagnostics.hh"
#include "pism/util/Context.hh"
#include "pism/util/Vars.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace stressbalance {

ShallowStressBalance::ShallowStressBalance(std::shared_ptr<const Grid> g)
  : Component(g),
    m_basal_sliding_law(NULL),
    m_flow_law(NULL),
    m_EC(g->ctx()->enthalpy_converter()),
    m_velocity(m_grid, "bar"),
    m_basal_frictional_heating(m_grid, "bfrict"),
    m_e_factor(1.0)
{

  if (m_config->get_flag("basal_resistance.pseudo_plastic.enabled")) {
    m_basal_sliding_law = new IceBasalResistancePseudoPlasticLaw(*m_config);
  } else if (m_config->get_flag("basal_resistance.regularized_coulomb.enabled")) {
    m_basal_sliding_law = new IceBasalResistanceRegularizedLaw(*m_config);
  } else {
    m_basal_sliding_law = new IceBasalResistancePlasticLaw(*m_config);
  }

  m_velocity.metadata(0)
      .long_name("thickness-advective ice velocity (x-component)")
      .units("m s^-1");
  m_velocity.metadata(1)
      .long_name("thickness-advective ice velocity (y-component)")
      .units("m s^-1");

  m_basal_frictional_heating.metadata(0)
      .long_name("basal frictional heating")
      .units("W m^-2")
      .output_units("mW m^-2");
}

ShallowStressBalance::~ShallowStressBalance() {
  delete m_basal_sliding_law;
}

void ShallowStressBalance::init() {
  this->init_impl();
}

void ShallowStressBalance::init_impl() {
  // empty
}

std::string ShallowStressBalance::stdout_report() const {
  return "";
}

std::shared_ptr<const rheology::FlowLaw> ShallowStressBalance::flow_law() const {
  return m_flow_law;
}

double ShallowStressBalance::flow_enhancement_factor() const {
  return m_e_factor;
}

EnthalpyConverter::Ptr ShallowStressBalance::enthalpy_converter() const {
  return m_EC;
}

const IceBasalResistancePlasticLaw* ShallowStressBalance::sliding_law() const {
  return m_basal_sliding_law;
}

//! \brief Get the thickness-advective 2D velocity.
const array::Vector1& ShallowStressBalance::velocity() const {
  return m_velocity;
}

//! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
const array::Scalar& ShallowStressBalance::basal_frictional_heating() {
  return m_basal_frictional_heating;
}


DiagnosticList ShallowStressBalance::diagnostics_impl() const {
  DiagnosticList result = {
    {"beta",     Diagnostic::Ptr(new SSB_beta(this))},
    {"taub",     Diagnostic::Ptr(new SSB_taub(this))},
    {"taub_mag", Diagnostic::Ptr(new SSB_taub_mag(this))},
    {"taud",     Diagnostic::Ptr(new SSB_taud(this))},
    {"taud_mag", Diagnostic::Ptr(new SSB_taud_mag(this))}
  };

  if(m_config->get_flag("output.ISMIP6")) {
    result["strbasemag"] = Diagnostic::Ptr(new SSB_taub_mag(this));
  }

  return result;
}


ZeroSliding::ZeroSliding(std::shared_ptr<const Grid> g)
  : ShallowStressBalance(g) {

  // Use the SIA flow law.
  rheology::FlowLawFactory ice_factory("stress_balance.sia.", m_config, m_EC);
  m_flow_law = ice_factory.create();
}

//! \brief Update the trivial shallow stress balance object.
void ZeroSliding::update(const Inputs &inputs, bool full_update) {
  (void) inputs;

  if (full_update) {
    m_velocity.set(0.0);
    m_basal_frictional_heating.set(0.0);
  }
}

//! \brief Compute the basal frictional heating.
/*!
  Ice shelves have zero basal friction heating.

  \param[in] V *basal* sliding velocity
  \param[in] tauc basal yield stress
  \param[in] mask (used to determine if floating or grounded)
  \param[out] result
 */
void ShallowStressBalance::compute_basal_frictional_heating(const array::Vector &V,
                                                            const array::Scalar &tauc,
                                                            const array::CellType &mask,
                                                            array::Scalar &result) const {

  array::AccessScope list{&V, &result, &tauc, &mask};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i,j)) {
      result(i,j) = 0.0;
    } else {
      const double
        C = m_basal_sliding_law->drag(tauc(i,j), V(i,j).u, V(i,j).v),
        basal_stress_x = - C * V(i,j).u,
        basal_stress_y = - C * V(i,j).v;
      result(i,j) = - basal_stress_x * V(i,j).u - basal_stress_y * V(i,j).v;
    }
  }
}


SSB_taud::SSB_taud(const ShallowStressBalance *m)
  : Diag<ShallowStressBalance>(m) {

  // set metadata:
  m_vars = { { m_sys, "taud_x" }, { m_sys, "taud_y" } };
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
std::shared_ptr<array::Array> SSB_taud::compute_impl() const {

  auto result = allocate<array::Vector>("taud");

  const array::Scalar *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const array::Scalar *surface   = m_grid->variables().get_2d_scalar("surface_altitude");

  double standard_gravity = m_config->get_number("constants.standard_gravity"),
         ice_density      = m_config->get_number("constants.ice.density");

  array::AccessScope list{ surface, thickness, result.get() };

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double pressure = ice_density * standard_gravity * (*thickness)(i, j);
    if (pressure <= 0.0) {
      (*result)(i, j).u = 0.0;
      (*result)(i, j).v = 0.0;
    } else {
      (*result)(i, j).u = -pressure * diff_x_p(*surface, i, j);
      (*result)(i, j).v = -pressure * diff_y_p(*surface, i, j);
    }
  }

  return result;
}

SSB_taud_mag::SSB_taud_mag(const ShallowStressBalance *m) : Diag<ShallowStressBalance>(m) {
  m_vars = { { m_sys, "taud_mag" } };
  m_vars[0]
      .long_name("magnitude of the gravitational driving stress at the base of ice")
      .units("Pa");
  m_vars[0]["comment"] = "this field is purely diagnostic (not used by the model)";
}

std::shared_ptr<array::Array> SSB_taud_mag::compute_impl() const {
  auto result = allocate<array::Scalar>("taud_mag");
  auto taud = array::cast<array::Vector>(SSB_taud(model).compute());

  compute_magnitude(*taud, *result);

  return result;
}

SSB_taub::SSB_taub(const ShallowStressBalance *m) : Diag<ShallowStressBalance>(m) {
  m_vars = { { m_sys, "taub_x" }, { m_sys, "taub_y" } };

  m_vars[0].long_name("X-component of the shear stress at the base of ice");
  m_vars[1].long_name("Y-component of the shear stress at the base of ice");

  for (auto &v : m_vars) {
    v.units("Pa");
    v["comment"] = "this field is purely diagnostic (not used by the model)";
  }
}


std::shared_ptr<array::Array> SSB_taub::compute_impl() const {

  auto result = allocate<array::Vector>("taub");

  const auto &velocity = model->velocity();
  const auto *tauc     = m_grid->variables().get_2d_scalar("tauc");
  const auto &mask     = *m_grid->variables().get_2d_cell_type("mask");

  const IceBasalResistancePlasticLaw *basal_sliding_law = model->sliding_law();

  array::AccessScope list{ tauc, &velocity, &mask, result.get() };
  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.grounded_ice(i, j)) {
      double beta     = basal_sliding_law->drag((*tauc)(i, j), velocity(i, j).u, velocity(i, j).v);
      (*result)(i, j) = -beta * velocity(i, j);
    } else {
      (*result)(i, j) = 0.0;
    }
  }

  return result;
}

SSB_taub_mag::SSB_taub_mag(const ShallowStressBalance *m) : Diag<ShallowStressBalance>(m) {

  auto ismip6 = m_config->get_flag("output.ISMIP6");

  m_vars = { { m_sys, ismip6 ? "strbasemag" : "taub_mag" } };
  m_vars[0]
      .long_name("magnitude of the basal shear stress at the base of ice")
      .standard_name("land_ice_basal_drag") // ISMIP6 "standard" name
      .units("Pa");
  m_vars[0]["comment"] = "this field is purely diagnostic (not used by the model)";
}

std::shared_ptr<array::Array> SSB_taub_mag::compute_impl() const {
  auto result = allocate<array::Scalar>("taub_mag");

  std::shared_ptr<array::Vector> taub = array::cast<array::Vector>(SSB_taub(model).compute());

  compute_magnitude(*taub, *result);

  return result;
}

/**
 * Shallow stress balance class that reads `u` and `v` fields from a
 * file and holds them constant.
 *
 * The only use I can think of right now is testing.
 */
PrescribedSliding::PrescribedSliding(std::shared_ptr<const Grid> g) : ZeroSliding(g) {
  // empty
}

void PrescribedSliding::update(const Inputs &inputs, bool full_update) {
  (void)inputs;
  if (full_update) {
    m_basal_frictional_heating.set(0.0);
  }
}

void PrescribedSliding::init_impl() {
  ShallowStressBalance::init_impl();

  auto input_filename = m_config->get_string("stress_balance.prescribed_sliding.file");

  if (input_filename.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION, "stress_balance.prescribed_sliding.file is required.");
  }

  m_velocity.regrid(input_filename, io::Default::Nil());
}

SSB_beta::SSB_beta(const ShallowStressBalance *m) : Diag<ShallowStressBalance>(m) {
  m_vars = { { m_sys, "beta" } };
  m_vars[0].long_name("basal drag coefficient").units("Pa s / m");
}

std::shared_ptr<array::Array> SSB_beta::compute_impl() const {
  auto result = allocate<array::Scalar>("beta");

  const array::Scalar *tauc = m_grid->variables().get_2d_scalar("tauc");

  const IceBasalResistancePlasticLaw *basal_sliding_law = model->sliding_law();

  const array::Vector &velocity = model->velocity();

  array::AccessScope list{tauc, &velocity, result.get()};
  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) =  basal_sliding_law->drag((*tauc)(i,j), velocity(i,j).u, velocity(i,j).v);
  }

  return result;
}

} // end of namespace stressbalance
} // end of namespace pism
