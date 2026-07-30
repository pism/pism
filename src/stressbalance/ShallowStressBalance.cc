// Copyright (C) 2010--2023, 2025, 2026 Constantine Khroulev and Ed Bueler
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
#include "pism/util/Context.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/IO_Flags.hh"

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
      .units("m s^-1")
      .output_units("m year^-1");
  m_velocity.metadata(1)
      .long_name("thickness-advective ice velocity (y-component)")
      .units("m s^-1")
      .output_units("m year^-1");

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
void ShallowStressBalance::update(const Inputs &inputs, bool full_update) {
  this->update_impl(inputs, full_update);
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

std::shared_ptr<EnthalpyConverter> ShallowStressBalance::enthalpy_converter() const {
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


DiagnosticList ShallowStressBalance::spatial_diagnostics_impl() const {
  DiagnosticList result = {
  };

  return {};
}


ZeroSliding::ZeroSliding(std::shared_ptr<const Grid> g)
  : ShallowStressBalance(g) {

  rheology::FlowLawFactory ice_factory(m_config, m_EC);
  // Use the SIA flow law.
  m_flow_law = ice_factory.create(m_config->get_string("stress_balance.sia.flow_law"),
                                  m_config->get_number("stress_balance.sia.Glen_exponent"));
}

//! \brief Update the trivial shallow stress balance object.
void ZeroSliding::update_impl(const Inputs &inputs, bool full_update) {
  (void) inputs;

  if (full_update) {
    m_velocity.set(0.0);
    m_basal_frictional_heating.set(0.0);
  }
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

void PrescribedSliding::update_impl(const Inputs &inputs, bool full_update) {
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

} // end of namespace stressbalance
} // end of namespace pism
