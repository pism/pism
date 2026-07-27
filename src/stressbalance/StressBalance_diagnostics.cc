// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2022, 2023, 2024, 2025, 2026 Constantine Khroulev
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

#include <cmath>                // std::pow, std::sqrt
#include <algorithm>            // std::max


#include "pism/stressbalance/StressBalance_diagnostics.hh"
#include "pism/stressbalance/SSB_Modifier.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Config.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/array/CellType.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace stressbalance {

DiagnosticList StressBalance::spatial_diagnostics_impl() const {
  DiagnosticList result = {
    {"bfrict",              Diagnostic::Ptr(new PSB_bfrict(this))},
    {"strainheat",          Diagnostic::Ptr(new PSB_strainheat(this))},
  };

  // add diagnostics from the shallow stress balance and the "modifier"
  result = pism::combine(result, m_shallow_stress_balance->spatial_diagnostics());
  result = pism::combine(result, m_modifier->spatial_diagnostics());

  return result;
}

TSDiagnosticList StressBalance::scalar_diagnostics_impl() const {
  return pism::combine(m_shallow_stress_balance->scalar_diagnostics(),
                       m_modifier->scalar_diagnostics());
}



PSB_bfrict::PSB_bfrict(const StressBalance *m) : Diag<StressBalance>(m) {
  m_vars = { { m_sys, "bfrict", *m_grid } };
  m_vars[0].long_name("basal frictional heating").units("W m^-2");
}

std::shared_ptr<array::Array> PSB_bfrict::compute_impl() const {

  auto result = allocate<array::Scalar>("bfrict");

  result->copy_from(model->basal_frictional_heating());

  return result;
}



PSB_strainheat::PSB_strainheat(const StressBalance *m) : Diag<StressBalance>(m) {
  m_vars = { { m_sys, "strainheat", *m_grid, m_grid->z() } };
  m_vars[0]
      .long_name("rate of strain heating in ice (dissipation heating)")
      .units("W m^-3")
      .output_units("mW m^-3");
}

std::shared_ptr<array::Array> PSB_strainheat::compute_impl() const {
  auto result = std::make_shared<array::Array3D>(m_grid, "strainheat",
                                                 array::WITHOUT_GHOSTS, m_grid->z());

  result->metadata() = m_vars[0];

  result->copy_from(model->volumetric_strain_heating());

  return result;
}


} // end of namespace stressbalance
} // end of namespace pism
