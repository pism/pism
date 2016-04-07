// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <gsl/gsl_math.h>       // M_PI

#include "base/iceModel.hh"
#include "iceEISModel.hh"

#include "base/stressbalance/PISMStressBalance.hh"
#include "base/stressbalance/ShallowStressBalance.hh"
#include "base/stressbalance/sia/SIAFD.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Context.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "coupler/ocean/POConstant.hh"
#include "coupler/surface/PS_EISMINTII.hh"
#include "earth/PISMBedDef.hh"

namespace pism {

IceEISModel::IceEISModel(IceGrid::Ptr g, Context::Ptr context, char experiment)
  : IceModel(g, context), m_experiment(experiment) {
  m_config->set_string("EISMINT_II_experiment", std::string(1, m_experiment));
  m_config->set_string("EISMINT_II_experiment_doc", "EISMINT II experiment name");

  // the following flag must be here in constructor because
  // IceModel::createVecs() uses it non-polythermal methods; can be
  // overridden by the command-line option "-energy enthalpy"
  m_config->set_boolean("do_cold_ice_methods", true);

  // see EISMINT II description; choose no ocean interaction,
  m_config->set_boolean("is_dry_simulation", true);

  // purely SIA, and E=1
  m_config->set_double("sia_enhancement_factor", 1.0);

  // none use bed smoothing & bed roughness parameterization
  m_config->set_double("bed_smoother_range", 0.0);

  // basal melt does not change computation of mass continuity or vertical velocity:
  m_config->set_boolean("include_bmr_in_continuity", false);

  // Make bedrock thermal material properties into ice properties.  Note that
  // zero thickness bedrock layer is the default, but we want the ice/rock
  // interface segment to have geothermal flux applied directly to ice without
  // jump in material properties at base.
  m_config->set_double("bedrock_thermal_density",
                       m_config->get_double("ice_density"));
  m_config->set_double("bedrock_thermal_conductivity",
                       m_config->get_double("ice_thermal_conductivity"));
  m_config->set_double("bedrock_thermal_specific_heat_capacity",
                       m_config->get_double("ice_specific_heat_capacity"));
}

//! \brief Decide which stress balance model to use.
void IceEISModel::allocate_stressbalance() {

  using namespace pism::stressbalance;

  if (m_stress_balance != NULL) {
    return;
  }

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  ShallowStressBalance *my_stress_balance = new ZeroSliding(m_grid, EC);
  SSB_Modifier *modifier = new SIAFD(m_grid, EC);

  // ~StressBalance() will de-allocate my_stress_balance and modifier.
  m_stress_balance = new StressBalance(m_grid, my_stress_balance, modifier);

}

void IceEISModel::allocate_couplers() {

  // Climate will always come from intercomparison formulas.
  if (m_surface == NULL) {
    m_surface = new surface::EISMINTII(m_grid, m_experiment);
  }

  if (m_ocean == NULL) {
    m_ocean = new ocean::Constant(m_grid);
  }
}

void IceEISModel::generateTroughTopography(IceModelVec2S &result) {
  // computation based on code by Tony Payne, 6 March 1997:
  // http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f

  const double b0    = 1000.0;  // plateau elevation
  const double L     = 750.0e3; // half-width of computational domain
  const double w     = 200.0e3; // trough width
  const double slope = b0/L;
  const double dx61  = (2*L) / 60; // = 25.0e3

  IceModelVec::AccessList list(result);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double nsd = i * m_grid->dx(), ewd = j * m_grid->dy();
    if ((nsd >= (27 - 1) * dx61) && (nsd <= (35 - 1) * dx61) &&
        (ewd >= (31 - 1) * dx61) && (ewd <= (61 - 1) * dx61)) {
      result(i,j) = 1000.0 - std::max(0.0, slope * (ewd - L) * cos(M_PI * (nsd - L) / w));
    } else {
      result(i,j) = 1000.0;
    }
  }
}


void IceEISModel::generateMoundTopography(IceModelVec2S &result) {
  // computation based on code by Tony Payne, 6 March 1997:
  // http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f

  const double slope = 250.0;
  const double w     = 150.0e3; // mound width

  IceModelVec::AccessList list(result);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double nsd = i * m_grid->dx(), ewd = j * m_grid->dy();
    result(i,j) = fabs(slope * sin(M_PI * ewd / w) + slope * cos(M_PI * nsd / w));
  }
}


//! Only executed if NOT initialized from file (-i).
void IceEISModel::set_vars_from_options() {

  // initialize from EISMINT II formulas
  m_log->message(2,
             "initializing variables from EISMINT II experiment %c formulas... \n",
             m_experiment);

  IceModelVec2S tmp;
  tmp.create(m_grid, "topg", WITHOUT_GHOSTS);

  // set bed topography
  {
    if (m_experiment == 'I' or m_experiment == 'J') {
      generateTroughTopography(tmp);
    } else if (m_experiment == 'K' or m_experiment == 'L') {
      generateMoundTopography(tmp);
    } else {
      tmp.set(0.0);
    }

    m_beddef->set_elevation(tmp);
  }

  // set bed uplift; no experiments have uplift at start
  {
    tmp.set(0.0);
    m_beddef->set_uplift(tmp);
  }

  m_basal_melt_rate.set(0.0);
  m_geothermal_flux.set(0.042); // EISMINT II value; J m-2 s-1
  m_ice_thickness.set(0.0); // start with zero ice

  // regrid 2D variables
  regrid(2);

  // this IceModel bootstrap method should do right thing because of
  // variable settings above and init of coupler above
  putTempAtDepth();

  // regrid 3D variables
  regrid(3);
}


} // end of namespace pism
