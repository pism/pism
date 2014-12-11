// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "IceGrid.hh"
#include "iceModel.hh"
#include "iceEISModel.hh"
#include "SIAFD.hh"
#include "SIA_Sliding.hh"
#include "PISMStressBalance.hh"
#include "pism_options.hh"
#include "POConstant.hh"
#include "PS_EISMINTII.hh"

#include "error_handling.hh"

namespace pism {

IceEISModel::IceEISModel(IceGrid &g, Config &conf, Config &conf_overrides)
  : IceModel(g, conf, conf_overrides) {
  m_experiment = 'A';

  // the following flag must be here in constructor because
  // IceModel::createVecs() uses it non-polythermal methods; can be
  // overridden by the command-line option "-energy enthalpy"
  config.set_flag("do_cold_ice_methods", true);

  // see EISMINT II description; choose no ocean interaction, 
  config.set_flag("is_dry_simulation", true);

  // purely SIA, and E=1
  config.set_double("sia_enhancement_factor", 1.0);

  // none use bed smoothing & bed roughness parameterization
  config.set_double("bed_smoother_range", 0.0);

  // basal melt does not change computation of mass continuity or vertical velocity:
  config.set_flag("include_bmr_in_continuity", false);

  // Make bedrock thermal material properties into ice properties.  Note that
  // zero thickness bedrock layer is the default, but we want the ice/rock
  // interface segment to have geothermal flux applied directly to ice without
  // jump in material properties at base.
  config.set_double("bedrock_thermal_density", config.get("ice_density"));
  config.set_double("bedrock_thermal_conductivity", config.get("ice_thermal_conductivity"));
  config.set_double("bedrock_thermal_specific_heat_capacity", config.get("ice_specific_heat_capacity"));
}

void IceEISModel::set_grid_defaults() {
  double Lx = 750e3;
  grid.set_extent(0.0, 0.0, Lx, Lx);

  grid.time->init();
}

void IceEISModel::setFromOptions() {

  // set experiment name using command-line options
  {
    std::string name = "A";
    char temp = m_experiment;
    bool EISIIchosen;
    OptionsString("-eisII", "EISMINT II experiment name",
                  name, EISIIchosen);

    if (EISIIchosen == true) {
      temp = (char)toupper(name.c_str()[0]);
      if ((temp >= 'A') && (temp <= 'L')) {
        m_experiment = temp;
      } else {
        throw RuntimeError::formatted("option -eisII must have value A, B, C, D, E, F, G, H, I, J, K, or L; got %c",
                                      temp);
      }
    }

    char tempstr[2] = {temp, 0};
    config.set_string("EISMINT_II_experiment", tempstr);
  }

  IceModel::setFromOptions();
}


//! \brief Decide which stress balance model to use.
void IceEISModel::allocate_stressbalance() {

  if (stress_balance == NULL) {
    ShallowStressBalance *my_stress_balance;

    SSB_Modifier *modifier = new SIAFD(grid, *EC, config);

    if (m_experiment == 'G' || m_experiment == 'H') {
      my_stress_balance = new SIA_Sliding(grid, *EC, config);
    } else {
      my_stress_balance = new ZeroSliding(grid, *EC, config);
    }
  
    // ~StressBalance() will de-allocate my_stress_balance and modifier.
    stress_balance = new StressBalance(grid, my_stress_balance,
                                           modifier, config);

    // Note that in PISM stress balance computations are diagnostic, i.e. do not
    // have a state that changes in time. This means that this call can be here
    // and not in model_state_setup() and we don't need to re-initialize after
    // the "diagnostic time step".
    stress_balance->init(variables);

    if (config.get_flag("include_bmr_in_continuity")) {
      stress_balance->set_basal_melt_rate(&basal_melt_rate);
    }
  }
  
}

void IceEISModel::allocate_couplers() {

  // Climate will always come from intercomparison formulas.
  if (surface == NULL) {
    surface = new PS_EISMINTII(grid, config, m_experiment);
  }

  if (ocean == NULL) {
    ocean = new POConstant(grid, config);
  }
}

void IceEISModel::generateTroughTopography() {
  // computation based on code by Tony Payne, 6 March 1997:
  // http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f
  
  const double b0    = 1000.0;  // plateau elevation
  const double L     = 750.0e3; // half-width of computational domain
  const double w     = 200.0e3; // trough width
  const double slope = b0/L;
  const double dx61  = (2*L) / 60; // = 25.0e3

  IceModelVec::AccessList list;
  list.add(bed_topography);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double nsd = i * grid.dx(), ewd = j * grid.dy();
    if ((nsd >= (27 - 1) * dx61) && (nsd <= (35 - 1) * dx61) &&
        (ewd >= (31 - 1) * dx61) && (ewd <= (61 - 1) * dx61)) {
      bed_topography(i,j) = 1000.0 - std::max(0.0, slope * (ewd - L) * cos(M_PI * (nsd - L) / w));
    } else {
      bed_topography(i,j) = 1000.0;
    }
  }
}


void IceEISModel::generateMoundTopography() {
  // computation based on code by Tony Payne, 6 March 1997:
  // http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f
  
  const double slope = 250.0;
  const double w     = 150.0e3; // mound width

  IceModelVec::AccessList list;
  list.add(bed_topography);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double nsd = i * grid.dx(), ewd = j * grid.dy();
    bed_topography(i,j) = fabs(slope * sin(M_PI * ewd / w) + slope * cos(M_PI * nsd / w));
  }
}


//! Only executed if NOT initialized from file (-i).
void IceEISModel::set_vars_from_options() {

  // initialize from EISMINT II formulas
  verbPrintf(2, grid.com,
             "initializing variables from EISMINT II experiment %c formulas... \n", 
             m_experiment);

  if ((m_experiment == 'I') || (m_experiment == 'J')) {
    generateTroughTopography();
  } 
  if ((m_experiment == 'K') || (m_experiment == 'L')) {
    generateMoundTopography();
  } 

  // communicate b in any case; it will be horizontally-differentiated
  bed_topography.update_ghosts();

  basal_melt_rate.set(0.0); 
  geothermal_flux.set(0.042); // EISMINT II value; J m-2 s-1
  bed_uplift_rate.set(0.0); // no experiments have uplift at start
  ice_thickness.set(0.0); // start with zero ice

  // regrid 2D variables
  regrid(2);
  
  // this IceModel bootstrap method should do right thing because of
  // variable settings above and init of coupler above
  putTempAtDepth();

  // regrid 3D variables
  regrid(3);
}


} // end of namespace pism
