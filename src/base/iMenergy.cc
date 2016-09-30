// Copyright (C) 2004-2011, 2013, 2014, 2015, 2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "iceModel.hh"

#include "base/energy/BedThermalUnit.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "earth/PISMBedDef.hh"
#include "enthalpyConverter.hh"
#include "base/util/Profiling.hh"
#include "base/util/pism_utilities.hh"

#include "base/hydrology/PISMHydrology.hh"
#include "base/stressbalance/PISMStressBalance.hh"

namespace pism {

static void check_input(const IceModelVec *ptr, const char *name) {
  if (ptr == NULL) {
    throw RuntimeError::formatted("energy balance model input %s was not provided", name);
  }
}
EnergyModelInputs::EnergyModelInputs() {
  basal_frictional_heating = NULL;
  basal_heat_flux          = NULL;
  cell_type                = NULL;
  ice_thickness            = NULL;
  surface_liquid_fraction  = NULL;
  shelf_base_temp          = NULL;
  surface_temp             = NULL;
  till_water_thickness     = NULL;

  strain_heating3          = NULL;
  u3                       = NULL;
  v3                       = NULL;
  w3                       = NULL;
}

void EnergyModelInputs::check() const {
  check_input(cell_type,                "cell_type");
  check_input(basal_frictional_heating, "basal_frictional_heating");
  check_input(basal_heat_flux,          "basal_heat_flux");
  check_input(ice_thickness,            "ice_thickness");
  check_input(surface_liquid_fraction,  "surface_liquid_fraction");
  check_input(shelf_base_temp,          "shelf_base_temp");
  check_input(surface_temp,             "surface_temp");
  check_input(till_water_thickness,     "till_water_thickness");

  check_input(strain_heating3, "strain_heating3");
  check_input(u3, "u3");
  check_input(v3, "v3");
  check_input(w3, "w3");
}

//! \file iMenergy.cc Methods of IceModel which address conservation of energy.
//! Common to enthalpy (polythermal) and temperature (cold-ice) methods.

//! Manage the solution of the energy equation, and related parallel communication.
/*!
  This method updates three fields:
  - IceModelVec3 m_ice_enthalpy
  - IceModelVec2 basal_melt_rate
  - IceModelVec2 vHmelt
  That is, energyStep() is in charge of calling other methods that actually update, and
  then it is in charge of doing the ghost communication as needed.  If
  energy.temperature_based == true, then energyStep() must also update this field
  - IceModelVec3 m_ice_temperature

  Normally calls the method enthalpyAndDrainageStep().  Calls temperatureStep() if
  energy.temperature_based == true.
*/
void IceModel::energyStep() {

  const Profiling &profiling = m_ctx->profiling();

  EnergyModelStats stats;

  // operator-splitting occurs here (ice and bedrock energy updates are split):
  //   tell BedThermalUnit* btu that we have an ice base temp; it will return
  //   the z=0 value of geothermal flux when called inside temperatureStep() or
  //   enthalpyAndDrainageStep()
  IceModelVec2S &bedtoptemp = m_work2d[0];
  get_bed_top_temp(bedtoptemp);

  profiling.begin("BTU");
  m_btu->update(bedtoptemp, t_TempAge, dt_TempAge);
  profiling.end("BTU");

  EnergyModelInputs inputs;
  {
    m_surface->ice_surface_temperature(m_ice_surface_temp);
    m_surface->ice_surface_liquid_water_fraction(m_liqfrac_surface);

    m_ocean->shelf_base_temperature(m_shelfbtemp);

    IceModelVec2S &till_water_thickness = m_work2d[0];
    m_subglacial_hydrology->till_water_thickness(till_water_thickness);

    inputs.basal_frictional_heating = &m_stress_balance->basal_frictional_heating();
    inputs.basal_heat_flux          = &m_btu->flux_through_top_surface(); // bedrock thermal layer
    inputs.cell_type                = &m_cell_type;                       // geometry
    inputs.ice_thickness            = &m_ice_thickness;                   // geometry
    inputs.shelf_base_temp          = &m_shelfbtemp;                      // ocean model
    inputs.surface_liquid_fraction  = &m_liqfrac_surface;                 // surface model
    inputs.surface_temp             = &m_ice_surface_temp;                // surface model
    inputs.till_water_thickness     = &till_water_thickness;              // hydrology model

    inputs.strain_heating3          = &m_stress_balance->volumetric_strain_heating();
    inputs.u3                       = &m_stress_balance->velocity_u();
    inputs.v3                       = &m_stress_balance->velocity_v();
    inputs.w3                       = &m_stress_balance->velocity_w();

    inputs.check();             // make sure all data members were set
  }

  if (m_config->get_boolean("energy.temperature_based")) {
    // new temperature values go in vTnew; also updates Hmelt:
    profiling.begin("temp step");
    temperatureStep(inputs, stats);
    profiling.end("temp step");

    m_work3d.update_ghosts(m_ice_temperature);

    // compute_enthalpy_cold() updates ghosts of m_ice_enthalpy using
    // update_ghosts(). Is not optimized because this
    // (energy.temperature_based) is a rare case.
    compute_enthalpy_cold(m_ice_temperature, m_ice_thickness, m_ice_enthalpy);

  } else {
    // new enthalpy values go in m_work3d

    profiling.begin("enth step");
    enthalpyAndDrainageStep(inputs, stats);
    profiling.end("enth step");

    m_work3d.update_ghosts(m_ice_enthalpy);

    stats.liquified_ice_volume = GlobalSum(m_grid->com, stats.liquified_ice_volume);
    if (stats.liquified_ice_volume > 0.0) {
      m_log->message(1,
                 "\n PISM WARNING: fully-liquified cells detected: volume liquified = %.3f km^3\n\n",
                 stats.liquified_ice_volume / 1.0e9);
    }
  }

  // always count CFL violations for sanity check (but can occur only if -skip N with N>1)
  m_CFL_violation_counter = countCFLViolations();

  stats.reduced_accuracy_counter = GlobalSum(m_grid->com, stats.reduced_accuracy_counter);
  if (stats.reduced_accuracy_counter > 0.0) { // count of when BOMBPROOF switches to lower accuracy
    const double bfsacrPRCNT = 100.0 * (stats.reduced_accuracy_counter / (m_grid->Mx() * m_grid->My()));
    const double BPSACR_REPORT_VERB2_PERCENT = 5.0; // only report if above 5%
    if (bfsacrPRCNT > BPSACR_REPORT_VERB2_PERCENT &&
        m_log->get_threshold() > 2) {
      char tempstr[50] = "";
      snprintf(tempstr,50, "  [BPsacr=%.4f%%] ", bfsacrPRCNT);
      m_stdout_flags = tempstr + m_stdout_flags;
    }
  }

  stats.bulge_counter = GlobalSum(m_grid->com, stats.bulge_counter);
  if (stats.bulge_counter > 0) {
    // count of when advection bulges are limited; frequently it is identically zero
    char tempstr[50] = "";
    snprintf(tempstr,50, " BULGE=%d ", stats.bulge_counter);
    m_stdout_flags = tempstr + m_stdout_flags;
  }
}

//! @brief Combine basal melt rate in grounded and floating areas.
/**
 * Grounded basal melt rate is computed as a part of the energy
 * (enthalpy or temperature) step; floating basal melt rate is
 * provided by an ocean model.
 *
 * This method updates IceModel::basal_melt_rate (in meters per second
 * ice-equivalent).
 *
 * The sub shelf mass flux provided by an ocean model is in [kg m-2
 * s-1], so we divide by the ice density to convert to [m second-1].
 */
void IceModel::combine_basal_melt_rate() {

  assert(m_ocean != NULL);
  m_ocean->shelf_base_mass_flux(m_shelfbmassflux);

  const bool sub_gl = (m_config->get_boolean("geometry.grounded_cell_fraction") and
                       m_config->get_boolean("energy.basal_melt.use_grounded_cell_fraction"));

  IceModelVec::AccessList list;

  if (sub_gl) {
    list.add(m_gl_mask);
  }

  double ice_density = m_config->get_double("constants.ice.density");

  list.add(m_cell_type);
  list.add(m_basal_melt_rate);
  list.add(m_shelfbmassflux);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double lambda = 1.0;      // 1.0 corresponds to the grounded case
    // Note: here we convert shelf base mass flux from [kg m-2 s-1] to [m s-1]:
    const double
      M_grounded   = m_basal_melt_rate(i,j),
      M_shelf_base = m_shelfbmassflux(i,j) / ice_density;

    // Use the fractional floatation mask to adjust the basal melt
    // rate near the grounding line:
    if (sub_gl) {
      lambda = m_gl_mask(i,j);
    } else if (m_cell_type.ocean(i,j)) {
      lambda = 0.0;
    }
    m_basal_melt_rate(i,j) = lambda * M_grounded + (1.0 - lambda) * M_shelf_base;
  }
}

//! \brief Extract from enthalpy field (m_ice_enthalpy) the temperature which the top of
//! the bedrock thermal layer will see.
void IceModel::get_bed_top_temp(IceModelVec2S &result) {
  double
    T0                     = m_config->get_double("constants.fresh_water.melting_point_temperature"),
    beta_CC_grad_sea_water = (m_config->get_double("constants.ice.beta_Clausius_Clapeyron") * m_config->get_double("constants.sea_water.density") *
                              m_config->get_double("constants.standard_gravity")); // K m-1

  // will need coupler fields in ice-free land and
  assert(m_surface != NULL);
  m_surface->ice_surface_temperature(m_ice_surface_temp);

  assert(m_ocean != NULL);
  double sea_level = m_ocean->sea_level_elevation();

  // start by grabbing 2D enthalpy field at z=0; converted to temperature if needed, below
  m_ice_enthalpy.getHorSlice(result, 0.0);

  const IceModelVec2S &bed_topography = m_beddef->bed_elevation();

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  IceModelVec::AccessList list;
  list.add(bed_topography);
  list.add(result);
  list.add(m_ice_thickness);
  list.add(m_cell_type);
  list.add(m_ice_surface_temp);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_cell_type.grounded(i,j)) {
        if (m_cell_type.ice_free(i,j)) { // no ice: sees air temp
          result(i,j) = m_ice_surface_temp(i,j);
        } else { // ice: sees temp of base of ice
          const double pressure = EC->pressure(m_ice_thickness(i,j));
          result(i,j) = EC->temperature(result(i,j), pressure);
        }
      } else { // floating: apply pressure melting temp as top of bedrock temp
        result(i,j) = T0 - (sea_level - bed_topography(i,j)) * beta_CC_grad_sea_water;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

} // end of namespace pism
