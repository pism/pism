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
#include "base/energy/EnergyModel.hh"
#include "base/energy/utilities.hh"

namespace pism {

//! \file iMenergy.cc Methods of IceModel which address conservation of energy.
//! Common to enthalpy (polythermal) and temperature (cold-ice) methods.

//! Manage the solution of the energy equation, and related parallel communication.
void IceModel::energyStep() {

  const Profiling &profiling = m_ctx->profiling();

  energy::EnergyModelStats stats;

  // operator-splitting occurs here (ice and bedrock energy updates are split):
  //   tell BedThermalUnit* btu that we have an ice base temp; it will return
  //   the z=0 value of geothermal flux when called inside temperatureStep() or
  //   enthalpyStep()

  IceModelVec2S &ice_surface_temperature = m_work2d[0];
  IceModelVec2S &bedtoptemp              = m_work2d[1];
  IceModelVec2S &basal_enthalpy          = m_work2d[2];
  m_energy_model->enthalpy().getHorSlice(basal_enthalpy, 0.0);
  m_surface->ice_surface_temperature(ice_surface_temperature);
  bedrock_surface_temperature(m_ocean->sea_level_elevation(),
                              m_cell_type,
                              m_beddef->bed_elevation(),
                              m_ice_thickness,
                              basal_enthalpy,
                              ice_surface_temperature,
                              bedtoptemp);

  profiling.begin("BTU");
  m_btu->update(bedtoptemp, t_TempAge, dt_TempAge);
  profiling.end("BTU");

  energy::EnergyModelInputs inputs;
  {
    IceModelVec2S &ice_surface_liquid_water_fraction = m_work2d[1];
    IceModelVec2S &till_water_thickness              = m_work2d[2];
    IceModelVec2S &shelf_base_temperature            = m_work2d[3];

    m_surface->ice_surface_liquid_water_fraction(ice_surface_liquid_water_fraction);

    m_ocean->shelf_base_temperature(shelf_base_temperature);

    m_subglacial_hydrology->till_water_thickness(till_water_thickness);

    inputs.basal_frictional_heating = &m_stress_balance->basal_frictional_heating();
    inputs.basal_heat_flux          = &m_btu->flux_through_top_surface(); // bedrock thermal layer
    inputs.cell_type                = &m_cell_type;                       // geometry
    inputs.ice_thickness            = &m_ice_thickness;                   // geometry
    inputs.shelf_base_temp          = &shelf_base_temperature;            // ocean model
    inputs.surface_liquid_fraction  = &ice_surface_liquid_water_fraction; // surface model
    inputs.surface_temp             = &ice_surface_temperature;           // surface model
    inputs.till_water_thickness     = &till_water_thickness;              // hydrology model

    inputs.strain_heating3          = &m_stress_balance->volumetric_strain_heating();
    inputs.u3                       = &m_stress_balance->velocity_u();
    inputs.v3                       = &m_stress_balance->velocity_v();
    inputs.w3                       = &m_stress_balance->velocity_w();

    inputs.check();             // make sure all data members were set
  }

  m_energy_model->update(t_TempAge, dt_TempAge, inputs);

  m_stdout_flags = m_energy_model->stdout_flags() + m_stdout_flags;
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

  IceModelVec2S &shelf_base_mass_flux = m_work2d[0];
  assert(m_ocean != NULL);
  m_ocean->shelf_base_mass_flux(shelf_base_mass_flux);

  const bool sub_gl = (m_config->get_boolean("geometry.grounded_cell_fraction") and
                       m_config->get_boolean("energy.basal_melt.use_grounded_cell_fraction"));

  IceModelVec::AccessList list;

  if (sub_gl) {
    list.add(m_gl_mask);
  }

  double ice_density = m_config->get_double("constants.ice.density");

  const IceModelVec2S &M_grounded = m_energy_model->basal_melt_rate();

  list.add(m_cell_type);
  list.add(M_grounded);
  list.add(shelf_base_mass_flux);
  list.add(m_basal_melt_rate);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double lambda = 1.0;      // 1.0 corresponds to the grounded case
    // Note: here we convert shelf base mass flux from [kg m-2 s-1] to [m s-1]:
    const double
      M_shelf_base = shelf_base_mass_flux(i,j) / ice_density;

    // Use the fractional floatation mask to adjust the basal melt
    // rate near the grounding line:
    if (sub_gl) {
      lambda = m_gl_mask(i,j);
    } else if (m_cell_type.ocean(i,j)) {
      lambda = 0.0;
    }
    m_basal_melt_rate(i,j) = lambda * M_grounded(i, j) + (1.0 - lambda) * M_shelf_base;
  }
}

//! @brief Compute the temperature seen by the top of the bedrock thermal layer.
void bedrock_surface_temperature(double sea_level,
                                 const IceModelVec2CellType &cell_type,
                                 const IceModelVec2S &bed_topography,
                                 const IceModelVec2S &ice_thickness,
                                 const IceModelVec2S &basal_enthalpy,
                                 const IceModelVec2S &ice_surface_temperature,
                                 IceModelVec2S &result) {

  IceGrid::ConstPtr grid  = result.get_grid();
  Config::ConstPtr config = grid->ctx()->config();

  const double
    T0                     = config->get_double("constants.fresh_water.melting_point_temperature"),
    beta_CC_grad_sea_water = (config->get_double("constants.ice.beta_Clausius_Clapeyron") *
                              config->get_double("constants.sea_water.density") *
                              config->get_double("constants.standard_gravity")); // K m-1

  EnthalpyConverter::Ptr EC = grid->ctx()->enthalpy_converter();

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(bed_topography);
  list.add(ice_thickness);
  list.add(ice_surface_temperature);
  list.add(basal_enthalpy);
  list.add(result);
  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.grounded(i,j)) {
        if (cell_type.ice_free(i,j)) { // no ice: sees air temp
          result(i,j) = ice_surface_temperature(i,j);
        } else { // ice: sees temp of base of ice
          const double pressure = EC->pressure(ice_thickness(i,j));
          result(i,j) = EC->temperature(basal_enthalpy(i,j), pressure);
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
