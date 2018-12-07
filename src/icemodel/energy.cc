// Copyright (C) 2004-2011, 2013, 2014, 2015, 2016, 2017, 2018 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "IceModel.hh"

#include "pism/energy/BedThermalUnit.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/Profiling.hh"

#include "pism/hydrology/Hydrology.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/energy/EnergyModel.hh"
#include "pism/energy/utilities.hh"

namespace pism {

//! \file energy.cc Methods of IceModel which address conservation of energy.
//! Common to enthalpy (polythermal) and temperature (cold-ice) methods.

void IceModel::bedrock_thermal_model_step() {

  const Profiling &profiling = m_ctx->profiling();

  IceModelVec2S &basal_enthalpy = m_work2d[2];

  m_energy_model->enthalpy().getHorSlice(basal_enthalpy, 0.0);

  bedrock_surface_temperature(m_geometry.sea_level_elevation,
                              m_geometry.cell_type,
                              m_geometry.bed_elevation,
                              m_geometry.ice_thickness,
                              basal_enthalpy,
                              m_surface->temperature(),
                              m_bedtoptemp);

  profiling.begin("btu");
  m_btu->update(m_bedtoptemp, t_TempAge, dt_TempAge);
  profiling.end("btu");
}

//! Manage the solution of the energy equation, and related parallel communication.
void IceModel::energy_step() {

  // operator-splitting occurs here (ice and bedrock energy updates are split):
  //   tell BedThermalUnit* btu that we have an ice base temp; it will return
  //   the z=0 value of geothermal flux when called inside temperatureStep() or
  //   enthalpyStep()
  bedrock_thermal_model_step();

  m_energy_model->update(t_TempAge, dt_TempAge, energy_model_inputs());

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
void IceModel::combine_basal_melt_rate(const Geometry &geometry,
                                       const IceModelVec2S &shelf_base_mass_flux,
                                       const IceModelVec2S &grounded_basal_melt_rate,
                                       IceModelVec2S &result) {

  const bool sub_gl = (m_config->get_boolean("geometry.grounded_cell_fraction") and
                       m_config->get_boolean("energy.basal_melt.use_grounded_cell_fraction"));

  IceModelVec::AccessList list{&geometry.cell_type,
      &grounded_basal_melt_rate, &shelf_base_mass_flux, &result};
  if (sub_gl) {
    list.add(geometry.cell_grounded_fraction);
  }

  double ice_density = m_config->get_double("constants.ice.density");

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double lambda = 1.0;      // 1.0 corresponds to the grounded case
    // Note: here we convert shelf base mass flux from [kg m-2 s-1] to [m s-1]:
    const double
      M_shelf_base = shelf_base_mass_flux(i,j) / ice_density;

    // Use the fractional floatation mask to adjust the basal melt
    // rate near the grounding line:
    if (sub_gl) {
      lambda = geometry.cell_grounded_fraction(i,j);
    } else if (geometry.cell_type.ocean(i,j)) {
      lambda = 0.0;
    }
    result(i,j) = lambda * grounded_basal_melt_rate(i, j) + (1.0 - lambda) * M_shelf_base;
  }
}

//! @brief Compute the temperature seen by the top of the bedrock thermal layer.
void bedrock_surface_temperature(const IceModelVec2S &sea_level,
                                 const IceModelVec2CellType &cell_type,
                                 const IceModelVec2S &bed_topography,
                                 const IceModelVec2S &ice_thickness,
                                 const IceModelVec2S &basal_enthalpy,
                                 const IceModelVec2S &ice_surface_temperature,
                                 IceModelVec2S &result) {

  IceGrid::ConstPtr grid  = result.grid();
  Config::ConstPtr config = grid->ctx()->config();

  const double
    T0                     = config->get_double("constants.fresh_water.melting_point_temperature"),
    beta_CC_grad_sea_water = (config->get_double("constants.ice.beta_Clausius_Clapeyron") *
                              config->get_double("constants.sea_water.density") *
                              config->get_double("constants.standard_gravity")); // K m-1

  EnthalpyConverter::Ptr EC = grid->ctx()->enthalpy_converter();

  IceModelVec::AccessList list{&cell_type, &bed_topography, &sea_level, &ice_thickness,
      &ice_surface_temperature, &basal_enthalpy, &result};
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
        result(i,j) = T0 - (sea_level(i, j) - bed_topography(i,j)) * beta_CC_grad_sea_water;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

} // end of namespace pism
