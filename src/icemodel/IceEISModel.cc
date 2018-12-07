// Copyright (C) 2004-2018 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "IceEISModel.hh"

#include "pism/util/Context.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"

#include "pism/coupler/ocean/Constant.hh"
#include "pism/coupler/ocean/Initialization.hh"

#include "pism/coupler/SeaLevel.hh"
#include "pism/coupler/ocean/sea_level/Initialization.hh"

#include "pism/coupler/surface/EISMINTII.hh"
#include "pism/coupler/surface/Initialization.hh"

#include "pism/earth/BedDef.hh"

namespace pism {

IceEISModel::IceEISModel(IceGrid::Ptr g, Context::Ptr context, char experiment)
  : IceModel(g, context), m_experiment(experiment) {

  // the following flag must be here in constructor because
  // IceModel::createVecs() uses it non-polythermal methods; can be
  // overridden by the command-line option "-energy enthalpy"
  m_config->set_boolean("energy.temperature_based", true);

  // see EISMINT II description; choose no ocean interaction,
  m_config->set_boolean("ocean.always_grounded", true);

  // purely SIA, and E=1
  m_config->set_double("stress_balance.sia.enhancement_factor", 1.0);

  // none use bed smoothing & bed roughness parameterization
  m_config->set_double("stress_balance.sia.bed_smoother.range", 0.0);

  // basal melt does not change computation of mass continuity or vertical velocity:
  m_config->set_boolean("geometry.update.use_basal_melt_rate", false);

  // Make bedrock thermal material properties into ice properties.  Note that
  // zero thickness bedrock layer is the default, but we want the ice/rock
  // interface segment to have geothermal flux applied directly to ice without
  // jump in material properties at base.
  m_config->set_double("energy.bedrock_thermal_density",
                       m_config->get_double("constants.ice.density"));
  m_config->set_double("energy.bedrock_thermal_conductivity",
                       m_config->get_double("constants.ice.thermal_conductivity"));
  m_config->set_double("energy.bedrock_thermal_specific_heat_capacity",
                       m_config->get_double("constants.ice.specific_heat_capacity"));

  // no sliding + SIA
  m_config->set_string("stress_balance.model", "sia");
}

void IceEISModel::allocate_couplers() {

  // Climate will always come from intercomparison formulas.
  if (not m_surface) {
    std::shared_ptr<surface::SurfaceModel> surface(new surface::EISMINTII(m_grid, m_experiment));
    m_surface.reset(new surface::InitializationHelper(m_grid, surface));
    m_submodels["surface process model"] = m_surface.get();
  }

  if (not m_ocean) {
    std::shared_ptr<ocean::OceanModel> ocean(new ocean::Constant(m_grid));
    m_ocean.reset(new ocean::InitializationHelper(m_grid, ocean));
    m_submodels["ocean model"] = m_ocean.get();
  }

  if (not m_sea_level) {
    using namespace ocean::sea_level;
    std::shared_ptr<SeaLevel> sea_level(new SeaLevel(m_grid));
    m_sea_level.reset(new InitializationHelper(m_grid, sea_level));
    m_submodels["sea level forcing"] = m_sea_level.get();
  }
}

void generate_trough_topography(IceModelVec2S &result) {
  // computation based on code by Tony Payne, 6 March 1997:
  // http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f

  IceGrid::ConstPtr grid = result.grid();

  const double
    b0    = 1000.0,  // plateau elevation
    L     = 750.0e3, // half-width of computational domain
    w     = 200.0e3, // trough width
    slope = b0 / L,
    dx61  = (2.0 * L) / 60; // = 25.0e3

  IceModelVec::AccessList list(result);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double nsd = i * grid->dx(), ewd = j * grid->dy();
    if ((nsd >= (27 - 1) * dx61) && (nsd <= (35 - 1) * dx61) &&
        (ewd >= (31 - 1) * dx61) && (ewd <= (61 - 1) * dx61)) {
      result(i,j) = 1000.0 - std::max(0.0, slope * (ewd - L) * cos(M_PI * (nsd - L) / w));
    } else {
      result(i,j) = 1000.0;
    }
  }
}

void generate_mound_topography(IceModelVec2S &result) {
  // computation based on code by Tony Payne, 6 March 1997:
  // http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f

  IceGrid::ConstPtr grid = result.grid();

  const double slope = 250.0;
  const double w     = 150.0e3; // mound width

  IceModelVec::AccessList list(result);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double nsd = i * grid->dx(), ewd = j * grid->dy();
    result(i,j) = fabs(slope * sin(M_PI * ewd / w) + slope * cos(M_PI * nsd / w));
  }
}

void IceEISModel::initialize_2d() {

  m_log->message(2,
                 "initializing variables from EISMINT II experiment %c formulas... \n",
                 m_experiment);

  // set bed topography
  if (m_experiment == 'I' or m_experiment == 'J') {
    generate_trough_topography(m_geometry.bed_elevation);
  } else if (m_experiment == 'K' or m_experiment == 'L') {
    generate_mound_topography(m_geometry.bed_elevation);
  } else {
    m_geometry.bed_elevation.set(0.0);
  }

  m_geometry.sea_level_elevation.set(0.0);

  // set uplift
  IceModelVec2S bed_uplift(m_grid, "uplift", WITHOUT_GHOSTS);
  bed_uplift.set(0.0);

  // start with zero ice
  m_geometry.ice_thickness.set(0.0);

  m_beddef->bootstrap(m_geometry.bed_elevation, bed_uplift, m_geometry.ice_thickness,
                      m_geometry.sea_level_elevation);
}

} // end of namespace pism
