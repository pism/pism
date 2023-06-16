// Copyright (C) 2004-2018, 2021, 2022, 2023 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "pism/util/IceGrid.hh"

#include "pism/coupler/ocean/Constant.hh"
#include "pism/coupler/ocean/Initialization.hh"

#include "pism/coupler/SeaLevel.hh"
#include "pism/coupler/ocean/sea_level/Initialization.hh"

#include "pism/coupler/surface/EISMINTII.hh"
#include "pism/coupler/surface/Initialization.hh"

#include "pism/earth/BedDef.hh"

namespace pism {

IceEISModel::IceEISModel(std::shared_ptr<IceGrid> g,
                         std::shared_ptr<Context> context,
                         char experiment)
  : IceModel(g, context), m_experiment(experiment) {
}

void IceEISModel::allocate_couplers() {

  // Climate will always come from intercomparison formulas.
  if (not m_surface) {
    auto surface = std::make_shared<surface::EISMINTII>(m_grid, m_experiment);
    m_surface.reset(new surface::InitializationHelper(m_grid, surface));
    m_submodels["surface process model"] = m_surface.get();
  }

  if (not m_ocean) {
    auto ocean = std::make_shared<ocean::Constant>(m_grid);
    m_ocean.reset(new ocean::InitializationHelper(m_grid, ocean));
    m_submodels["ocean model"] = m_ocean.get();
  }

  if (not m_sea_level) {
    auto sea_level = std::make_shared<ocean::sea_level::SeaLevel>(m_grid);
    m_sea_level.reset(new ocean::sea_level::InitializationHelper(m_grid, sea_level));
    m_submodels["sea level forcing"] = m_sea_level.get();
  }
}

void generate_trough_topography(array::Scalar &result) {
  // computation based on code by Tony Payne, 6 March 1997:
  // http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f

  auto grid = result.grid();

  const double
    b0    = 1000.0,  // plateau elevation
    L     = 750.0e3, // half-width of computational domain
    w     = 200.0e3, // trough width
    slope = b0 / L,
    dx61  = (2.0 * L) / 60; // = 25.0e3

  array::AccessScope list(result);
  for (auto p = grid->points(); p; p.next()) {
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

void generate_mound_topography(array::Scalar &result) {
  // computation based on code by Tony Payne, 6 March 1997:
  // http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f

  auto grid = result.grid();

  const double slope = 250.0;
  const double w     = 150.0e3; // mound width

  array::AccessScope list(result);
  for (auto p = grid->points(); p; p.next()) {
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
  switch (m_experiment) {
  case 'I':
  case 'J':
    generate_trough_topography(m_geometry.bed_elevation);
    break;
  case 'K':
  case 'L':
    generate_mound_topography(m_geometry.bed_elevation);
    break;
  default:
    m_geometry.bed_elevation.set(0.0);
    break;
  }

  m_geometry.sea_level_elevation.set(0.0);

  // set uplift
  array::Scalar bed_uplift(m_grid, "uplift");
  bed_uplift.set(0.0);

  // start with zero ice
  m_geometry.ice_thickness.set(0.0);

  m_beddef->bootstrap(m_geometry.bed_elevation, bed_uplift, m_geometry.ice_thickness,
                      m_geometry.sea_level_elevation);
}

void IceEISModel::bootstrap_2d(const File &input_file) {
  (void) input_file;
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "EISMINT II mode does not support bootstrapping");
}


} // end of namespace pism
