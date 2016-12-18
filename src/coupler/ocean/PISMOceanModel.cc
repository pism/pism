/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "coupler/PISMOcean.hh"
#include "base/util/iceModelVec.hh"

namespace pism {
namespace ocean {
OceanModel::OceanModel(IceGrid::ConstPtr g)
  : Component_TS(g), m_sea_level(0) {
  // empty
}

OceanModel::~OceanModel() {
  // empty
}

void OceanModel::init() {
  this->init_impl();
}

double OceanModel::sea_level_elevation() const {
  double result;
  this->sea_level_elevation_impl(result);
  return result;
}

void OceanModel::shelf_base_temperature(IceModelVec2S &result) const {
  this->shelf_base_temperature_impl(result);
}

void OceanModel::shelf_base_mass_flux(IceModelVec2S &result) const {
  this->shelf_base_mass_flux_impl(result);
}

void OceanModel::melange_back_pressure_fraction(IceModelVec2S &result) const {
  this->melange_back_pressure_fraction_impl(result);
}

/** Set `result` to the melange back pressure fraction.
 *
 * This default implementation sets `result` to 0.0.
 *
 * @param[out] result back pressure fraction
 *
 * @return 0 on success
 */
void OceanModel::melange_back_pressure_fraction_impl(IceModelVec2S &result) const {
  result.set(0.0);
}

std::map<std::string, Diagnostic::Ptr> OceanModel::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"sea_level",                      Diagnostic::Ptr(new PO_sea_level(this))},
    {"shelfbtemp",                     Diagnostic::Ptr(new PO_shelf_base_temperature(this))},
    {"shelfbmassflux",                 Diagnostic::Ptr(new PO_shelf_base_mass_flux(this))},
    {"melange_back_pressure_fraction", Diagnostic::Ptr(new PO_melange_back_pressure_fraction(this))}
  };
  return result;
}

PO_sea_level::PO_sea_level(const OceanModel *m)
  : Diag<OceanModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "sea_level")};

  set_attrs("sea level elevation, relative to the geoid", "",
            "meters", "meters", 0);
}

IceModelVec::Ptr PO_sea_level::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "sea_level", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  result->set(model->sea_level_elevation());

  return result;
}

PO_shelf_base_temperature::PO_shelf_base_temperature(const OceanModel *m)
  : Diag<OceanModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "shelfbtemp")};

  set_attrs("ice temperature at the basal surface of ice shelves", "",
            "Kelvin", "Kelvin", 0);
}

IceModelVec::Ptr PO_shelf_base_temperature::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "shelfbtemp", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->shelf_base_temperature(*result);

  return result;
}

PO_shelf_base_mass_flux::PO_shelf_base_mass_flux(const OceanModel *m)
  : Diag<OceanModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "shelfbmassflux")};

  set_attrs("mass flux at the basal surface of ice shelves", "",
            "kg m-2 s-1", "kg m-2 s-1", 0);
}

IceModelVec::Ptr PO_shelf_base_mass_flux::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "shelfbmassflux", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->shelf_base_mass_flux(*result);

  return result;
}

PO_melange_back_pressure_fraction::PO_melange_back_pressure_fraction(const OceanModel *m)
  : Diag<OceanModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "melange_back_pressure_fraction")};

  set_attrs("dimensionless pressure fraction at calving fronts due to presence of melange ", "",
            "1", "1", 0);
}

IceModelVec::Ptr PO_melange_back_pressure_fraction::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "melange_back_pressure_fraction", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->melange_back_pressure_fraction(*result);

  return result;
}

} // end of namespace ocean
} // end of namespace pism
