/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include <gsl/gsl_math.h>       // GSL_NAN

#include "pism/coupler/OceanModel.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {
namespace ocean {
OceanModel::OceanModel(IceGrid::ConstPtr g)
  : Component(g), m_sea_level(0) {

  m_melange_back_pressure_fraction.create(m_grid,
                                          "melange_back_pressure_fraction", WITHOUT_GHOSTS);
  m_melange_back_pressure_fraction.set_attrs("diagnostic",
                                             "melange back pressure fraction",
                                             "1", "");
  m_melange_back_pressure_fraction.set(0.0);

  m_shelf_base_temperature.create(m_grid, "shelfbtemp", WITHOUT_GHOSTS);
  m_shelf_base_temperature.set_attrs("diagnostic",
                                     "ice temperature at the bottom of floating ice",
                                     "Kelvin", "");

  m_shelf_base_mass_flux.create(m_grid, "m_shelf_base_mass_flux", WITHOUT_GHOSTS);
  m_shelf_base_mass_flux.set_attrs("diagnostic",
                                   "shelf base mass flux",
                                   "kg m-2 s-1", "");
  m_shelf_base_mass_flux.metadata().set_string("glaciological_units", "kg m-2 year-1");
}

OceanModel::~OceanModel() {
  // empty
}

void OceanModel::init() {
  // every re-init restarts the clock
  m_t  = GSL_NAN;
  m_dt = GSL_NAN;

  this->init_impl();
}

void OceanModel::update(double t, double dt) {
  // stop if this time step was taken care of by an earlier call
  if ((fabs(t - m_t) < 1e-12) and (fabs(dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = t;
  m_dt = dt;

  this->update_impl(t, dt);
}

const IceModelVec2S& OceanModel::shelf_base_mass_flux() const {
  return m_shelf_base_mass_flux;
}

double OceanModel::sea_level_elevation() const {
  return m_sea_level;
}

const IceModelVec2S& OceanModel::shelf_base_temperature() const {
  return m_shelf_base_temperature;
}

const IceModelVec2S& OceanModel::melange_back_pressure_fraction() const {
  return m_melange_back_pressure_fraction;
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

IceModelVec::Ptr PO_sea_level::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "sea_level", WITHOUT_GHOSTS));
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

IceModelVec::Ptr PO_shelf_base_temperature::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "shelfbtemp", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->shelf_base_temperature());

  return result;
}

PO_shelf_base_mass_flux::PO_shelf_base_mass_flux(const OceanModel *m)
  : Diag<OceanModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "shelfbmassflux")};

  set_attrs("mass flux at the basal surface of ice shelves", "",
            "kg m-2 s-1", "kg m-2 s-1", 0);
}

IceModelVec::Ptr PO_shelf_base_mass_flux::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "shelfbmassflux", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->shelf_base_mass_flux());

  return result;
}

PO_melange_back_pressure_fraction::PO_melange_back_pressure_fraction(const OceanModel *m)
  : Diag<OceanModel>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "melange_back_pressure_fraction")};

  set_attrs("dimensionless pressure fraction at calving fronts due to presence of melange ", "",
            "1", "1", 0);
}

IceModelVec::Ptr PO_melange_back_pressure_fraction::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "melange_back_pressure_fraction", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->melange_back_pressure_fraction());

  return result;
}

} // end of namespace ocean
} // end of namespace pism
