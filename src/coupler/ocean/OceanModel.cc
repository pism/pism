/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021 PISM Authors
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
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {

namespace ocean {

IceModelVec2S::Ptr OceanModel::allocate_shelf_base_temperature(IceGrid::ConstPtr g) {
  IceModelVec2S::Ptr result(new IceModelVec2S(g, "shelfbtemp", WITHOUT_GHOSTS));
  result->set_attrs("diagnostic",
                    "ice temperature at the bottom of floating ice",
                    "Kelvin", "Kelvin", "", 0);
  return result;
}

IceModelVec2S::Ptr OceanModel::allocate_shelf_base_mass_flux(IceGrid::ConstPtr g) {
  IceModelVec2S::Ptr result(new IceModelVec2S(g, "shelfbmassflux", WITHOUT_GHOSTS));

  result->set_attrs("diagnostic", "shelf base mass flux",
                    "kg m-2 s-1", "kg m-2 year-1", "", 0);

  return result;
}

IceModelVec2S::Ptr OceanModel::allocate_water_column_pressure(IceGrid::ConstPtr g) {
  IceModelVec2S::Ptr result(new IceModelVec2S(g,
                                              "average_water_column_pressure",
                                              WITHOUT_GHOSTS));
  result->set_attrs("diagnostic",
                    "vertically-averaged water column pressure",
                    "Pa", "Pa", "", 0);

  return result;
}

// "modifier" constructor
OceanModel::OceanModel(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> input)
  : Component(g), m_input_model(input) {

  if (not input) {
    m_water_column_pressure = allocate_water_column_pressure(g);
  }
}

// "model" constructor
OceanModel::OceanModel(IceGrid::ConstPtr g)
  : OceanModel(g, std::shared_ptr<OceanModel>()) {
  // empty
}

void OceanModel::init(const Geometry &geometry) {
  this->init_impl(geometry);
}

void OceanModel::init_impl(const Geometry &geometry) {
  if (m_input_model) {
    m_input_model->init(geometry);
  } else {
    double
      ice_density   = m_config->get_number("constants.ice.density"),
      water_density = m_config->get_number("constants.sea_water.density"),
      g             = m_config->get_number("constants.standard_gravity");

    compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                             *m_water_column_pressure);
  }
}

void OceanModel::update(const Geometry &geometry, double t, double dt) {
  this->update_impl(geometry, t, dt);
}


const IceModelVec2S& OceanModel::shelf_base_mass_flux() const {
  return shelf_base_mass_flux_impl();
}

const IceModelVec2S& OceanModel::shelf_base_temperature() const {
  return shelf_base_temperature_impl();
}

const IceModelVec2S& OceanModel::average_water_column_pressure() const {
  return average_water_column_pressure_impl();
}

// pass-through default implementations for "modifiers"

void OceanModel::update_impl(const Geometry &geometry, double t, double dt) {
  if (m_input_model) {
    m_input_model->update(geometry, t, dt);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

MaxTimestep OceanModel::max_timestep_impl(double t) const {
  if (m_input_model) {
    return m_input_model->max_timestep(t);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void OceanModel::define_model_state_impl(const File &output) const {
  if (m_input_model) {
    return m_input_model->define_model_state(output);
  } else {
    // no state to define
  }
}

void OceanModel::write_model_state_impl(const File &output) const {
  if (m_input_model) {
    return m_input_model->write_model_state(output);
  } else {
    // no state to write
  }
}

const IceModelVec2S& OceanModel::shelf_base_temperature_impl() const {
  if (m_input_model) {
    return m_input_model->shelf_base_temperature();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

const IceModelVec2S& OceanModel::shelf_base_mass_flux_impl() const {
  if (m_input_model) {
    return m_input_model->shelf_base_mass_flux();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

const IceModelVec2S& OceanModel::average_water_column_pressure_impl() const {
  if (m_input_model) {
    return m_input_model->average_water_column_pressure();
  } else {
    return *m_water_column_pressure;
  }
}

namespace diagnostics {

/*! @brief Shelf base temperature. */
class PO_shelf_base_temperature : public Diag<OceanModel>
{
public:
  PO_shelf_base_temperature(const OceanModel *m)
    : Diag<OceanModel>(m) {

    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "shelfbtemp")};

    set_attrs("ice temperature at the basal surface of ice shelves", "",
              "Kelvin", "Kelvin", 0);
  }
protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "shelfbtemp", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];

    result->copy_from(model->shelf_base_temperature());

    return result;
  }
};



/*! @brief Shelf base mass flux. */
class PO_shelf_base_mass_flux : public Diag<OceanModel>
{
public:
  PO_shelf_base_mass_flux(const OceanModel *m)
    : Diag<OceanModel>(m) {

    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "shelfbmassflux")};

    set_attrs("mass flux at the basal surface of ice shelves", "",
              "kg m-2 s-1", "kg m-2 s-1", 0);
  }
protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "shelfbmassflux", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];

    result->copy_from(model->shelf_base_mass_flux());

    return result;
  }
};

} // end of namespace diagnostics

DiagnosticList OceanModel::diagnostics_impl() const {
  using namespace diagnostics;
  DiagnosticList result = {
    {"shelfbtemp",                     Diagnostic::Ptr(new PO_shelf_base_temperature(this))},
    {"shelfbmassflux",                 Diagnostic::Ptr(new PO_shelf_base_mass_flux(this))}
  };

  if (m_input_model) {
    return combine(m_input_model->diagnostics(), result);
  } else {
    return result;
  }
}

TSDiagnosticList OceanModel::ts_diagnostics_impl() const {
  if (m_input_model) {
    return m_input_model->ts_diagnostics();
  } else {
    return {};
  }
}

void OceanModel::compute_average_water_column_pressure(const Geometry &geometry,
                                                          double ice_density,
                                                          double water_density,
                                                          double g,
                                                          IceModelVec2S &result) {

  auto grid = result.grid();

  const IceModelVec2S
    &bed = geometry.bed_elevation,
    &H   = geometry.ice_thickness,
    &z_s = geometry.sea_level_elevation;

  IceModelVec::AccessList l{&bed, &H, &z_s, &result};

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      result(i, j) = pism::average_water_column_pressure(H(i, j), bed(i, j), z_s(i, j),
                                                         ice_density, water_density, g);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


}


} // end of namespace ocean
} // end of namespace pism
