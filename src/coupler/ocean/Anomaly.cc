// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021, 2022, 2023 PISM Authors
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

#include "pism/coupler/ocean/Anomaly.hh"
#include "pism/util/Grid.hh"
#include "pism/coupler/util/options.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace ocean {

Anomaly::Anomaly(std::shared_ptr<const Grid> g, std::shared_ptr<OceanModel> in)
  : OceanModel(g, in) {

  ForcingOptions opt(*m_grid->ctx(), "ocean.anomaly");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_shelf_base_mass_flux_anomaly = std::make_shared<array::Forcing>(m_grid,
                                                                  file,
                                                                  "shelf_base_mass_flux_anomaly",
                                                                  "", // no standard name
                                                                  buffer_size,
                                                                  opt.periodic);
  }

  m_shelf_base_mass_flux_anomaly->metadata(0)
      .intent("climate_forcing")
      .long_name("anomaly of the shelf base mass flux rate")
      .units("kg m-2 s-1")
      .glaciological_units("kg m-2 year-1");

  m_shelf_base_mass_flux = allocate_shelf_base_mass_flux(g);


}

void Anomaly::init_impl(const Geometry &geometry) {

  if (m_input_model) {
    m_input_model->init(geometry);
  }

  m_log->message(2,
                 "* Initializing the '-ocean ...,anomaly' modifier...\n");

  ForcingOptions opt(*m_grid->ctx(), "ocean.anomaly");

  m_log->message(2,
                 "    reading anomalies from %s ...\n", opt.filename.c_str());

  m_shelf_base_mass_flux_anomaly->init(opt.filename, opt.periodic);
}

void Anomaly::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_shelf_base_mass_flux_anomaly->update(t, dt);

  m_shelf_base_mass_flux_anomaly->average(t, dt);

  m_input_model->shelf_base_mass_flux().add(1.0, *m_shelf_base_mass_flux_anomaly,
                                 *m_shelf_base_mass_flux);
}

const array::Scalar &Anomaly::shelf_base_mass_flux_impl() const {
  return *m_shelf_base_mass_flux;
}


} // end of namespace ocean
} // end of namespace pism
