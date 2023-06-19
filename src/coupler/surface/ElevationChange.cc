// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023 PISM Authors
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

#include "ElevationChange.hh"
#include "pism/coupler/util/options.hh"
#include "pism/coupler/util/lapse_rates.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace surface {

ElevationChange::ElevationChange(std::shared_ptr<const Grid> g, std::shared_ptr<SurfaceModel> in)
  : SurfaceModel(g, in) {

  {
    m_smb_lapse_rate = m_config->get_number("surface.elevation_change.smb.lapse_rate",
                                            "(m / s) / m");
    // convert from [m s-1 / m] to [kg m-2 s-1 / m]
    m_smb_lapse_rate *= m_config->get_number("constants.ice.density");

    m_smb_exp_factor = m_config->get_number("surface.elevation_change.smb.exp_factor");
  }

  {
    auto method = m_config->get_string("surface.elevation_change.smb.method");
    m_smb_method = method == "scale" ? SCALE : SHIFT;
  }

  m_temp_lapse_rate = m_config->get_number("surface.elevation_change.temperature_lapse_rate",
                                           "K / m");

  {
    ForcingOptions opt(*m_grid->ctx(), "surface.elevation_change");

    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_reference_surface = std::make_shared<array::Forcing>(m_grid,
                                                      file,
                                                      "usurf",
                                                      "", // no standard name
                                                      buffer_size,
                                                      opt.periodic,
                                                      LINEAR);
    m_reference_surface->set_attrs("climate_forcing", "ice surface elevation",
                                   "m", "m", "surface_altitude", 0);
  }

  m_mass_flux    = allocate_mass_flux(g);
  m_temperature  = allocate_temperature(g);
  m_accumulation = allocate_accumulation(g);
  m_melt         = allocate_melt(g);
  m_runoff       = allocate_runoff(g);
}

void ElevationChange::init_impl(const Geometry &geometry) {
  using units::convert;

  m_input_model->init(geometry);

  m_log->message(2,
                 "  [using temperature and mass balance lapse corrections]\n");

  m_log->message(2,
                 "   ice upper-surface temperature lapse rate: %3.3f K per km\n",
                 convert(m_sys, m_temp_lapse_rate, "K / m", "K / km"));

  if (m_smb_method == SHIFT) {
    double ice_density = m_config->get_number("constants.ice.density");
    m_log->message(2,
                   "   ice-equivalent surface mass balance lapse rate: %3.3f m year-1 per km\n",
                   convert(m_sys, m_smb_lapse_rate, "kg / (m2 second)", "kg / (m2 year)") / ice_density);
  } else {
    m_log->message(2,
                   "   surface mass balance scaling factor with temperature: %3.3f Kelvin-1\n",
                   m_smb_exp_factor);
  }

  ForcingOptions opt(*m_grid->ctx(), "surface.elevation_change");
  m_reference_surface->init(opt.filename, opt.periodic);
}

void ElevationChange::update_impl(const Geometry &geometry, double t, double dt) {

  m_input_model->update(geometry, t, dt);

  m_reference_surface->update(t, dt);
  m_reference_surface->interp(t + 0.5*dt);

  const array::Scalar &surface = geometry.ice_surface_elevation;

  m_temperature->copy_from(m_input_model->temperature());
  lapse_rate_correction(surface, *m_reference_surface,
                        m_temp_lapse_rate, *m_temperature);

  m_mass_flux->copy_from(m_input_model->mass_flux());

  switch (m_smb_method) {
  case SCALE:
    {
      array::AccessScope list{&surface, m_reference_surface.get(), m_mass_flux.get()};

      for (auto p = m_grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();

        double dT = -m_temp_lapse_rate * (surface(i, j) - (*m_reference_surface)(i, j));

        (*m_mass_flux)(i, j) *= exp(m_smb_exp_factor * dT);
      }
    }
    break;
  default:
  case SHIFT:
    {
      lapse_rate_correction(surface, *m_reference_surface,
                            m_smb_lapse_rate, *m_mass_flux);
    }
    break;
  }

  // This modifier changes m_mass_flux, so we need to compute accumulation, melt, and
  // runoff.
  dummy_accumulation(*m_mass_flux, *m_accumulation);
  dummy_melt(*m_mass_flux, *m_melt);
  dummy_runoff(*m_mass_flux, *m_runoff);

}

const array::Scalar &ElevationChange::mass_flux_impl() const {
  return *m_mass_flux;
}

const array::Scalar &ElevationChange::temperature_impl() const {
  return *m_temperature;
}

const array::Scalar &ElevationChange::accumulation_impl() const {
  return *m_accumulation;
}

const array::Scalar &ElevationChange::melt_impl() const {
  return *m_melt;
}

const array::Scalar &ElevationChange::runoff_impl() const {
  return *m_runoff;
}


} // end of namespace surface
} // end of namespace pism
