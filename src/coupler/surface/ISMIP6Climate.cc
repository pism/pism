// Copyright (C) 2019, 2021, 2022, 2023 PISM Authors
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

#include "pism/coupler/surface/ISMIP6Climate.hh"

#include "pism/util/Grid.hh"
#include "pism/coupler/util/options.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace surface {

ISMIP6::ISMIP6(std::shared_ptr<const Grid> grid, std::shared_ptr<atmosphere::AtmosphereModel> input)
  : SurfaceModel(grid),
    m_mass_flux_reference(m_grid, "climatic_mass_balance"),
    m_temperature_reference(m_grid, "ice_surface_temp"),
    m_surface_reference(m_grid, "usurf")
{
  (void) input;

  // allocate model outputs
  m_temperature = allocate_temperature(m_grid);
  m_mass_flux   = allocate_mass_flux(m_grid);

  // set metadata of reference fields
  {
    m_mass_flux_reference.metadata(0)
        .intent("climate_forcing")
        .long_name("reference surface mass balance rate")
        .units("kg m-2 s-1")
        .glaciological_units("kg m-2 year-1")
        .standard_name("land_ice_surface_specific_mass_balance_flux");

    auto smb_max = m_config->get_number("surface.given.smb_max", "kg m-2 second-1");
    m_mass_flux_reference.metadata()["valid_range"] = { -smb_max, smb_max };
    m_mass_flux_reference.set_time_independent(true);

    m_surface_reference.metadata(0)
        .intent("climate_forcing")
        .long_name("reference surface altitude")
        .units("m")
        .standard_name("surface_altitude");

    m_surface_reference.metadata()["valid_range"] = { 0.0, m_grid->Lz() };
    m_surface_reference.set_time_independent(true);

    m_temperature_reference.metadata(0)
        .intent("climate_forcing")
        .long_name("reference temperature")
        .units("Kelvin");

    m_temperature_reference.metadata()["valid_range"] = { 0.0, 373.15 };
    m_temperature_reference.set_time_independent(true);
  }

  // allocate storage for time-dependent inputs
  ForcingOptions opt(*m_grid->ctx(), "surface.ismip6");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    {
      m_mass_flux_anomaly =
          std::make_shared<array::Forcing>(m_grid, file, "climatic_mass_balance_anomaly",
                                           "", // no standard name
                                           buffer_size, opt.periodic);

      m_mass_flux_anomaly->metadata(0)
          .intent("climate_forcing")
          .long_name("surface mass balance rate anomaly")
          .units("kg m-2 s-1")
          .glaciological_units("kg m-2 year-1");
    }

    {
      m_mass_flux_gradient =
          std::make_shared<array::Forcing>(m_grid, file, "climatic_mass_balance_gradient",
                                           "", // no standard name
                                           buffer_size, opt.periodic);

      m_mass_flux_gradient->metadata(0)
          .intent("climate_forcing")
          .long_name("surface mass balance rate elevation lapse rate")
          .units("kg m-2 s-1 m-1")
          .glaciological_units("kg m-2 year-1 m-1");
    }

    {
      m_temperature_anomaly =
          std::make_shared<array::Forcing>(m_grid, file, "ice_surface_temp_anomaly",
                                           "", // no standard name
                                           buffer_size, opt.periodic);

      m_temperature_anomaly->metadata(0)
          .intent("climate_forcing")
          .long_name("ice surface temperature anomaly")
          .units("Kelvin");
    }

    {
      m_temperature_gradient =
          std::make_shared<array::Forcing>(m_grid, file, "ice_surface_temp_gradient",
                                           "", // no standard name
                                           buffer_size, opt.periodic);

      m_temperature_gradient->metadata(0)
          .intent("climate_forcing")
          .long_name("ice surface temperature elevation lapse rate")
          .units("Kelvin m-1");
    }
  }
}

void ISMIP6::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2, "* Initializing the ISMIP6 surface model...\n");

  {
    // File with reference surface elevation, temperature, and climatic mass balance
    auto reference_filename = m_config->get_string("surface.ismip6.reference_file");
    File reference_file(m_grid->com, reference_filename, io::PISM_GUESS, io::PISM_READONLY);

    m_mass_flux_reference.regrid(reference_file, io::CRITICAL);
    m_surface_reference.regrid(reference_file, io::CRITICAL);
    m_temperature_reference.regrid(reference_file, io::CRITICAL);
  }

  {
    ForcingOptions opt(*m_grid->ctx(), "surface.ismip6");

    m_mass_flux_anomaly->init(opt.filename, opt.periodic);
    m_mass_flux_gradient->init(opt.filename, opt.periodic);

    m_temperature_anomaly->init(opt.filename, opt.periodic);
    m_temperature_gradient->init(opt.filename, opt.periodic);
  }
}

void ISMIP6::update_impl(const Geometry &geometry, double t, double dt) {

  // inputs
  const array::Scalar &h       = geometry.ice_surface_elevation;
  const array::Scalar &h_ref   = m_surface_reference;
  const array::Scalar &T_ref   = m_temperature_reference;
  const array::Scalar &SMB_ref = m_mass_flux_reference;

  array::Forcing &dTdz   = *m_temperature_gradient;
  array::Forcing &dSMBdz = *m_mass_flux_gradient;
  array::Forcing &aT     = *m_temperature_anomaly;
  array::Forcing &aSMB   = *m_mass_flux_anomaly;

  // outputs
  array::Scalar &T   = *m_temperature;
  array::Scalar &SMB = *m_mass_flux;

  // get time-dependent input fields at the current time
  {
    aT.update(t, dt);
    aSMB.update(t, dt);
    dTdz.update(t, dt);
    dSMBdz.update(t, dt);

    aT.average(t, dt);
    aSMB.average(t, dt);
    dTdz.average(t, dt);
    dSMBdz.average(t, dt);
  }

  // From http://www.climate-cryosphere.org/wiki/index.php?title=ISMIP6-Projections-Greenland:
  // SMB(x,y,t) = SMB_ref(x,y) + aSMB(x,y,t) + dSMBdz(x,y,t) * [h(x,y,t) - h_ref(x,y)]

  array::AccessScope list{ &h, &h_ref, &SMB, &SMB_ref, &aSMB, &dSMBdz, &T, &T_ref, &aT, &dTdz };

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    SMB(i, j) = SMB_ref(i, j) + aSMB(i, j) + dSMBdz(i, j) * (h(i, j) - h_ref(i, j));
    T(i, j)   = T_ref(i, j) + aT(i, j) + dTdz(i, j) * (h(i, j) - h_ref(i, j));
  }

  dummy_accumulation(SMB, *m_accumulation);
  dummy_melt(SMB, *m_melt);
  dummy_runoff(SMB, *m_runoff);
}

MaxTimestep ISMIP6::max_timestep_impl(double t) const {
  using std::min;

  auto dt = m_temperature_anomaly->max_timestep(t);
  dt      = min(dt, m_temperature_gradient->max_timestep(t));
  dt      = min(dt, m_mass_flux_anomaly->max_timestep(t));
  dt      = min(dt, m_mass_flux_gradient->max_timestep(t));

  if (dt.finite()) {
    return MaxTimestep(dt.value(), "surface ISMIP6");
  }
  return MaxTimestep("surface ISMIP6");
}

const array::Scalar &ISMIP6::mass_flux_impl() const {
  return *m_mass_flux;
}

const array::Scalar &ISMIP6::temperature_impl() const {
  return *m_temperature;
}

const array::Scalar &ISMIP6::accumulation_impl() const {
  return *m_accumulation;
}

const array::Scalar &ISMIP6::melt_impl() const {
  return *m_melt;
}

const array::Scalar &ISMIP6::runoff_impl() const {
  return *m_runoff;
}

} // end of namespace surface
} // end of namespace pism
