// Copyright (C) 2019, 2021, 2022, 2023, 2024, 2025, 2026 PISM Authors
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

#include "pism/coupler/surface/ISMIP7Climate.hh"

#include "pism/util/Grid.hh"
#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/Logger.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace surface {

ISMIP7::ISMIP7(std::shared_ptr<const Grid> grid, std::shared_ptr<atmosphere::AtmosphereModel> input)
  : SurfaceModel(grid),
    m_surface_reference(m_grid, "usurf")
{
  (void) input;

  // allocate model outputs
  m_temperature = allocate_temperature(m_grid);
  m_mass_flux   = allocate_mass_flux(m_grid);
  m_runoff      = allocate_runoff(m_grid);


  // allocate storage for time-dependent inputs
  ForcingOptions opt(*m_grid->ctx(), "surface.ismip7");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);
  // set metadata of reference fields
    
    {
      m_mass_flux_reference =
          std::make_shared<array::Forcing>(m_grid, file, "climatic_mass_balance_gradient",
                                           "", // no standard name
                                           buffer_size, opt.periodic);

      m_mass_flux_reference->metadata(0)
          .long_name("reference surface mass balance rate")
          .units("kg m^-2 s^-1")
          .output_units("kg m^-2 year^-1");
    }

    {
      m_mass_flux_gradient =
          std::make_shared<array::Forcing>(m_grid, file, "climatic_mass_balance_gradient",
                                           "", // no standard name
                                           buffer_size, opt.periodic);

      m_mass_flux_gradient->metadata(0)
          .long_name("surface mass balance rate elevation lapse rate")
          .units("kg m^-2 s^-1 m^-1")
          .output_units("kg m^-2 year^-1 m^-1");
    }

    {
      m_temperature_reference =
          std::make_shared<array::Forcing>(m_grid, file, "reference_temperature",
                                           "", // no standard name
                                           buffer_size, opt.periodic);

      m_temperature_reference->metadata(0)
          .long_name("reference temperature")
          .units("kelvin");
    }

    {
      m_temperature_gradient =
          std::make_shared<array::Forcing>(m_grid, file, "ice_surface_temp_gradient",
                                           "", // no standard name
                                           buffer_size, opt.periodic);

      m_temperature_gradient->metadata(0)
          .long_name("ice surface temperature elevation lapse rate")
          .units("kelvin m^-1");
    }

  }
}

void ISMIP7::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2, "* Initializing the ISMIP7 surface model...\n");

  {
    // File with reference surface elevation, temperature, and climatic mass balance
    auto reference_filename = m_config->get_string("surface.ismip7.reference_file");
    File reference_file(m_grid->com, reference_filename, io::PISM_GUESS, io::PISM_READONLY);

  }

  {
    ForcingOptions opt(*m_grid->ctx(), "surface.ismip7");

    m_mass_flux_gradient->init(opt.filename, opt.periodic);
    
    m_runoff_gradient->init(opt.filename, opt.periodic);

    m_temperature_gradient->init(opt.filename, opt.periodic);
  }
}

void ISMIP7::update_impl(const Geometry &geometry, double t, double dt) {

  // inputs
  const array::Scalar &h       = geometry.ice_surface_elevation;
  const array::Scalar &h_ref   = m_surface_reference;

  array::Forcing &T_ref   = *m_temperature_reference;
  array::Forcing &SMB_ref = *m_mass_flux_reference;
  array::Forcing &runoff_ref = *m_runoff_reference;
  array::Forcing &dTdz   = *m_temperature_gradient;
  array::Forcing &dSMBdz = *m_mass_flux_gradient;
  array::Forcing &drunoffdz = *m_runoff_gradient;

  // outputs
  array::Scalar &T   = *m_temperature;
  array::Scalar &SMB = *m_mass_flux;

  // get time-dependent input fields at the current time
  {
    drunoffdz.update(t, dt);
    dTdz.update(t, dt);
    dSMBdz.update(t, dt);

    drunoffdz.average(t, dt);
    dTdz.average(t, dt);
    dSMBdz.average(t, dt);
  }

  // From http://www.climate-cryosphere.org/wiki/index.php?title=ISMIP7-Projections-Greenland:
  // SMB(x,y,t) = SMB_ref(x,y,t) + dSMBdz(x,y,t) * [h(x,y,t) - h_ref(x,y)]

  array::AccessScope list{ &h, &h_ref, &SMB, &SMB_ref, &dSMBdz, &T, &T_ref, &dTdz, &runoff, &runoff_ref, &drunoffdz };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    SMB(i, j) = SMB_ref(i, j) + dSMBdz(i, j) * (h(i, j) - h_ref(i, j));
    T(i, j)   = T_ref(i, j) + dTdz(i, j) * (h(i, j) - h_ref(i, j));
    runoff(i, j)   = runoff_ref(i, j) + drunoffdz(i, j) * (h(i, j) - h_ref(i, j));
  }

  dummy_accumulation(SMB, *m_accumulation);
  dummy_melt(SMB, *m_melt);
}

MaxTimestep ISMIP7::max_timestep_impl(double t) const {
  using std::min;

  auto dt = m_temperature_reference->max_timestep(t);
  dt      = min(dt, m_temperature_gradient->max_timestep(t));
  dt      = min(dt, m_mass_flux_reference->max_timestep(t));
  dt      = min(dt, m_mass_flux_gradient->max_timestep(t));

  if (dt.finite()) {
    return MaxTimestep(dt.value(), "surface ISMIP7");
  }
  return MaxTimestep("surface ISMIP7");
}

const array::Scalar &ISMIP7::mass_flux_impl() const {
  return *m_mass_flux;
}

const array::Scalar &ISMIP7::temperature_impl() const {
  return *m_temperature;
}

const array::Scalar &ISMIP7::accumulation_impl() const {
  return *m_accumulation;
}

const array::Scalar &ISMIP7::melt_impl() const {
  return *m_melt;
}

const array::Scalar &ISMIP7::runoff_impl() const {
  return *m_runoff;
}

} // end of namespace surface
} // end of namespace pism
