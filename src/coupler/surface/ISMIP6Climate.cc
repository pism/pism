// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019 PISM Authors
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

#include "ISMIP6Climate.hh"

#include "pism/util/Time.hh"
#include "pism/util/IceGrid.hh"
#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace surface {

ISMIP6::ISMIP6(IceGrid::ConstPtr grid, std::shared_ptr<atmosphere::AtmosphereModel> input)
  : SurfaceModel(grid)
{
  (void) input;

  ForcingOptions opt(*m_grid->ctx(), "surface.ismip6");

  {
    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    
    // We need two files; the reference file with reference SMB and
    // elevation (time-independent) and
    // the anomaly/gradient file (time-dependent)
    std::string reference_filename = m_config->get_string("surface.ismip6.reference_file");
    
    // File with anomalies and gradients of temperature and climatic mass balance
    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

    m_mass_flux_anomaly = IceModelVec2T::ForcingField(m_grid,
                                              file,
                                              "aSMB",
                                              "", // no standard name
                                              buffer_size,
                                              evaluations_per_year,
                                              periodic);

    m_mass_flux_gradient = IceModelVec2T::ForcingField(m_grid,
                                              file,
                                              "dSMBdz",
                                              "", // no standard name
                                              buffer_size,
                                              evaluations_per_year,
                                              periodic);

    m_mass_flux_reference->create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);

    m_surface_reference->create(m_grid, "surface_elevation", WITHOUT_GHOSTS);

    m_temperature_anomaly = IceModelVec2T::ForcingField(m_grid,
                                                file,
                                                "ice_surface_temp_anomaly",
                                                "", // no standard name
                                                buffer_size,
                                                evaluations_per_year,
                                                periodic);

    m_temperature_gradient = IceModelVec2T::ForcingField(m_grid,
                                                file,
                                                "ice_surface_temp_gradient",
                                                "", // no standard name
                                                buffer_size,
                                                evaluations_per_year,
                                                periodic);

    m_temperature_reference->create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);

  }

  const double smb_max = m_config->get_double("surface.given.smb_max", "kg m-2 second-1");


  m_mass_flux_anomaly->set_attrs("climate_forcing",
                         "surface mass balance (accumulation/ablation) rate",
                                 "kg m-2 s-1", "kg m-2 year-1",
                                 "land_ice_surface_specific_mass_balance_flux", 0);
  
  m_mass_flux_gradient->set_attrs("climate_forcing",
                                  "surface mass balance (accumulation/ablation) rate",
                                  "kg m-2 s-1 m-1", "kg m-2 year-1 m-1",
                                  "land_ice_surface_specific_mass_balance_flux", 0);
  
  m_mass_flux_reference->set_attrs("climate_forcing",
                                   "surface mass balance (accumulation/ablation) rate",
                                   "kg m-2 s-1", "kg m-2 year-1",
                                   "land_ice_surface_specific_mass_balance_flux", 0);

  m_mass_flux_reference->metadata().set_doubles("valid_range", {-smb_max, smb_max});
  m_mass_flux_reference->set_time_independent(true);

  m_surface_reference->set_attrs("surface_altitude",
                                 "reference surface altitude",
                                 "m", "m", "surface_altitude", 0);
  m_surface_reference->set_time_independent(true);

  m_temperature->set_attrs("climate_forcing",
                           "temperature of the ice at the ice surface but below firn processes",
                           "Kelvin", "Kelvin", "", 0);
  m_temperature->metadata().set_doubles("valid_range", {0.0, 323.15}); // [0C, 50C]


  m_temperature_anomaly->set_attrs("climate_forcing_anomaly",
                                   "temperature of the ice at the ice surface but below firn processes",
                                   "Kelvin", "Kelvin", "", 0);

  m_temperature_gradient->set_attrs("climate_forcing_gradient",
                                    "temperature of the ice at the ice surface but below firn processes",
                                    "Kelvin m-1", "Kelvin m-1", "", 0);
}

ISMIP6::~ISMIP6() {
  // empty
}

void ISMIP6::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing the ISMIP6 surface model\n");
  ForcingOptions opt(*m_grid->ctx(), "surface.ismip6");

  // File with reference surface elevation, temperature, and climatic mass balance
  PIO reference_file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

  m_mass_flux_reference->regrid(reference_file, CRITICAL);
  m_surface_reference->regrid(reference_file, CRITICAL);
  m_temperature_reference->regrid(reference_file, CRITICAL);
  
  m_mass_flux_anomaly->init(opt.filename, opt.period, opt.reference_time);
  m_mass_flux_gradient->init(opt.filename, opt.period, opt.reference_time);

  m_temperature_anomaly->init(opt.filename, opt.period, opt.reference_time);
  m_temperature_gradient->init(opt.filename, opt.period, opt.reference_time);
}

void ISMIP6::update_impl(const Geometry &geometry, double t, double dt) {

  dummy_accumulation(*m_mass_flux, *m_accumulation);
  dummy_melt(*m_mass_flux, *m_melt);
  dummy_runoff(*m_mass_flux, *m_runoff);
  
  // From http://www.climate-cryosphere.org/wiki/index.php?title=ISMIP6-Projections-Greenland:
  // SMB(x,y,t) = SMB_ref(x,y) + aSMB(x,y,t) + dSMBdz(x,y) * [h(x,y,t) - h_ref(x,y)]
  
  const IceModelVec2S &surface = geometry.ice_surface_elevation;

  IceModelVec::AccessList list
    {&surface, m_mass_flux_reference.get(), m_mass_flux_anomaly.get(), m_mass_flux_gradient.get(), m_surface_reference.get(),  m_temperature_anomaly.get(), m_temperature_gradient.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    
    (*m_mass_flux)(i, j) = (*m_mass_flux_reference)(i, j) + (*m_mass_flux_anomaly)(i, j) + (*m_mass_flux_gradient)(i, j) * (surface(i, j) - (*m_surface_reference)(i, j));

    (*m_temperature)(i, j) = (*m_temperature_reference)(i, j) + (*m_temperature_anomaly)(i, j) + (*m_temperature_gradient)(i, j) * (surface(i, j) - (*m_surface_reference)(i, j));
  }

}

const IceModelVec2S &ISMIP6::mass_flux_impl() const {
  
  return *m_mass_flux;
}

const IceModelVec2S &ISMIP6::temperature_impl() const {
  return *m_temperature;
}

const IceModelVec2S &ISMIP6::accumulation_impl() const {
  return *m_accumulation;
}

const IceModelVec2S &ISMIP6::melt_impl() const {
  return *m_melt;
}

const IceModelVec2S &ISMIP6::runoff_impl() const {
  return *m_runoff;
}

void ISMIP6::define_model_state_impl(const PIO &output) const {
  m_mass_flux->define(output);
  m_temperature->define(output);
}

void ISMIP6::write_model_state_impl(const PIO &output) const {
  m_mass_flux->write(output);
  m_temperature->write(output);
}

} // end of namespace surface
} // end of namespace pism
