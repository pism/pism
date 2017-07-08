// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include <gsl/gsl_math.h>

#include "PSLapseRates.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace surface {

LapseRates::LapseRates(IceGrid::ConstPtr g, SurfaceModel* in)
  : PLapseRates<SurfaceModel,SurfaceModifier>(g, in) {
  m_smb_lapse_rate = 0;
  m_option_prefix = "-surface_lapse_rate";
}

LapseRates::~LapseRates() {
  // empty
}

void LapseRates::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_input_model->init();

  m_log->message(2,
             "  [using temperature and mass balance lapse corrections]\n");

  init_internal();

  m_smb_lapse_rate = options::Real("-smb_lapse_rate",
                                   "Elevation lapse rate for the surface mass balance,"
                                   " in m year-1 per km",
                                   m_smb_lapse_rate);

  m_log->message(2,
             "   ice upper-surface temperature lapse rate: %3.3f K per km\n"
             "   ice-equivalent surface mass balance lapse rate: %3.3f m year-1 per km\n",
             m_temp_lapse_rate, m_smb_lapse_rate);

  m_temp_lapse_rate = units::convert(m_sys, m_temp_lapse_rate, "K/km", "K/m");

  // convert from [m year-1 / km] to [kg m-2 year-1 / km]
  m_smb_lapse_rate *= m_config->get_double("constants.ice.density");
  m_smb_lapse_rate = units::convert(m_sys, m_smb_lapse_rate,
                                    "(kg m-2) year-1 / km", "(kg m-2) second-1 / m");

}

void LapseRates::mass_flux_impl(IceModelVec2S &result) const {
  m_input_model->mass_flux(result);
  lapse_rate_correction(result, m_smb_lapse_rate);
}

void LapseRates::temperature_impl(IceModelVec2S &result) const {
  m_input_model->temperature(result);
  lapse_rate_correction(result, m_temp_lapse_rate);
}

SpatialSMBGradients::SpatialSMBGradients(IceGrid::ConstPtr g, SurfaceModel* in)
  : PLapseRates<SurfaceModel,SurfaceModifier>(g, in) {
  m_smb_lapse_rate = 0;
  m_option_prefix = "-spatial_smb_gradients";
}

SpatialSMBGradients::~SpatialSMBGradients() {
  // empty
}

void SpatialSMBGradients::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_input_model->init();

  m_log->message(2,
             "  [using temperature and mass balance lapse corrections]\n");

  init_internal();

  m_log->message(2,
             "   ice upper-surface temperature lapse rate: %3.3f K per km\n",
             m_temp_lapse_rate);

  m_temp_lapse_rate = units::convert(m_sys, m_temp_lapse_rate, "K/km", "K/m");

  options::RealList option("-spatial_smb_gradients",
                                   "Elevation lapse rate for the surface mass balance,"
                                   " in kg / (m2 yr)");
  double
    smb_grad_negative_south  = 0.0,
    smb_grad_positive_south  = 0.0,
    smb_grad_negative_north  = 0.0,
    smb_grad_positive_north  = 0.0,
    north_south_lat = 0;
  
  if (not option.is_set()) {
    smb_grad_negative_south  = m_config->get_double("surface.lapse_rate.smb_grad_negative_south");
    smb_grad_positive_south  = m_config->get_double("surface.lapse_rate.smb_grad_positive_south");
    smb_grad_negative_north  = m_config->get_double("surface.lapse_rate.smb_grad_negative_north");
    smb_grad_positive_north  = m_config->get_double("surface.lapse_rate.smb_grad_positive_north");
    north_south_lat          = m_config->get_double("surface.lapse_rate.north_south_lat");
  } else {
    smb_grad_negative_south  = option[0];
    smb_grad_positive_south  = option[1];
    smb_grad_negative_north  = option[2];
    smb_grad_positive_north  = option[3];
    north_south_lat          = option[4];

  }

  if (option.is_set() && option->size() != 5) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid -smb_lapse_rate_spatial arguments: has to be a list"
                                  " of 5 numbers, got %d", (int)option->size());
  }

  m_log->message(2,
                 "  SMB gradients after Edwards et al :\n"
                 "  SMB<0 north of %.f: %.f\n"
                 "  SMB>0 north of %.f: %.f\n"
                 "  SMB<0 south of %.f: %.f\n"
                 "  SMB>0 south of %.f: %.f\n",
                 smb_grad_negative_north, north_south_lat,
                 smb_grad_positive_north, north_south_lat,
                 smb_grad_negative_south, north_south_lat,
                 smb_grad_positive_south, north_south_lat);

  m_log->message(2,
                 "   ice upper-surface temperature lapse rate: %3.3f K per km\n",
                 m_temp_lapse_rate);

  m_temp_lapse_rate = units::convert(m_sys, m_temp_lapse_rate, "K/km", "K/m");


  const IceModelVec2S
    &climatic_mass_balance = *m_grid->variables().get_2d_scalar("latitude"),
    &latitude = *m_grid->variables().get_2d_scalar("latitude");

  IceModelVec::AccessList list{&climatic_mass_balance, &latitude};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (latitude(i, j) < north_south_lat) {
      if (climatic_mass_balance(i, j) < 0.) {
        m_smb_lapse_rate = smb_grad_negative_south;
      } else {
        m_smb_lapse_rate = smb_grad_positive_south;
      }
    } else {
      if (climatic_mass_balance(i, j) > 0.) {
        m_smb_lapse_rate = smb_grad_negative_north;
      } else {
        m_smb_lapse_rate = smb_grad_positive_north;
      }
    }
  }

  m_smb_lapse_rate = units::convert(m_sys, m_smb_lapse_rate,
                                    "(kg m-2) year-1 / km", "(kg m-2) second-1 / m");

}

void SpatialSMBGradients::mass_flux_impl(IceModelVec2S &result) const {
  m_input_model->mass_flux(result);
  lapse_rate_correction(result, m_smb_lapse_rate);
}

void SpatialSMBGradients::temperature_impl(IceModelVec2S &result) const {
  m_input_model->temperature(result);
  lapse_rate_correction(result, m_temp_lapse_rate);
}

} // end of namespace surface
} // end of namespace pism
