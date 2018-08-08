// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "OrographicPrecipitation.hh"

#include "pism/coupler/util/options.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/geometry/Geometry.hh"
#include "OrographicPrecipitationSerial.hh"

namespace pism {
namespace atmosphere {

OrographicPrecipitation::OrographicPrecipitation(IceGrid::ConstPtr g)
  : AtmosphereModel(g, nullptr), m_surface(g, "ice_surface_elevation", WITHOUT_GHOSTS)
{
  ForcingOptions opt(*m_grid->ctx(), "atmosphere.orographic_precipitation");
  
  {
    
    m_snow_temp_july_day = m_config->get_double("atmosphere.fausto_air_temp.summer_peak_day");
    
    // Allocate internal IceModelVecs:
    m_air_temp_mean_annual.create(m_grid, "air_temp_mean_annual", WITHOUT_GHOSTS);
    m_air_temp_mean_annual.set_attrs("diagnostic",
                                     "mean annual near-surface air temperature (without sub-year time-dependence or forcing)",
                                     "K",
                                     "");  // no CF standard_name ??
    
    m_air_temp_mean_july.create(m_grid, "air_temp_mean_july", WITHOUT_GHOSTS);
    m_air_temp_mean_july.set_attrs("diagnostic",
                                   "mean July near-surface air temperature (without sub-year time-dependence or forcing)",
                                   "Kelvin",
                                   "");  // no CF standard_name ??
    
    m_precipitation.create(m_grid, "precipitation", WITHOUT_GHOSTS);
    m_precipitation.set_attrs("model_state", "precipitation rate",
                              "kg m-2 second-1", "precipitation_flux", 0);
    m_precipitation.metadata(0).set_string("glaciological_units", "kg m-2 year-1");
    m_precipitation.set_time_independent(true);
    
  }

  m_work0 = m_surface.allocate_proc0_copy();

  const int
    Mx = m_grid->Mx(),
    My = m_grid->My(),
    Z  = m_config->get_double("bed_deformation.lc.grid_size_factor"),
    Nx = Z*(Mx - 1) + 1,
    Ny = Z*(My - 1) + 1;

  const double
    Lx = Z * (m_grid->x0() - m_grid->x(0)),
    Ly = Z * (m_grid->y0() - m_grid->y(0));

  m_extended_grid = IceGrid::Shallow(m_grid->ctx(),
                                     Lx, Ly,
                                     m_grid->x0(), m_grid->y0(),
                                     Nx, Ny, CELL_CORNER, NOT_PERIODIC);
  
  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      m_serial_model.reset(new OrographicPrecipitationSerial(*m_config,
                                           Mx, My,
                                           m_grid->dx(), m_grid->dy(),
                                           Nx, Ny));
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}

OrographicPrecipitation::~OrographicPrecipitation() {
  // empty
}

//! Copies the stored mean annual near-surface air temperature field into result.
const IceModelVec2S& OrographicPrecipitation::mean_annual_temp_impl() const {
  return m_air_temp_mean_annual;
}

//! Copies the stored mean July near-surface air temperature field into result.
const IceModelVec2S& OrographicPrecipitation::mean_july_temp() const {
  return m_air_temp_mean_july;
}


void OrographicPrecipitation::init_impl(const Geometry &geometry) {
  (void) geometry;
  
  m_log->message(2,
             "* Initializing the atmosphere model computing precipitation using the\n"
             "  Linear Theory of Orographic Precipitation model with scalar wind speeds\n"
             "  and reading near-surface air temperature from a file...\n");

  m_reference =
    "R. B. Smith and I. Barstad, 2004.\n"
    "A Linear Theory of Orographic Precipitation. J. Atmos. Sci. 61, 1377-1391.";

  m_precipitation.metadata().set_string("source", m_reference);

}

void OrographicPrecipitation::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  m_log->message(2,
                 "* UPDATING\n");


  // make a copy of the surface elevation so that it is available in methods computing
  // precipitation
  m_surface.copy_from(geometry.ice_surface_elevation);
  m_surface.put_on_proc0(*m_work0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {  // only processor zero does the step
      PetscErrorCode ierr = 0;

      m_serial_model->step(*m_work0);

      ierr = VecCopy(m_serial_model->orographic_precipitation(), *m_work0);
      PISM_CHK(ierr, "VecCopy");
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  m_precipitation.get_from_proc0(*m_work0);

}

void OrographicPrecipitation::init_timeseries_impl(const std::vector<double> &ts) const {

  // constants related to the standard yearly cycle
  const double
    julyday_fraction = m_grid->ctx()->time()->day_of_the_year_to_day_fraction(m_snow_temp_july_day);

  size_t N = ts.size();

  m_ts_times.resize(N);
  m_cosine_cycle.resize(N);
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    double tk = m_grid->ctx()->time()->year_fraction(ts[k]) - julyday_fraction;

    m_ts_times[k] = ts[k];
    m_cosine_cycle[k] = cos(2.0 * M_PI * tk);
  }
  
}

void OrographicPrecipitation::temp_time_series_impl(int i, int j, std::vector<double> &result) const {

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] = m_air_temp_mean_annual(i,j) + (m_air_temp_mean_july(i,j) - m_air_temp_mean_annual(i,j)) * m_cosine_cycle[k];
  }
}

void OrographicPrecipitation::begin_pointwise_access_impl() const {
  m_precipitation.begin_access();
  m_air_temp_mean_annual.begin_access();
  m_air_temp_mean_july.begin_access();
  m_surface.begin_access();
}

void OrographicPrecipitation::end_pointwise_access_impl() const {
  m_precipitation.end_access();
  m_air_temp_mean_annual.end_access();
  m_air_temp_mean_july.end_access();
  m_surface.end_access();
}

DiagnosticList OrographicPrecipitation::diagnostics_impl() const {
  DiagnosticList result = AtmosphereModel::diagnostics_impl();

  result["air_temp_mean_july"] = Diagnostic::Ptr(new PA_mean_july_temp_op(this));

  return result;
}

PA_mean_july_temp_op::PA_mean_july_temp_op(const OrographicPrecipitation *m)
  : Diag<OrographicPrecipitation>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "air_temp_mean_july")};

  set_attrs("mean July near-surface air temperature used in the cosine yearly cycle", "",
            "Kelvin", "Kelvin", 0);
}

IceModelVec::Ptr PA_mean_july_temp_op::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "air_temp_mean_july", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->mean_july_temp());

  return result;
}

} // end of namespace atmosphere
} // end of namespace pism
